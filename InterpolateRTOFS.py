# This is a a set of routines that interpolate wind forecasts to an unstructured 
# mesh from regular or curvilinear grids (using ESMPY). In addition to interpolation 
# field the out file contains a spatially variable error variance estimate based on the 
# distance to the forecast boundary. This error variance is used later to update 
# dispirate forecasts in a bayesian manner to yield a spatially smooth estimate 
# incorporating all forecasts. The understanding is that the more local forecast 
# products are of higher accuracy than the coarser broader scale forecasts

import numpy as np
import os

import datetime
import netCDF4 as nc
import sys
import re
import InterpNWPSUtility as nwps

import xarray as xr
import esmpy
import scipy.sparse as sp

import GeoComputeMeshToMeshInterpWeights as mshint
import InterpNWPSUtility as  mshutil
def QuickDistance(lat1, lon1, lats2, lons2):
    deg2kmY=111.
    deg2kmX=np.cos( np.pi * lat1 / 180.)*deg2kmY
    d= np.min(  np.sqrt( (  (lat1-lats2)*deg2kmY)**2 + ((lon1-lons2)*deg2kmX)**2 )  )
    return d

def VarianceLinearDistanceToBndy(InteriorNodeList, DistanceToBoundary, InteriorVariance, VarianceOnBoundary, LengthScale):
    Variance=np.zeros(len(DistanceToBoundary))+np.inf
    SpatialFunction=DistanceToBoundary/LengthScale
    j=np.where(SpatialFunction>1.)
    SpatialFunction[j]=1.
    Variance[InteriorNodeList] = VarianceOnBoundary + ( InteriorVariance - VarianceOnBoundary ) * SpatialFunction[InteriorNodeList]
    return Variance

def VarianceInverseDistanceToBndy(InteriorNodeList, DistanceToBoundary, InteriorVariance, LengthScale):
    Variance=np.zeros(len(DistanceToBoundary))+np.inf
    SpatialFunction=LengthScale / DistanceToBoundary
    j=np.where(SpatialFunction>1.)
    SpatialFunction[j]=1.
    Variance[InteriorNodeList] = InteriorVariance  * SpatialFunction[InteriorNodeList]
    return Variance

def CurvilinearGridCreateInterpWeights(xi,yi,x1,y1, weights_file):
# Compute interpolation weights to interpolate from curvilinear grid (x1,y1) to points (xi,yi)
# and store in netcdf file using ESMPY
    nx,ny=x1.shape
    nn=len(xi)
    n1=nx*ny
    src_lon = x1
    src_lat = y1
    dst_lon=np.zeros((1,nn))
    dst_lat=np.zeros((1,nn))
    print("nn="+str(nn))
    dst_lon[0,:] = xi[:]
    dst_lat[0,:] = yi[:]

    esmpy.Manager()

    src_grid = esmpy.Grid(
      max_index=np.array([src_lon.shape[0], src_lon.shape[1]]),
      staggerloc=esmpy.StaggerLoc.CENTER,
      coord_sys=esmpy.CoordSys.SPH_DEG
    )

    dst_grid = esmpy.Grid(
      max_index=np.array([nn, 1]),
      staggerloc=esmpy.StaggerLoc.CENTER,
      coord_sys=esmpy.CoordSys.SPH_DEG
    )

    # 4. Populate the coordinate data
    src_lon_ptr = src_grid.get_coords(0)
    src_lat_ptr = src_grid.get_coords(1)
    src_lon_ptr[...] = src_lon
    src_lat_ptr[...] = src_lat

    dst_lon_ptr = dst_grid.get_coords(0)
    dst_lat_ptr = dst_grid.get_coords(1)
    dst_lon_ptr[...] = dst_lon.T
    dst_lat_ptr[...] = dst_lat.T

    # 5. Create esmpy Fields
    src_field = esmpy.Field(src_grid, name="src_field")
    dst_field = esmpy.Field(dst_grid, name="dst_field")
    src_field.data[...]=np.sqrt(x1/180)/(90+y1) # arbitrary function of x,y

    np.savetxt('F.txt', src_field.data[...])
    np.savetxt('X.txt', x1)
    np.savetxt('Y.txt', y1)
    print("Creating weights: "+ weights_file)
    regrid = esmpy.Regrid(
      src_field,
      dst_field,
      filename=weights_file,
      regrid_method=esmpy.RegridMethod.BILINEAR,
      ignore_degenerate=True, # <--- Add this parameter
      unmapped_action=esmpy.UnmappedAction.IGNORE # Optional: Ignores missing/masked points
    )
#Add number of rows and columns to weights file for clarity when constructing sparse matrix for interpolation
    with nc.Dataset(weights_file, mode="a") as ds:
        ds.Nrows = nn
        ds.Ncols = n1
    
    np.savetxt('Fi.txt', dst_field.data[...])
    np.savetxt('xi.txt', xi)
    np.savetxt('yi.txt', yi)
    
    return

def CalculateDistanceToBoundary(xi,yi,x1,y1):
# Distance to boundary calculation for use when interpolation envelope corresponds with 
# interior of curvilinear grid boundary.
#
# Inputs:
#   xi (nn): longitude of unstructured mesh nodes 
#   yi (nn): latitude of unstructured mesh nodes 
#   x1 (nx x ny): longitude for interpoltated field
#   y1 (nx x ny): latitude for interpoltated field
#
# Outputs:
#   dist2bnd (nn) : distance to edge of interpolation envelope.  
    nx=x1.shape[0]
    ny=x1.shape[1]
    xb=np.hstack((x1[1,:],x1[:,ny-1].T,x1[nx-1,:],x1[:,1].T))
    yb=np.hstack((y1[1,:],y1[:,ny-1].T,y1[nx-1,:],y1[:,1].T))
    np.savetxt('xbyb.txt', np.vstack((xb,yb)))
    dist2bnd=np.zeros(nn)
    for k in range(nn):
        dist2bnd[k]=QuickDistance(yi[k],xi[k],yb,xb)
        if k%10000==0:
            print("calculating distance to boundary, "+str(k)+":"+ str(nn)+":"+str(k/nn) )
    return dist2bnd

def CalculateDistanceToInterpEnvelope(xi,yi,fi,SearchWidth):
# Alternative distance to boundary calculation for use when interpolation envelope is 
# distinctly interior to curvilinear grid boundary as happens for RRFS NA grid 
#
# Inputs:
#   xi (nn): longitude of unstructured mesh nodes 
#   yi (nn): latitude of unstructured mesh nodes 
#   fi (nn): interpolated field on mesh nodes with 'nan' values outside of interpolation envelope
#   SearchWidth: Computational speed up to remove extra search points in distance to boundary
#
# Outputs:
#   dist2bnd (nn) : distance to edge of interpolation envelope.  dist2bnd[k]=0 if (xi[k],yi[k]) is outside
#                   of the interpolation envelope
    dist2bnd=np.zeros(nn)
    jin = np.where(~np.isnan(fi))[0].tolist()#points inside interpolation envelope
    jout = np.where(np.isnan(fi))[0].tolist() #points outside interpolation envelope
    xin=xi[jin]
    yin=yi[jin]
    xout=xi[jout]
    yout=yi[jout]
    jxU=np.where( xout < np.max(xin)+SearchWidth )[0].tolist()
    jxD=np.where( xout > np.min(xin)-SearchWidth )[0].tolist()
    jyU=np.where( yout < np.max(yin)+SearchWidth )[0].tolist()
    jyD=np.where( yout > np.min(yin)-SearchWidth )[0].tolist()
    j=list( set(jxU) & set(jxD)  & set(jyU)   & set(jyD)  )
    xout=xout[j]
    yout=yout[j]
    din=np.zeros(len(jin))
    print(len(xin))
    print(len(xout))
    for k in range(len(xin)):
        din[k]=QuickDistance(yin[k],xin[k],yout,xout)
        if k%10000==0:
            print("calculating distance to boundary, "+str(k)+":"+ str(nn)+":"+str(k/nn) )
    dist2bnd[jin]=din
    return dist2bnd


# Main program

nargin = len(sys.argv) - 1

if nargin < 5 :
    print("missing argument for specifying spatial forecast variance")
    print("call as: InterpWindToMesh.crvln.py flin mesh flout variance lambda [Extrapolate-optional]")

flin=sys.argv[1]
mshfl=sys.argv[2]
flout=sys.argv[3]
ForecastErrorVariance=float(sys.argv[4]) #(m m /s /s)
LambdaBndyErrorVariance=float(sys.argv[5]) #(km)

# Don't use nearest neighbor interpolation unless 6th positive integer argument present.
# This may be needed on boundary of RWPS mesh if node alignment is outside NBM OC domain
Extrapolate=False
if nargin >5:
    if int(sys.argv[6])>0:
        print("using nearest neighbor to extrapolate wind field beyond coverage")
        Extrapolate=True

xi, yi, ei, zi = nwps.loadWW3Mesh(mshfl)
nn=len(xi)

if np.mean(xi)<0:
	xi=xi+360.

data = nc.Dataset(flin,"r")
#read spaital dimensions and determine if input mesh is curvilinear or regular
x1=np.asarray(data["Longitude"][:])
y1=np.asarray(data["Latitude"][:])
# last RTOFS row is degenerate
x1=x1[1:-1,1:-1]
y1=y1[1:-1,1:-1]
print(y1.shape)
if (len(x1.shape)==2 and len(y1.shape)==2):
    IsCrvLn=True
elif (len(x1.shape)==1 and len(y1.shape)==1):
    IsCrvLn=False
else:
    print("input file spatial dimension is not right. ending program")
    sys.exit()
# shift to common longitude
if np.mean(np.mean(x1))<0:
    x1=x1+360.

#######################################
# === Create weights or resuse it  ===#
#######################################

meshslash=mshfl.rfind('/')+1

flinslash=flin.rfind('/')+1
flint=flin[flinslash : ]
flins=flint.split("_")
model=flins[0]
dom=flins[1]
weights_file = "InterpolationWeights."+mshfl[meshslash:-3]+model+"."+dom+".nc"
print("interpolation weights file = "+ weights_file)

if IsCrvLn:
    nx=x1.shape[0]
    ny=x1.shape[1]
else:
    nx=len(x1)
    ny=len(y1)

if not os.path.isfile(weights_file):
    print("No existing interpolation weights. Computing weights and saving to file: "+ weights_file)
    if not IsCrvLn:
        x1 = np.tile(x1,(ny,1))
        y1 = np.tile(y1,(nx,1)).T
    CurvilinearGridCreateInterpWeights(xi, yi, x1, y1, weights_file)

n1=nx*ny

##################################################################################
# START: Interpolate wind data using predefined weights
##################################################################################

#######################################
# === resuse weights  ===#
#######################################
with xr.open_dataset(weights_file) as ds_s:
   # Standard sparse storage uses 'row', 'col', and 'data' variables
   row = ds_s['row'].values
   col = ds_s['col'].values
   weights = ds_s['S'].values
   Nrows=ds_s.attrs.get('Nrows')
   Ncols=ds_s.attrs.get('Ncols')
   
print("nn = "+str(nn)+": Nrows = "+str(Nrows))
print("n1 = "+str(n1)+": Ncols = "+str(Ncols))
if not ((nn==Nrows) and (n1==Ncols)):
    print("Wrong matrix weights: number of rows from "+ mshfl +" = "+str(nn)+
    " but number of rows in "+ weights_file +" = "+str(Nrows)+ 
    ", number of spatial points in "+ flin +" = "+str(n1)+ 
    " but number of columns in "+ weights_file +" = "+str(Ncols)  )
    print("You probably need to remove file "+ weights_file +" and rerun to generate appropriate weights")

matrix = sp.coo_matrix((weights, (row-1, col-1)), shape=(nn,n1)).tocsr()
print("sparse interpolation matrix")
print(matrix)
######################################
t1=mshutil.ConvertTimeToUnixTime(flin,"MT")
#t1=np.asarray(data["MT"][:])
nt=len(t1)
nan=float("nan")

# read and reshape input data for interpolation
Up=np.zeros((nt,n1))+nan
Vp=np.zeros((nt,n1))+nan
fill_value0=data["u_velocity"]._FillValue
print("fill value="+str(fill_value0))
for k in range(nt):
    U=data["u_velocity"][k,0,1:-1,1:-1]
    Up[k,:]=np.transpose(U).reshape(n1)
    V=data["v_velocity"][k,0,1:-1,1:-1]
    Vp[k,:]=np.transpose(V).reshape(n1)

j=np.where( Up==fill_value0 )
Up[j]=nan
Vp[j]=nan
j=np.where( Vp==fill_value0 )
Up[j]=nan
Vp[j]=nan

#Carry out interpolation with sparse matrix multiplication
u = matrix @ Up.T
v = matrix @ Vp.T

#Use fill value where the interpolator has no coverage
row_sum = matrix.sum(axis=1)
j0=np.where(row_sum==0)
j0=j0[0].tolist()
#np.savetxt('rowsum.txt', row_sum)
#np.savetxt('j0.txt', j0)
u[j0,:]=nan
v[j0,:]=nan
##################################################################################
# FINISHED: Interpolate wind data using predefined weights
##################################################################################

##################################################################################
# START: Compute distance to boundary for each node in mesh:
##################################################################################
dist2bnd_file = "DistToBndy."+mshfl[meshslash:-3]+model+"."+dom+".nc"
if os.path.isfile(dist2bnd_file):
    print(f"Reusing distance to file: {dist2bnd_file}")
    #dist2bnd = np.loadtxt(dist2bnd_file)
    with xr.open_dataset(dist2bnd_file) as ds_s:
        dist2bnd = ds_s['dist2bnd'].values
else:
    print("using CalculateDistanceToInterpEnvelop")
        # The ESMPY interpolator doesn't cover a significant part of the "na" domain's interior
    u0=u[:,0]; # grab first interpolation estimate to find valid footprint of interpolator
    dist2bnd=CalculateDistanceToInterpEnvelope(xi,yi,u0, 1.)
    with nc.Dataset(dist2bnd_file, 'w', format='NETCDF4') as ncout:
        ncout.createDimension('node' , nn)
        d_var=ncout.createVariable('dist2bnd', 'f4', ('node',))
        d_var.long_name     = 'distance to boundary'
        d_var.units         = 'km'
        d_var.standard_name = 'distance to boundary'
        d_var[:]=dist2bnd[:]
##################################################################################
# FINISHED: Compute distance to boundary for each node in mesh:
##################################################################################

##################################################################################
# START: Assign forecast error variance based on distance to boundary and inputs
##################################################################################
if LambdaBndyErrorVariance > 0. :
    um0=u[:,0]
    InteriorNodeList=np.where(um0**2>=0) #find points that are valid floats
    Variance0=VarianceLinearDistanceToBndy(InteriorNodeList, dist2bnd, ForecastErrorVariance, 10.*ForecastErrorVariance,LambdaBndyErrorVariance )
else :
    #Assign spatially constant error variance if LambdaBndyErrorVariance<=0
    Variance0=0.*xi + ForecastErrorVariance 
Variance=np.zeros((nt,nn))

# Need To add some depth dependence opposite to what is in for STOFS
# Estimates are bad in shallow water because of missing tides and other coastal procs
VarianceDeep=ForecastErrorVariance 
VarianceShallow=25.*ForecastErrorVariance
VarianceZ = VarianceShallow + (VarianceDeep-VarianceShallow)*(zi-50.)/(250.-50.)
js=np.where(VarianceZ>VarianceShallow)
VarianceZ[js]=VarianceShallow
jd=np.where(VarianceZ<VarianceDeep)
VarianceZ[jd]=VarianceDeep

#variance is maximum of distance to boundary and depth based term
VarianceZandD2B=np.maximum(Variance0,VarianceZ)

for k in range(nt):
    Variance[k,:]=VarianceZandD2B


##################################################################################
# FINISHED: Assign forecast error variance based on distance to boundary and inputs
##################################################################################

ne=ei.shape[0]

with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:

    ncout.createDimension('level' , 1)  
    ncout.createDimension('node' , nn)
    ncout.createDimension('element' , ne)
    ncout.createDimension('time', nt)
    ncout.createDimension('noel', 3)

    lon_var=ncout.createVariable('longitude', 'f8', ('node',))
    lon_var.units         = 'degree_east'
    lon_var.long_name     = 'longitude'
    lon_var.standard_name = 'longitude'
    lon_var.axis          = 'X'
    lon_var[:]=xi[:]

    lat_var=ncout.createVariable('latitude', 'f8', ('node',))
    lat_var.units         = 'degree_north'
    lat_var.long_name     = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.axis          = 'Y'
    lat_var[:]=yi[:]

    time_var=ncout.createVariable('time', 'f8', ('time',))
    time_var.units         = 'seconds since 1970-01-01 00:00:00.0 0:00'
    time_var.long_name     = 'verification time generated by wgrib2 function verftime()'
    time_var.standard_name = 'time'
    time_var.axis          = 'T'
    time_var.reference_time = 1777334400
    time_var.reference_date = '2026.04.28 00:00:00 UTC'
    time_var[:]=t1[:]

    tri_var=ncout.createVariable('tri', 'i4', ('noel','element'))
    tri_var.long_name     = 'element list'
    tri_var.standard_name = 'element list'
    tri_var[:]=np.transpose(ei)

    d_var=ncout.createVariable('dist2bnd', 'f4', ('node',))
    d_var.long_name     = 'distance to boundary'
    d_var.units         = 'km'
    d_var.standard_name = 'distance to boundary'
    d_var[:]=dist2bnd[:]

    u_var=ncout.createVariable('u-vel', 'f4', ('time','node'),fill_value    = fill_value0)
    u_var.standard_name     = 'eastward_sea_water_velocity'
    u_var.units         = 'm/s'
    u_var.long_name = ' u-veloc. [93.1H]'
    u_var.level = 'layer 1 from:'+flin
    u_var[:,:]=u[:,:].T

    v_var=ncout.createVariable('v-vel', 'f4', ('time','node'),fill_value    = fill_value0)
    v_var.standard_name     = 'northward_sea_water_velocity'
    v_var.units         = 'm/s'
    v_var.long_name = ' v-veloc. [93.1H]''northward_wind'
    v_var.level = 'layer 1'
    v_var[:,:]=v[:,:].T

    ErrorVariance_var=ncout.createVariable('ErrorVariance', 'f4', ('time','node'),fill_value    = fill_value0)
    ErrorVariance_var.long_name     = 'forecast error variance'
    ErrorVariance_var.units         = 'm m /s /s'
    ErrorVariance_var.standard_name = 'variance'
    ErrorVariance_var.level = 'layer 1 from'+flin
    ErrorVariance_var[:,:]=Variance[:,:]

    ncout.close
standard_name = 'northward_sea_water_velocity'
long_name     = ' v-veloc. [93.1H]'
units         = 'm/s'
coordinates   = 'Latitude Longitude'
