import numpy as np
import netCDF4 as nc
import sys
import GeoComputeMeshToMeshInterpWeights as mshint
import InterpNWPSUtility as  mshutil
import xarray as xr
import scipy.sparse as sp
from scipy.interpolate import NearestNDInterpolator
import datetime

# Engine for interpolating to WW3 unstructured mesh using precomputed interpolation weights from netcdf files with forecasts
#
# to call:
# python InterpolateSTOFS.py input_file meshpath outputfile variable1:variable2:variable3 ExtrapMethod
#
# example:
# python InterpolateSTOFS.py stofs.20260608.00/stofs.cwl.vel.nc meshes/RWPS.V0a.small.msh tesdtoZ.vel.nc u-vel:v-vel 2
# or:
# python InterpolateSTOFS.py stofs.20260608.00/stofs.cwl.nc meshes/RWPS.V0a.small.msh tesdtoZ.vel.nc zeta 1
#
# ExtrapMethod =-2 missing values in source field replaced with 0
# ExtrapMethod =-1 no extrapolation, NaN's potentially in output where source field is dry
# ExtrapMethod = 0 NaN values in interpolated field replaced with 0.0
# ExtrapMethod = 1 Nearest Neighbor extrapolation from valid source values
# ExtrapMethod = 2 Nearest Neighbor extrapolation from valid interpolated values
# ExtrapMethod = 3 Nearest Neighbor extrapolation from interpolated nodes which allways have valid values (faster)

UseUnixTime=True
nargin = len(sys.argv) - 1

flin=sys.argv[1]

IsVel=(flin[-6:-3]=="vel")

mshfl=sys.argv[2]
meshslash=mshfl.rfind('/')+1
weights_file="STOFS.wght."+mshfl[meshslash:len(mshfl)-4]+".nc"

flout=sys.argv[3]
varname0=sys.argv[4]
varname=varname0.split(":")
#
ExtrapMethod=0 # no extrapolation
if nargin>4:
    ExtrapMethod=int(sys.argv[5])
    
if ExtrapMethod==-2:
    print("missing values in source field replaced with 0 and remaining missing values in interpolated field (outside of source footprint) with value 0 ")
if ExtrapMethod==-1:
    print("no extrapolation, nan left in place in output")
if ExtrapMethod==0:
    print("Fill missing values in interpolated field with value 0")
if ExtrapMethod==1:
    print("extrapolation from nearest valid point in source- can be slow if source mesh is much larger than destination mesh")
if ExtrapMethod==2:
    print("extrapolation from nearest valid point in destination (interpolated field)")

with xr.open_dataset(weights_file) as ds_s:
   # Standard sparse storage uses 'row', 'col', and 'data' variables
   row = ds_s['row'].values
   col = ds_s['col'].values
   weights = ds_s['S'].values
   Nrows=ds_s.attrs.get('Nrows')
   Ncols=ds_s.attrs.get('Ncols')
matrix = sp.coo_matrix((weights, (row-1, col-1)), shape=(Nrows,Ncols)).tocsr()
print("sparse interpolation matrix")
print(matrix)

xi, yi, ei = mshutil.loadWW3Mesh(mshfl)
nni=len(xi)

data = nc.Dataset(flin,"r")
if UseUnixTime:
    time=mshutil.ConvertTimeToUnixTime(flin)
else:
    time=np.asarray(data["time"][:])

x=np.asarray(data["x"][:])
#shift to RWPS convention
j=np.where(x>90)
x[j]=x[j]-360.
y=np.asarray(data["y"][:])

n1=len(x)
nt=len(time)

#nt=4
#time=time[0:nt]

print(time)

nvar=len(varname)
vari=np.zeros((nvar,nt,nni))

if ( ExtrapMethod>0 or  ExtrapMethod==-2):
    IsExtrap=np.zeros((nvar,nt,nni),dtype=int)
    
if ExtrapMethod==3:
    AnyExtrap=np.zeros((nvar,nni),dtype=int)

nan=float("nan")
    
for jv in range(nvar):
    fill_value0=data[varname[jv]]._FillValue
    print("fill value="+str(fill_value0))
    for k in range(nt):
        print("interpolating for time step = "+str(k)+" of "+str(nt))
        var=np.asarray(data[varname[jv]][k,:])
        #replace fill with nan to avoid interpolating fill
        j=np.where(var==fill_value0)
        var[j]=nan
        if ExtrapMethod==-2:
            j=np.where(np.isnan(var))
            var[j]=0.
        vari[jv,k,:] = matrix @ var
        if ExtrapMethod>0 and ExtrapMethod<3:
            jd=np.where(np.isnan(vari[jv,k,:]))
            dstp=np.array((xi[jd],yi[jd]))
            if ExtrapMethod==1:
            #extrapolate using nearest neighbor of source with valid value
                js=np.where(~np.isnan(var))
                srcp=np.array((x[js],y[js]))
                srcv=var[js]
            if ExtrapMethod==2:
            #extrapolate using nearest neighbor of interpolated field with valid value
                js=np.where(~np.isnan(vari[jv,k,:]))
                srcp=np.array((xi[js],yi[js]))
                tmp=vari[jv,k,js]
                srcv=tmp.flatten()
            interp = NearestNDInterpolator(srcp.T,srcv)
            ExtrapVals = interp( dstp.T )
            vari[jv,k,jd]=ExtrapVals
            IsExtrap[jv,k,jd]=1
        if ExtrapMethod==3: # Fast posthoc extrapolator
            jd=np.where(np.isnan(vari[jv,k,:]))
            AnyExtrap[jv,jd]=1.

if ( ExtrapMethod==0 or ExtrapMethod==-2 ):
    jd=np.where(np.isnan(vari))
    vari[jd]==0.
    IsExtrap[jd]=1
    
if ExtrapMethod==3: #posthoc extrapolation from points which are valid at all times
    for jv in range(nvar):
        jd=np.where(AnyExtrap[jv,:]==1) # nodes that have some "nan" intrepolated values
        js=np.where(AnyExtrap[jv,:]==0) # nodes that have no "nan" intrepolated values
        srcp = np.array((xi[js],yi[js])).T
        srcv = vari[jv,0,js] #dummy input field
        dstp = np.array((xi[jd],yi[jd])).T
        srcv=srcv[0,:]
        interpolator = NearestNDInterpolator(srcp, srcv)
        distances, jsrc = interpolator.tree.query(dstp)
        jd=jd[0]
        for k in range(nt):
            jdk=np.where(np.isnan(vari[jv,k,jd]))
            jdk=jdk[0]
            vari[jv,k,jd[jdk]]=vari[jv,k,jsrc[jdk]]
            IsExtrap[jv,k,jd[jdk]]=1

            
            
print("nn(target mesh) = "+str(nni)+": Nrows = "+str(Nrows))
print("nn(source mesh) = "+str(n1)+": Ncols = "+str(Ncols))
if not ((nni==Nrows) and (n1==Ncols)):
    print("WARNING: Wrong matrix weights: number of rows from "+ mshfl +" = "+str(nni)+
    " but number of rows in "+ weights_file +" = "+str(Nrows)+ 
    ", number of spatial points in "+ flin +" = "+str(n1)+ 
    " but number of columns in "+ weights_file +" = "+str(Ncols)  )
    print("  You may need to regnerate file "+ weights_file +" with appropriate weights")


ne=ei.shape[0]
with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:

    ncout.createDimension('level' , 1)  
    ncout.createDimension('node' , nni)
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
    
    if UseUnixTime:
        units = 'seconds since 1970-01-01 00:00:00.0 0:00'
        long_name = 'verification time generated by wgrib2 function verftime()'
        standard_name = 'time'
    else:
        varin = data["time"]
        units = varin.units 
        standard_name = varin.standard_name
        long_name = varin.long_name
    time_var=ncout.createVariable('time', 'f8', ('time',))
    time_var.units         = units
    time_var.long_name     = long_name 
    time_var.standard_name = standard_name
    time_var[:]=time[:]

    tri_var=ncout.createVariable('tri', 'i4', ('noel','element'))
    tri_var.long_name     = 'element list'
    tri_var.standard_name = 'element list'
    tri_var[:]=np.transpose(ei)

    for jv in range(nvar):
        print("writing output for :"+varname[jv])
        varin = data[varname[jv]]
        units = varin.units 
        standard_name = varin.standard_name
        long_name = varin.long_name
        location = 'node' 

        F_var=ncout.createVariable(varname[jv], 'f8', ('time','node'),fill_value = fill_value0)
        F_var.long_name     = long_name
        F_var.units         = units
        F_var.standard_name = standard_name
        F_var.location      = location
        F_var[:,:]          = vari[jv,:,:]
        
        if ExtrapMethod >= 0 :
            xtrp_var=ncout.createVariable(varname[jv]+'IsExtrap', 'i1', ('time','node'))
            xtrp_var.long_name     = '==1 if the interpolated value extrapolated. 0 if interpolated'
            xtrp_var.standard_name = 'is extrapolated'
            xtrp_var.location      = 'node'
            if ExtrapMethod == 0:
                xtrp_var.method        = 'Interpolated nan values replaced with 0'
            if ExtrapMethod == 1:
                xtrp_var.method        = 'nearest valid neighbor in source field'
            if ExtrapMethod == 2:
                xtrp_var.method        = 'nearest valid neighbor in interpolated field'
            xtrp_var[:,:]          = IsExtrap[jv,:,:]

    
    ncout.close




"""


Source:
           /mnt/sda/keston/STOFSInterp/stofs.20260608.00/stofs.cwl.nc
Format:
           netcdf4_classic
Global Attributes:
           _FillValue        = -99999
           model             = 'ADCIRC'
           version           = 'noaa.stofs.2d.glo.v2.1.0r1.v55.12'
           git_hash          = '23947fbd9683d0ef48f12e6ce62d45d18bc27ff3'
           grid_type         = 'Triangular'
           description       = '2026060800 :-6 hr nowcast and +180 hr forecast ! 32 CHARACTER ALPHANUMERIC RUN D'
           agrid             = 'OceanMesh2D'
           rundes            = '2026060800 :-6 hr nowcast and +180 hr forecast ! 32 CHARACTER ALPHANUMERIC RUN D'
           runid             = 'STOFS 2D GLOBAL v5.6.5     ! 24 CHARACTER ALPHANUMERIC RUN IDENTIFICATION'
           title             = 'STOFS_2D_GLOBAL.V2.1.0     ! NCPROJ - PROJECT TITLE'
           institution       = 'NOS/OCS/CSDL/CMMB          ! NCINST - PROJECT INSTITUTION'
           source            = 'Dogwood/Cactus             ! NCSOUR - PROJECT SOURCE'
           history           = 'PRODUCTION                 ! NCHIST - PROJECT HISTORY'
           references        = 'http://www.adcirc.org      ! NCREF  - PROJECT REFERENCES'
           comments          = 'STOFS_2D_GLOBAL.V2.1.0     ! NCCOM  - PROJECT COMMENTS'
           host              = 'NOS/OCS/CSDL/CMMB          ! NCHOST - PROJECT HOST'
           convention        = 'CF-1.0                     ! NCCONV - CONVENTIONS'
           Conventions       = 'UGRID-0.9.0'
           contact           = 'Yuji.Funaoshi@noaa.gov     ! NCCONT - CONTACT INFORMATION'
           creation_date     = '2026-06-08  3:55:15  00:00'
           modification_date = '2026-06-08  3:55:15  00:00'
           fort.15           = '==== Input File Parameters (below) ===='
           dt                = 6
           ihot              = 568
           ics               = 22
           nolibf            = 1
           nolifa            = 2
           nolica            = 1
           nolicat           = 1
           nwp               = 7
           ncor              = 1
           ntip              = 2
           nws               = 10
           nramp             = 1
           tau0              = 0.053333
           statim            = 0
           reftim            = 0
           rnday             = 794.5
           dramp             = 6.75
           a00               = 0.8
           b00               = 0.2
           c00               = 0
           h0                = 0.1
           slam0             = 0
           sfea0             = 45
           cf                = 0.0005
           eslm              = -0.2
           cori              = 0
           ntif              = 8
           nbfr              = 0
Dimensions:
           time      = 186   (UNLIMITED)
           node      = 12785004
           nele      = 24875336
           nvertex   = 3
           nbou      = 262
           nvel      = 12421
           max_nvell = 1772
           mesh      = 1
Variables:
    time       
           Size:       186x1
           Dimensions: time
           Datatype:   double
           Attributes:
                       long_name     = 'model time'
                       standard_name = 'time'
                       units         = 'seconds since 2024-04-04 12:00:00        ! NCDASE - BASE_DAT'
                       base_date     = '2024-04-04 12:00:00        ! NCDASE - BASE_DATE'
    x          
           Size:       12785004x1
           Dimensions: node
           Datatype:   double
           Attributes:
                       long_name     = 'longitude'
                       standard_name = 'longitude'
                       units         = 'degrees_east'
                       positive      = 'east'
    y          
           Size:       12785004x1
           Dimensions: node
           Datatype:   double
           Attributes:
                       long_name     = 'latitude'
                       standard_name = 'latitude'
                       units         = 'degrees_north'
                       positive      = 'north'
    element    
           Size:       3x24875336
           Dimensions: nvertex,nele
           Datatype:   int32
           Attributes:
                       long_name   = 'element'
                       cf_role     = 'face_node_connectivity'
                       start_index = 1
                       units       = 'nondimensional'
    adcirc_mesh
           Size:       1x1
           Dimensions: mesh
           Datatype:   int32
           Attributes:
                       long_name              = 'mesh_topology'
                       cf_role                = 'mesh_topology'
                       topology_dimension     = 2
                       node_coordinates       = 'x y'
                       face_node_connectivity = 'element'
    nvel       
           Size:       1x1
           Dimensions: 
           Datatype:   int32
           Attributes:
                       long_name = 'total number of normal flow specified boundary nodes including both the front and back nodes on internal barrier boundaries'
                       units     = 'nondimensional'
    nvell      
           Size:       262x1
           Dimensions: nbou
           Datatype:   int32
           Attributes:
                       long_name = 'number of nodes in each normal flow specified boundary segment'
                       units     = 'nondimensional'
    max_nvell  
           Size:       1x1
           Dimensions: 
           Datatype:   int32
    ibtype     
           Size:       262x1
           Dimensions: nbou
           Datatype:   int32
           Attributes:
                       long_name = 'type of normal flow (discharge) boundary'
                       units     = 'nondimensional'
    nbvv       
           Size:       12421x1
           Dimensions: nvel
           Datatype:   int32
           Attributes:
                       long_name = 'node numbers on normal flow boundary segment'
                       units     = 'nondimensional'
    depth      
           Size:       12785004x1
           Dimensions: node
           Datatype:   double
           Attributes:
                       long_name     = 'distance  below geoid'
                       standard_name = 'depth below geoid'
                       coordinates   = 'time y x'
                       location      = 'node'
                       mesh          = 'adcirc_mesh'
                       units         = 'm'
    zeta       
           Size:       12785004x186
           Dimensions: node,time
           Datatype:   double
           Attributes:
                       long_name     = 'water surface elevation above geoid'
                       standard_name = 'sea_surface_height_above_geoid'
                       coordinates   = 'time y x'
                       location      = 'node'
                       mesh          = 'adcirc_mesh'
                       units         = 'm'
                       _FillValue    = -99999
                       
                       
"""                       
