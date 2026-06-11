import numpy as np
import netCDF4 as nc
import sys
import GeoComputeMeshToMeshInterpWeights as mshint
import InterpNWPSUtility as  mshutil
import xarray as xr
import scipy.sparse as sp
from scipy.interpolate import NearestNDInterpolator

def copy_netcdf(source_path, target_path):
    # Open the source file in read mode and create the target file in write mode
    with nc.Dataset(source_path, 'r') as src, nc.Dataset(target_path, 'w') as trg:
        
        # 1. Copy structural dimensions
        for name, dim in src.dimensions.items():
            # Check if the dimension is unlimited (like time) to maintain its behavior
            dim_len = len(dim) if not dim.isunlimited() else None
            trg.createDimension(name, dim_len)

        # 2. Copy global file attributes
        trg.setncatts({attr: src.getncattr(attr) for attr in src.ncattrs()})

        # 3. Copy variables with their metadata and data arrays
        for name, src_var in src.variables.items():
            # Create the variable structure in the new file
            trg_var = trg.createVariable(
                name, 
                src_var.datatype, 
                src_var.dimensions,
                fill_value=getattr(src_var, '_FillValue', None) # Safely handles fill values
            )
            
            # Copy specific variable attributes (e.g., units, long_name)
            var_attrs = {attr: src_var.getncattr(attr) for attr in src_var.ncattrs() if attr != '_FillValue'}
            trg_var.setncatts(var_attrs)
            
            # Copy the actual numerical/data arrays
            trg_var[:] = src_var[:]

def NearestPoint(x,y,xi,yi):
#return index of nearest (x,y) for each (xi.yi)
    deg2kmY=111.
    n=len(xi)
    indx=np.zeros(n,dtype=int)
    for k in range(n):
        deg2kmX=np.cos( np.pi * yi[k] / 180.)*deg2kmY
        distances=np.abs( deg2kmX*(xi[k] - x) + 1j*deg2kmY*(yi[k] - y))
        indx[k]=np.argmin(distances)
    return indx
    

flin=sys.argv[1]

mshfl=sys.argv[2]
meshslash=mshfl.rfind('/')+1
weights_file="STOFS.wght."+mshfl[meshslash:len(mshfl)-4]+".nc"

flout=sys.argv[3]
varname=sys.argv[4]


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

#var=np.asarray(data[varname][:,:])
time=np.asarray(data["time"][:])
x=np.asarray(data["x"][:])
j=np.where(x>90)
x[j]=x[j]-360.
y=np.asarray(data["y"][:])
n1=len(x)
nt=len(time)

vari=np.zeros((nt,nni))
IsExtrap=np.zeros((nt,nni),dtype=int)
#for i in range(0, 50, 7):
fill_value0=data[varname]._FillValue
print("fill value="+str(fill_value0))
nan=float("nan")
extrapolate=False
extrapolate=True
for k in range(nt):
    print("interpolating for time step = "+str(k)+" of "+str(nt))
    var=np.asarray(data[varname][k,:])
    j=np.where(var==fill_value0)
    var[j]=nan
    #_FillValue    = -99999
    print(var.shape)
    print(vari.shape)
    print(matrix.shape)
    vari[k,:] = matrix @ var
    if extrapolate:
        #extrapolate using nearest neighbor
        jd=np.where(np.isnan(vari[k,:]))
        js=np.where(~np.isnan(var))
        srcp=np.array((x[js],y[js]))
        srcv=var[js]
        dstp=np.array((xi[jd],yi[jd]))
        interp = NearestNDInterpolator(srcp.T,srcv)
        ExtrapVals = interp( dstp.T )
        vari[k,jd]=ExtrapVals
        IsExtrap[k,jd]=1

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

    time_var=ncout.createVariable('time', 'f8', ('time',))
    time_var.units         = 'seconds since 2024-04-04 12:00:00        ! NCDASE - BASE_DAT'
    time_var.long_name     = 'model time'
    time_var.standard_name = 'time'
    time_var.axis          = 'T'
    time_var[:]=time[:]

    tri_var=ncout.createVariable('tri', 'i4', ('noel','element'))
    tri_var.long_name     = 'element list'
    tri_var.standard_name = 'element list'
    tri_var[:]=np.transpose(ei)

    zeta_var=ncout.createVariable(varname, 'f8', ('time','node'),fill_value    = fill_value0)
    zeta_var.long_name     = 'water surface elevation above geoid'
    zeta_var.units         = 'm'
    zeta_var.standard_name = 'sea_surface_height_above_geoid'
    zeta_var.location = 'node'
    zeta_var[:,:]=vari[:,:]
    
    xtrp_var=ncout.createVariable('IsExtrap', 'i4', ('time','node'),fill_value    = fill_value0)
    xtrp_var.long_name     = '==1 if the interpolated value extrapolated'
    xtrp_var.standard_name = 'is extrapolated'
    xtrp_var.location = 'node'
    xtrp_var[:,:]=IsExtrap[:,:]
    
    
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
