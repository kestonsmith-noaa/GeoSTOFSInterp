import numpy as np
import netCDF4 as nc
import sys
import GeoComputeMeshToMeshInterpWeights as mshint
import InterpNWPSUtility as  mshutil

flin=sys.argv[1]
mshfl=sys.argv[2]
#flout=sys.argv[3]
#flout=TextWeightFl[1:len(mshfl)-4]+".nc"
#mshfl="meshes/RWPS.V0a.small.msh"
#flin="stofs.20260608.00/stofs.cwl.nc"
meshslash=mshfl.rfind('/')+1

TextWeightFl="STOFS.wght."+mshfl[meshslash:len(mshfl)-4]+".txt"
flout="STOFS.wght."+mshfl[meshslash:len(mshfl)-4]+".nc"

xi, yi, ei = mshutil.loadWW3Mesh(mshfl)
nn_dst=len(xi)
data = nc.Dataset(flin,"r")
#read spaital dimensions and determine if input mesh is curvilinear or regular
x=np.asarray(data["x"][:])
nn_src=len(x)

F=np.loadtxt(TextWeightFl)

node_dst=F[:,0]
ele_num_src=F[:,1]
node_src=F[:,2:5]
DistToCenterKM=F[:,5]
weights=F[:,6:9]

node_dst=node_dst.astype(int)
node_src=node_src.astype(int)
ele_num_src=ele_num_src.astype(int)

nnz=F.shape[0]
n_s=3*nnz
row=np.zeros(n_s,dtype=int)
col=np.zeros(n_s,dtype=int)
val=np.zeros(n_s)

n=0
for k in range(nnz):
    for j in range(3):
        row[n]=node_dst[k]
        col[n]=node_src[k,j]
        val[n]=weights[k,j]
        n=n+1
print("n="+str(n))
print("n_s="+str(n_s))
n_s=n

with nc.Dataset(flout, 'w', format='NETCDF4') as ncout:
    ncout.createDimension('n_s' , n_s)
    ncout.setncattr("Nrows", nn_dst)
    ncout.setncattr("Ncols", nn_src)
    
    r_var=ncout.createVariable('row', 'i4', ('n_s',))
    r_var.long_name     = 'row index'
    r_var[:]=row[:]
    
    c_var=ncout.createVariable('col', 'i4', ('n_s',))
    c_var.long_name     = 'column index'
    c_var[:]=col[:]
    
    s_var=ncout.createVariable('S', 'f4', ('n_s',))
    s_var.long_name     = 'matrix value'
    s_var[:]=val[:]


"""
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
    print("  You probably need to remove file "+ weights_file +" and rerun to generate appropriate weights")

matrix = sp.coo_matrix((weights, (row-1, col-1)), shape=(nn,n1)).tocsr()
print("sparse interpolation matrix")
print(matrix)
"""
