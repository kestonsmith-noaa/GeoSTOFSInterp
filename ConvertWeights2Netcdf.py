import numpy as np
import netCDF4 as nc
import sys
import GeoComputeMeshToMeshInterpWeights as mshint
import InterpNWPSUtility as  mshutil

nargin = len(sys.argv) - 1

flin=sys.argv[1]
mshfl=sys.argv[2]
meshslash=mshfl.rfind('/')+1

AddExtrapolationSupport=False
if nargin==3:
    if int(sys.argv[3])>0:
        AddExtrapolationSupport=True

TextWeightFl="STOFS.wght."+mshfl[meshslash:len(mshfl)-4]+".txt"
flout="STOFS.wght."+mshfl[meshslash:len(mshfl)-4]+".nc"

xi, yi, ei = mshutil.loadWW3Mesh(mshfl)
nn_dst=len(xi)
data = nc.Dataset(flin,"r")
#read spaital dimensions and determine if input mesh is curvilinear or regular
x=np.asarray(data["x"][:])
nn_src=len(x)

#Read in weights from cat of text output(rows are not in order)
F=np.loadtxt(TextWeightFl)

node_dst=F[:,0] # row number (destination node)
ele_num_src=F[:,1] # Source element number (not used)
node_src=F[:,2:5] # Source node
DistToCenterKM=F[:,5] # distance of destination node from source element center (not used) 
weights=F[:,6:9] # interpolation weights from source nodes to destination node

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
    
    s_var=ncout.createVariable('S', 'f8', ('n_s',))
    s_var.long_name     = 'matrix value'
    s_var[:]=val[:]

#Consider adding (x,y) destination and (x,y) for source. This is usefull for extrapolation
    if AddExtrapolationSupport:
        x=np.asarray(data["x"][:])
        y=np.asarray(data["y"][:])
        ncout.createDimension('nn_src' , nn_src)
        ncout.createDimension('nn_dst' , nn_dst)
        
        xsrc_var=ncout.createVariable('x_src', 'f8', ('nn_src',))
        xsrc_var.long_name     = 'interpolation source node longitude'
        xsrc_var[:]=x[:]
        
        ysrc_var=ncout.createVariable('y_src', 'f8', ('nn_src',))
        ysrc_var.long_name     = 'interpolation source node latitude'
        ysrc_var[:]=y[:]
        
        xdst_var=ncout.createVariable('x_dst', 'f8', ('nn_dst',))
        xdst_var.long_name     = 'interpolation destination node longitude'
        xdst_var[:]=xi[:]
    
        ydst_var=ncout.createVariable('y_dst', 'f8', ('nn_dst',))
        ydst_var.long_name     = 'interpolation destination node latitude'
        ydst_var[:]=yi[:]
   


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
