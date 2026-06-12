import numpy as np
import time
from datetime import datetime


import numpy as np
#import matplotlib.tri as mtri
"""
def compute_and_save_weights(x, y, triangles, x_query, y_query):
   # Computes and extracts interpolation weights and triangle indices.
    
   # :param x, y: 1D arrays of mesh vertex coordinates.
   # :param triangles: (M, 3) array of triangle connectivity.
   # :param x_query, y_query: 1D arrays of query points.
   # :return: triangle_indices, weights, valid_mask

    # 1. Initialize the triangulation and interpolator
    triang = mtri.Triangulation(x, y, triangles)
    interp = mtri.LinearTriInterpolator(triang, np.zeros_like(x))
    
    # 2. Extract the default TriFinder from the interpolator
    trifinder = interp.get_trifinder()
    
    # 3. Find which triangle contains each query point (-1 means outside mesh)
    triangle_indices = trifinder(x_query, y_query)
    valid_mask = triangle_indices != -1
    
    # Filter points that actually fall inside the mesh
    tri_ids = triangle_indices[valid_mask]
    xq = x_query[valid_mask]
    yq = y_query[valid_mask]
    
    # 4. Extract mesh vertices for the containing triangles
    # node_indices shape: (num_valid_points, 3)
    node_indices = triang.triangles[tri_ids]
    
    x0, y0 = x[node_indices[:, 0]], y[node_indices[:, 0]]
    x1, y1 = x[node_indices[:, 1]], y[node_indices[:, 1]]
    x2, y2 = x[node_indices[:, 2]], y[node_indices[:, 2]]
    
    # 5. Calculate Barycentric Coordinates (Weights)
    # Exploit the transformation matrix determinant
    denom = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
    
    # Handle collinear/degenerate triangles safely
    denom = np.where(denom == 0, 1e-15, denom)
    
    w0 = ((y1 - y2) * (xq - x2) + (x2 - x1) * (yq - y2)) / denom
    w1 = ((y2 - y0) * (xq - x2) + (x0 - x2) * (yq - y2)) / denom
    w2 = 1.0 - w0 - w1
    
    # Pack weights into a (num_valid_points, 3) array
    weights = np.column_stack((w0, w1, w2))
    
    return triangle_indices, weights, valid_mask


def apply_saved_weights(triangles, field_values, triangle_indices, weights, valid_mask):
   # Instantly applies saved weights to interpolate a new field.
    result = np.full(len(triangle_indices), np.nan)
    
    if not np.any(valid_mask):
        return result
        
    tri_ids = triangle_indices[valid_mask]
    node_indices = triangles[tri_ids]
    
    # Gather field values at the 3 vertices of each containing triangle
    v_values = field_values[node_indices]
    
    # Multiply elements and sum across the 3 nodes
    result[valid_mask] = np.sum(weights * v_values, axis=1)
    
    return result

"""

def IsInElement(x,y,xp,yp):
    IsIn=False
    c1 = (x[1] - x[0]) * (yp - y[0]) - (y[1] - y[0]) * (xp - x[0])
    c2 = (x[2] - x[1]) * (yp - y[1]) - (y[2] - y[1]) * (xp - x[1])
    c3 = (x[0] - x[2]) * (yp - y[2]) - (y[0] - y[2]) * (xp - x[2])
    if ( ( c1 > 0 and c2 > 0 and c3 > 0) or ( c1 < 0 and c2 < 0 and c3 < 0) ):
        IsIn=True
    if (c1*c2*c3 == 0): # include points on triangle
        IsIn=True
    return IsIn

def FindElement(x,y,e,xi,yi):
    nd=len(e.shape)
    e=np.squeeze(e)
    j=-9999
    ne=e.shape[0]
    #print(nd)
    #print(ne)
    #print(e)
    if nd > 1: 
        xc = np.squeeze(np.mean(x[e], axis=1))
        yc = np.squeeze(np.mean(y[e], axis=1))
        DistanceToElements = np.abs((xi + 1j*yi) - (xc + 1j*yc))
        n=0
        while (j < 0 and n < ne)  :
            n=n+1
            j = np.argmin( DistanceToElements )
            xl = np.squeeze(x[e[j,:]])
            yl = np.squeeze(y[e[j,:]])
            IsIn=IsInElement(xl,yl,xi,yi)
            if IsIn :
                return j
            else:
                DistanceToElements[j]=float('inf')
                j=-9999
    else:
        xl = np.squeeze(x[e])
        yl = np.squeeze(y[e])
        IsIn=IsInElement(xl,yl,xi,yi)
        if IsIn:
            j=0
            return j
        else:
            j=-9999
            return j
        
    return j
"""
import matplotlib.tri as mtri
def compute_mesh_to_mesh_interp_weights(x, y, e, xi, yi, flout):

# Assume x, y are vertex coords (N,), triangles is (M, 3) connectivity array
    e0=e-1
    triang = mtri.Triangulation(x, y, e0)

    # z are the known values at each vertex (N,)
    linear_interp = mtri.LinearTriInterpolator(triang, z)

    # Evaluate at new points (callable directly)
    z_new = linear_interp(x_new, y_new)
    
"""    

def universal_len(obj):
    try:
        return len(obj)
    except TypeError:
        return 1


def compute_mesh_to_mesh_interp_weights(x, y, e, xi, yi):
    """
    Compute mesh-to-mesh interpolation weights using barycentric coordinates.
    
    Parameters:
    -----------
    x : array-like
        X coordinates of source mesh nodes
    y : array-like
        Y coordinates of source mesh nodes
    e : array-like
        Element connectivity matrix (n_elements x 3 for triangular elements)
        (nodes indexed from 1 to nn)
    xi : array-like
        X coordinates of target points
    yi : array-like
        Y coordinates of target points
        
    Returns:
    --------
    weights : ndarray
        Interpolation weights (n_points x 3)
    nodes : ndarray
        Node indices for each point (n_points x 3)
    elenum : ndarray
        Element number for each point (n_points,)
    Dist2EleCenter : ndarray
        Distance to element center for each element
    """
    
    # Convert to numpy arrays
    
    deg2kmY=111.
    UseNearestEle = False #if True just use closest element center
    N = 12 # search N nearest elements (nearest by element center to target node distance)
            
    print(xi)
    x = np.asarray(x)
    y = np.asarray(y)
    e = np.asarray(e)
    xi = np.asarray(xi)
    yi = np.asarray(yi)
    print("e.shape")
    print(e.shape)
    # convert node indexs to 0 .. nn-1
    e0=e-1
    

#    triangle_indices, weights, valid_mask=compute_and_save_weights(x, y, e0, xi, yi)
#    np.savetxt("triangle_indices.txt",triangle_indices)
#    np.savetxt("weights.txt",weights)
#    np.savetxt("mask.txt",valid_mask)

    xc = np.mean(x[e0], axis=1)
    yc = np.mean(y[e0], axis=1)
    
    # Initialize output arrays
    n_points = len(xi)
    n_elements = len(e)
    weights = np.zeros((n_points, 3))
    nodes = np.zeros((n_points, 3), dtype=int)
    elenum = np.zeros(n_points, dtype=int)
    Dist2EleCenter = np.zeros(n_points)
    
    t0 = time.time()

    for k in range(n_points):
        x0=xi[k]
        y0=yi[k]
        deg2kmX=np.cos( np.pi * y0 / 180.)*deg2kmY
        #distances = np.abs((x0 + 1j*y0) - (xc + 1j*yc))
        distances = np.abs( deg2kmX*(x0 - xc) + 1j*deg2kmY*(y0 - yc))
        if UseNearestEle:
            jg=np.argmin(distances)
        else:
            raw_indices = np.argpartition(distances, N)[:N]
            sorted_sub_indices = np.argsort(distances[raw_indices])
            jg = raw_indices[sorted_sub_indices]
        if UseNearestEle :
            j=jg
        else:
            j=FindElement(x,y,e0[jg,:],x0,y0)
            if j < 0:
                print("couldn't find element for point : "+ str(k))
                print(str(x0)+" : "+str(y0))
                j = jg[0]
                print("using closest element at distance : "+ str(distances[j]))
            else:
                j=jg[j]

        elenum[k] = j
        Dist2EleCenter[k] = distances[j]
        xl = x[e0[j, :]]
        yl = y[e0[j, :]]
            
        # Compute barycentric coordinates using shoelace formula for areas
        # Area of triangle formed by vertices 1, 2, and point (a3)
        xt = np.array([xl[0], xl[1], x0, xl[0]])
        yt = np.array([yl[0], yl[1], y0, yl[0]])
        a3 = -np.dot(xt[1:4] - xt[0:3], yt[0:3] + yt[1:4]) / 2
        # Area of triangle formed by vertices 3, 1, and point (a2)
        xt = np.array([xl[2], xl[0], x0, xl[2]])
        yt = np.array([yl[2], yl[0], y0, yl[2]])
        a2 = -np.dot(xt[1:4] - xt[0:3], yt[0:3] + yt[1:4]) / 2
        # Area of triangle formed by point, vertices 2, 3 (a1)
        xt = np.array([x0, xl[1], xl[2], x0])
        yt = np.array([y0, yl[1], yl[2], y0])
        a1 = -np.dot(xt[1:4] - xt[0:3], yt[0:3] + yt[1:4]) / 2

        # Normalize to get barycentric weights
        total_area = a1 + a2 + a3
        weights[k, :] = [a1, a2, a3] / total_area

        nodes[k, :] = e[j, :] #<- indexed 1 .. nn
        # Progress reporting every 100 iterations
        if (k + 1) % 10 == 0:
            t1 = time.time()
            time_per_iter = (t1 - t0) / (k + 1)
            time_remaining = (n_points - k - 1) * time_per_iter / 60
            print(f"Progress: {k+1}/{n_points}")
            print(f"Time remaining: {time_remaining:.2f} minutes")
    return weights, nodes, elenum, Dist2EleCenter


def InterpolateField2Nodes(nodes,weights, f):
    fi=np.zeros(weights.shape[0])
    for k in range(weights.shape[0]):
        fl=f[nodes[k,:]]
        wl=weights[k,:]
        fi[k]=np.dot(wl,fl)
        if (not np.abs(fi[k]) > 0.):
           fi[k]=0.
        if ( fi[k]>np.max(fl)  ):
           fi[k]=np.max(fl)
        if ( fi[k]<np.min(fl)  ):
           fi[k]=np.min(fl)
    return fi



def WriteInterpJobscript(fl,flin,mshfl,Njobs, ComputeNodes):
    
    meshslash=mshfl.rfind('/')+1
    TmpOutDir="STOFSInterpWeights."+mshfl[meshslash:len(mshfl)-4]
    WghtFl="STOFS.wght."+mshfl[meshslash:len(mshfl)-4]+".txt"
    WghtFlNetCDF="STOFS.wght."+mshfl[meshslash:len(mshfl)-4]+".nc"
    with open(fl, 'w') as f:
        #yi[k]=float(values[4])
        f.write("#!/bin/bash \n")
        f.write("#SBATCH --job-name=STOFS_interp_masterscript \n")
#        f.write("#SBATCH --ntasks="+str(N)+" \n")
        f.write("#SBATCH --ntasks=1 \n") # ntasks per interpolation
        f.write("#SBATCH --time=08:00:00 \n") 
        f.write("#SBATCH --output=mpi_test_%j.log \n")
        f.write("#SBATCH --error=%j.err \n")
        f.write("#SBATCH --account=marine-cpu \n")
        f.write("#SBATCH --nodes="+str(ComputeNodes)+" \n")
        f.write("#SBATCH --ntasks-per-core=1"+" \n")
        f.write("#SBATCH --array=0-"+str(Njobs-1)+" \n")

        f.write(" \n")

        f.write("module purge \n")
        f.write("module use /scratch4/NCEPDEV/marine/Ali.Salimi/Hera_Data/HR4-OPT/FromJessica/Keston/ICunstructuredRuns15km-implicit-450s/global-workflow/sorc/ufs_model.fd/modulefiles \n")
        f.write("module load ufs_ursa.intel \n")
        f.write("module load py-scipy/1.14.1 \n")
        f.write("module load py-netcdf4/1.7.1.post2 \n")
        f.write("pip list \n")
        f.write("# calculate interpolation weights in parallel geographically \n")
        f.write("srun python GeoSubsetInterpolateSTOFS.py "+flin+" "+mshfl+" $SLURM_ARRAY_TASK_ID " + str(Njobs)+" > InterpJob.$SLURM_ARRAY_TASK_ID.out \n")
        f.write("wait\n")
        f.write("# concatonate different parts of the mesh to common text file \n")
        f.write("cat "+TmpOutDir+"/Part.IntrpWghts.*.txt > "+WghtFl)
        f.write("# convert output weights to netcdf file \n")
        f.write("python ConvertWeights2NetCDF.py "+flin+" "+mshfl+" \n")

            
    flintrp="STOFS.to."+mshfl[meshslash:len(mshfl)-4]+".sh"
    flout=flin[0:-3]+".rwps.nc"
    flinuv=flin[0:-3]+".vel.nc"
    floutu=flinuv[0:-3]+".rwps.u.nc"
    floutv=flinuv[0:-3]+".rwps.v.nc"
    with open(flintrp, 'w') as f:
        #yi[k]=float(values[4])
        f.write("#!/bin/bash \n")
        f.write("#SBATCH --job-name=STOFS_interp_masterscript \n")
#        f.write("#SBATCH --ntasks="+str(N)+" \n")
        f.write(" \n")
        f.write("module purge \n")
        f.write("module use /scratch4/NCEPDEV/marine/Ali.Salimi/Hera_Data/HR4-OPT/FromJessica/Keston/ICunstructuredRuns15km-implicit-450s/global-workflow/sorc/ufs_model.fd/modulefiles \n")
        f.write("module load ufs_ursa.intel \n")
        f.write("module load py-scipy/1.14.1 \n")
        f.write("module load py-netcdf4/1.7.1.post2 \n")
        f.write("pip list \n")
        f.write("# calculate interpolation weights in parallel geographically \n")
        f.write("python InterpolateSTOFS.py "+flin+" "+mshfl+" "+flout+" zeta 2\n")
        f.write("python InterpolateSTOFS.py "+flinuv+" "+mshfl+" "+floutu+" u-vel:v-vel 2\n")

#-rw-r----- 1 keston keston 12G Jun  9 14:45 stofs.20260608.00/stofs.cwl.nc
#-rw-r----- 1 keston keston 25G Jun  9 14:58 stofs.20260608.00/stofs.cwl.vel.nc
