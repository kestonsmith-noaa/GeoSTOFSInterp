import time
import numpy as np
import netCDF4 as nc
import sys
from scipy.interpolate import RegularGridInterpolator
import GeoComputeMeshToMeshInterpWeights as mshint
import InterpNWPSUtility as  mshutil

import os
nargin = len(sys.argv) - 1

flin=sys.argv[1]
mshfl=sys.argv[2]

meshslash=mshfl.rfind('/')+1
TmpOutDir="STOFSInterpWeights."+mshfl[meshslash:len(mshfl)-4]

#if nargin==3 then write the parallel jobcard that carries out the interpolation weight calculation
if nargin ==3 :
    Njobs=int(sys.argv[3])
    os.makedirs(TmpOutDir, exist_ok=True)
    try:
        os.remove(TmpOutDir+"/*.txt")
    except:
        print("directory "+TmpOutDir+" is alread empty")
    mshint.WriteInterpJobscript("jobcardGeoInterpSTOFS",flin,mshfl,Njobs, 1)
    print("Made parallel jobcard to create interpolation weights with "+str(Njobs)+" processes. Next step:")
    print("sbatch jobcardGeoInterpSTOFS")
    sys.exit()

#if nargin==4 then do the part of the mesh for jobID of Njobs
if nargin == 4:
    jobID=int(sys.argv[3])
    Njobs=int(sys.argv[4])

#if nargin==6 then do the interpolation weight calculation for a specified window
else:
    lonW=int(sys.argv[3])
    lonE=int(sys.argv[4])
    latS=int(sys.argv[5])
    latN=int(sys.argv[6])

#latS=-40
#latN=81

meshslash=mshfl.rfind('/')+1
if nargin <  5:
    weights_file = TmpOutDir+"/Part.IntrpWghts."+str(jobID)+".txt"
else:
    weights_file = TmpOutDir+"/Part.IntrpWghts.W"+str(lonW)+".E"+str(lonE)+".S"+str(latS)+".N"+str(latN)+"."+mshfl[meshslash:len(mshfl)-3]+"txt"


print("saving output to:"+weights_file)

xi, yi, ei = mshutil.loadWW3Mesh(mshfl)
nni=len(xi)

# divide domain from east west based on job id=0 ... Njobs
# and make North-South window contain the full target domain
# NOTE: for efficiency the sections of the domain should have approximately
# the same number of destination nodes and/or the same number of source elements
# will implement as function shortly
if nargin <  5:
    Wdom=np.floor(np.min(xi))-1
    Edom=np.ceil( np.max(xi))+1
    lonW=Wdom+(Edom-Wdom)*jobID/(Njobs-1)
    lonE=Wdom+(Edom-Wdom)*(jobID+1)/(Njobs-1)

    latS=np.floor(np.min(yi))-1
    latN=np.ceil(np.max(yi))+1
    

data = nc.Dataset(flin,"r")
#read spaital dimensions and determine if input mesh is curvilinear or regular
x=np.asarray(data["x"][:])
y=np.asarray(data["y"][:])
e=np.asarray(data["element"][:,:])

#transform to match RWPS coordinates [-230W to +10E]
xp=x
j=np.where(xp>90.)
xp[j]=xp[j]-360

e0=e-1
xc = np.mean(x[e0], axis=1)
yc = np.mean(y[e0], axis=1)

nn=len(x)

SearchWidth=.5 # search rectangle extended beyond window to this number

print(e.shape)
ne=e.shape[0]
print(ne)

xwin=[lonW,lonE]
ywin=[latS,latN]
print("calculating interpolation weights for [W"+str(lonW)+": E "+str(lonE)+": S "+str(latS)+": N "+str(latN)+"]")
#Find target mesh nodes in window
jxUi=np.where( xi <= np.max(xwin))[0].tolist()
jxDi=np.where( xi >= np.min(xwin))[0].tolist()
jyUi=np.where( yi <= np.max(ywin))[0].tolist()
jyDi=np.where( yi >= np.min(ywin))[0].tolist()
ji=list( set(jxUi) & set(jxDi)  & set(jyUi)   & set(jyDi)  )

#Find source elements in and near window
SearchWidth=.5
jxU=np.where( xc < np.max(xwin)+SearchWidth )[0].tolist()
jxD=np.where( xc > np.min(xwin)-SearchWidth )[0].tolist()
jyU=np.where( yc < np.max(ywin)+SearchWidth )[0].tolist()
jyD=np.where( yc > np.min(ywin)-SearchWidth )[0].tolist()
je=list( set(jxU) & set(jxD)  & set(jyU)   & set(jyD)  )

print("number of target nodes for this region = "+str(len(ji)))
print("number of source elements in this region = "+str(len(je)))

weights, nodes, elenum, Dist2EleCenter=mshint.compute_mesh_to_mesh_interp_weights(x, y, e[je,:], xi[ji], yi[ji])

#Move back to global index's
for k in range(len(elenum)):
    elenum[k]=je[elenum[k]]
    ji[k]=ji[k]

#return to convention numbering from 1 ...
elenum = list(np.array(elenum) + 1)
ji = list(np.array(ji) + 1)
#nodes = np.array(nodes) + 1 # as writen using the correct 1 .. nn node numbering 

# File format 
# ["taget node number" "source element number" "source n1" "source n2 "source n3" "distance from target node to element center" "weight1" "weight2" "weight3"]

Fout=np.vstack((np.array(ji),elenum,nodes[:,0],nodes[:,1],nodes[:,2],Dist2EleCenter,weights[:,0],weights[:,1],weights[:,2]))
np.savetxt(weights_file,Fout.T, fmt='%d %d %d %d %d %.6f %.6f %.6f %.6f')


