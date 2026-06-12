import os
import argparse
import time
import numpy as np
import netCDF4 as nc
import sys
from scipy.interpolate import RegularGridInterpolator
import re
import ComputeMeshToMeshInterpWeights as mshint
import csv


def loadWW3MeshCoords(fl):
    f=open(fl, 'r')
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() # number of nodes
    nn = re.findall(r'\d+', header)
    nn=int(nn[0])
    #print(nn)
    xi=np.zeros(nn)
    yi=np.zeros(nn)
    k=0
    for i in range(nn):
        A = f.readline()
        #print(A)
        values = A.split(" ")
        #print(values)
        xi[k]=values[2]
        yi[k]=values[4]
        k=k+1
    return xi, yi


def WriteNodeByTimeMat(fl,F):
    nn=F.shape[0]
    nt=F.shape[1]
    f=open(fl,'w')
    for n in range(nn):
        rowstr=" "
        for k in range(nt):
            rowstr =rowstr + " " + "{:.4f}".format(F[n,k])
        f.write(rowstr  + "\n")
    f.close


def ReadWeights(fl):
    nodes=[]
    weights=[]
    with open(fl, 'r', newline='') as csvfile:
        csv_reader = csv.reader(csvfile)
        for row in csv_reader:
            vals=np.array(row)
            nloc=[ int(vals[0])   ,int(vals[1])   ,int(vals[2]) ]
            wloc=[ float(vals[3]) ,float(vals[4]) ,float(vals[5]) ]
            nodes.append(nloc)
            weights.append(wloc)
    nodes=np.array(nodes)
    weights=np.array(weights)
    return nodes, weights



N=int(sys.argv[1])
#File locations and parameter specifications
mesh="global_1deg_unstr"
OutDir=mesh+".files/"
flmsh="meshes/"+mesh+".msh"
flinz="20241201/stofs_2d_glo.t00z.fields.htp.nc"
floutw=OutDir+"InterpWeights.global_1deg_unstr."+str(N)
floutz=OutDir+"zeta.global_1deg_unstr."+str(N)
ComputeWeights = False
#Only search elements within MaxDist x MaxDist box of interpolation point
#This is only used to speed up finding element the point is within
MaxDist=.5 #(degrees of lat and lon)

xi, yi =loadWW3MeshCoords(flmsh)
xi=xi % 360

print("processing partition: ",str(N))

fln=OutDir+'NodeList.'+str(N)+'.txt'
f=open(fln,"r")
header = f.readline()
header = f.readline() # number of nodes
nn = re.findall(r'\d+', header)
nn=int(nn[0])
#print(nn)
xil=np.zeros(nn)
yil=np.zeros(nn)
for k in range(nn):
    j = int(f.readline())
    xil[k]=xi[j]
    yil[k]=yi[j]
    print(str(k)+" " +str(j)+" "+str(xil[k]))

data = nc.Dataset(flinz,"r")
x=np.asarray(data["x"][:])
y=np.asarray(data["y"][:])

x=x % 360
xil=xil % 360

e=np.asarray(data["element"][:])
#print(e)

if ComputeWeights:
    print("computing interpolation weights to be stored in: " +floutw+".csv")
    mshint.compute_mesh_to_mesh_interp_weights(x, y, e, xil, yil, floutw, MaxDist)
else:
    print("using precomputed interpolation weights in: " +"InterpWeights.global_1deg_unstr."+str(N)+".csv")

#InterpWeights=np.load(floutw+".npz")
#weights=InterpWeights["weights"][:,:]
#nodes=InterpWeights["nodes"][:,:]

nodes, weights = ReadWeights(floutw+".csv",)
time=np.asarray(data.variables['time'][:])
nt=0
nt=time.shape[0]
npoints=0
npoints=weights.shape[0]
zetai=np.zeros((nt,npoints))
for k in range(nt):
    zeta=data.variables["zeta"][k,:]
    zetai[k,:]=mshint.InterpolateField2Nodes(nodes,weights, zeta)
    print("on time point "+str(k)+" of "+str(nt))

#np.savez_compressed(floutz+".npz",zetai=zeta) #the .npz files were very large relative text files
WriteNodeByTimeMat(floutz+".txt",zetai)

# Interpolate velocities
"""
data = nc.Dataset("stofs_2d_glo.t00z.fields.cwl.vel.nc","r")
ui=np.zeros((nt,npoints))
vi=np.zeros((nt,npoints))
for k in range(nt):
    u=data.variables["u-vel"][k,:]
    ui[k,:]=mshint.InterpolateField2Nodes(nodes,weights, u)
    v=data.variables["v-vel"][k,:]
    vi[k,:]=mshint.InterpolateField2Nodes(nodes,weights, v)
    print("on time point "+str(k)+" of "+str(nt))
#np.savez_compressed("u_vel.global_1deg_unstr."+str(N)+".npz",ui=ui)
#np.savez_compressed("v_vel.global_1deg_unstr."+str(N)+".npz",vi=vi)
WriteNodeByTimeMat("test.u."+str(N)+".txt",ui)
WriteNodeByTimeMat("test.v."+str(N)+".txt",vi)
"""

