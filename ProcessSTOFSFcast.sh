#!/bin/bash
#PBS -N ESMPy
#PBS -j oe
#PBS -S /bin/bash
#PBS -q dev
#PBS -A NWPS-DEV
#PBS -l walltime=00:05:00
#PBS -l select=1:ncpus=1:mem=8G
#PBS -l place=excl
#PBS -l debug=true

module reset
module load PrgEnv-intel/8.5.0
module load intel/19.1.3.304
module load craype/2.7.17
module load cray-mpich/8.1.19
module load hdf5-C/1.14.0
module load netcdf-C/4.9.2
module load esmf-C/8.6.0
module load ve/hafs/2.1

pip list -v


#Script interpolates STOFS velocity and elevation to RWPS mesh for forcing


##date="20260527"
##cycl="00"
date=$1
cycl=$2

mesh=meshes/RWPS.V0a.msh
outdir=rwps.V0a.$date.$cycl

stofsele=stofs.$date.$cycl/stofs_2d_glo.t"$cycl"z.fields.cwl.nc
stofsvel=stofs.$date.$cycl/stofs_2d_glo.t"$cycl"z.fields.cwl.vel.nc
rwpsele=$outdir/rwps.stofs.$date.$cycl.cwl.nc
rwpsvel=$outdir/rwps.stofs.$date.$cycl.cwl.vel.nc



mkdir $outdir

python3 InterpolateSTOFS.py $stofsele $mesh $rwpsele zeta 2
python3 InterpolateSTOFS.py $stofsvel $mesh $rwpsvel u-vel:v-vel 2



#stofs.20260610.00/stofs_2d_glo.t00z.fields.cwl.nc
#stofs.20260610.00/stofs_2d_glo.t00z.fields.cwl.vel.nc
#meshes/RWPS.V0a.msh

