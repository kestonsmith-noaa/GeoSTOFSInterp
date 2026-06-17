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


##Script interpolates STOFS velocity and elevation to RWPS mesh for forcing


##date="20260527"
##cycl="00"
##mesh="../meshes/RWPS.V0a.small.msh"

date=$1
cycl=$2
mesh=$3

meshname="${mesh##*/}"
meshname="${meshname: 0: -4}"

outdir=$meshname.$date.$cycl

stofsele="stofs.$date.$cycl/stofs_2d_glo.t"$cycl"z.fields.cwl.nc"
stofsvel="stofs.$date.$cycl/stofs_2d_glo.t"$cycl"z.fields.cwl.vel.nc"
rwpsele="$outdir/stofs.$date.$cycl.cwl.nc"
rwpsvel="$outdir/stofs.$date.$cycl.cwl.vel.nc"

mkdir $outdir

echo "stofs  outdir=$outdir"
echo "stofs    mesh=$mesh"
echo "stofs rwpsele=$rwpsele"
echo "stofs rwpsvel=$rwpsvel"

python InterpolateSTOFS.py $stofsele $mesh $rwpsele zeta -2
python InterpolateSTOFS.py $stofsvel $mesh $rwpsvel u-vel:v-vel -2

