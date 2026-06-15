#/bin/bash

date=$1
cycl=$2



# /lfs/h1/ops/prod/com/stofs/v2.1/stofs_2d_glo.20260608/stofs_2d_glo.t06z.fields.cwl.vel.nc
# /lfs/h1/ops/prod/com/stofs/v2.1/stofs_2d_glo.20260608/stofs_2d_glo.t00z.fields.cwl.nc

mkdir stofs.$date.$cycl

cp /lfs/h1/ops/prod/com/stofs/v2.1/stofs_2d_glo.$date/stofs_2d_glo.t"$cycl"z.fields.cwl.nc  stofs.$date.$cycl/
cp /lfs/h1/ops/prod/com/stofs/v2.1/stofs_2d_glo.$date/stofs_2d_glo.t"$cycl"z.fields.cwl.vel.nc  stofs.$date.$cycl/

