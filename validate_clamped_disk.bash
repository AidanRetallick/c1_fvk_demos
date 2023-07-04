#! /bin/bash

if [ -e RESLT ]; then
    echo "RESLT already exists; please delete"
    exit
fi
if [ -e RESLT_clamped_2d ]; then
    echo "RESLT_clamped_2d already exists; please delete"
    exit
fi
if [ -e RESLT_clamped_axisym ]; then
    echo "RESLT_clamped_axisym already exists; please delete"
    exit
fi

mkdir RESLT
./circular_disc --use_clamped_bc
mv RESLT RESLT_clamped_2d

mkdir RESLT
./axisym_displ_based_fvk 
mv RESLT RESLT_clamped_axisym
cd RESLT_clamped_axisym
../pad_axisymmetric_to_cartesian.bash soln0.dat 

cd ..
tecplot validate_clamped_disk.lay

   
