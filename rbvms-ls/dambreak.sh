#!/bin/bash

# Set variables
exe=rbvms-ls
mesh=dambreak.msh

# Clean
rm -rf solution/* output_*.dat

# Create functions
mpicc -shared -o dambreak.so -fPIC dambreak.c

# Create mesh -- if required
if [ ! -f $mesh ]; then
   gmsh -2 dambreak.geo -format msh22 -clscale 0.04
fi

# Actual run 
#valgrind --tool=memcheck 
nohup mpirun -n 4 --oversubscribe $exe \
-l dambreak.so -m $mesh --normal-bdr "1" \
-s 45 -dt 0.005 --dt_vis 0.01 \
--linear-tolerance 1e-4 --linear-itermax 250 \
--newton-tolerance 1e-3 --newton-itermax 10 > log &

tail -111f log
