# RBVMS on Delft Blue

Delft Blue is a TU Delft resource, see
[Delft Blue](https://www.tudelft.nl/dhpc/system)
for general information.


## Account

To be able to run RBVMS on Delft Blue, you will need a
Delf blue [account](https://doc.dhpc.tudelft.nl/delftblue/Accounting-and-shares/)

Furhermore, it helps to be familiar with Linux and the 
specific policies of the cluster.
See the
[Delft Blue manuals](https://doc.dhpc.tudelft.nl/delftblue)
for more information.

## Install

The RBVMS code can be downloaded using the following command:

```
git clone git@github.com:IdoAkkerman/rbvms.git

```
or
```
git clone https://github.com/IdoAkkerman/rbvms.git

```
This will create a directory `rbvms` with the code in it.
 Before we can compile the code we need to make sure the correct modules are loaded.
This can be done using the following command:

```
source rbvms/config/delftblue-modules.sh
```
After this, the code can be compiled using the same commands as in the local case:
```
cd rbvms
mkdir build
cd build
cmake ..
make -j 16
```


## Run

On the cluster only small admin jobs, such as compilation and mesh generation can be done on the login node. The actual computations need to be put in the **queue**, so the cluster can schedule your job and run it on the assigned compute nodes.

To submit a job you need to go to the appropriate directory and run:
```
sbatch  batch.sh
```
where batch.sh is the batch file that includes all the run specifics. An example can be found at
[delft-blue-modules.sh](config/delft-blue-modules.sh).

It will look something like:
```
#!/bin/sh
#
#SBATCH --job-name="rbvms-von-karman"
#SBATCH --time=01:00:00
#SBATCH --ntasks=16
#SBATCH --partition=compute
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G

DIR=$HOME/rbvms
EXE=$DIR/build/bin/rbvms

source $DIR/config/delft-blue-modules.sh

srun $EXE \
--mesh von-karman-nurbs.mesh --order 2 --refine 4 \
--lib libvkvs.so --strong-bdr "1 2" --outflow-bdr 4  --weak-bdr 3 \
--ode-solver 45 --dt 0.1 --t-final 1000 \
--dt_vis 0.5 --restart-interval 10 \
--linear-tolerance 1e-3 --newton-itermax 20
````

The `#SBATCH` lines provide information to the batch system. YOu should not be affraid to change these, or add additonal information.

The `source` line makes sure the appropraite libraries are loaded. The `srun` command is the rbvms call with the associated cli options. 
Please note that the options
```
--mesh von-karman-nurbs.mesh
--lib libvkvs.so 
```
assume the library and mesh file are in place.


## Managing the queue

To see the current queue:
```
squeue
```

The list might be to long to be useful. You can also run for instance:
```
queue --me
```
which will only show your own jobs in the queue.

If a job needs to be cancelled you can do that using:

```
scancel  #id
```
where `#id` is the job-id. So:
```
scancel 5531858
```

TO SAVE COMPUTATIONAL RESOURCES: PLEASE CANCEL JOBS IF THEY ARE NO LONGER NEEDED!!!!!

## Additional info

This document is not intended to be complete.
For more information see for instance:
[Delft Blue manuals](https://doc.dhpc.tudelft.nl/delftblue)

Note that the Delft Blue is uses slurm for the queueing system.