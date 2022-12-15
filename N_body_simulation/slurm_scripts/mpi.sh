#!/bin/bash
#SBATCH --job-name=your_job_name # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=40                   # number of processes = 40
#SBATCH --cpus-per-task=1      # Number of CPU cores allocated to each process (please use 1 here, in comparison with pthread)
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)
cd ..
ulimit -s 1024000000
for num in 200 400 800 1000 2000
do
    for cores in 1 2 4 10 20 30 40
    do  
    mpirun -np $cores ./mpi $num 100 
    done
done
# cd /nfsmnt/119010355/CSC4005_2022Fall_Demo/project3_template/
# mpirun -np 4 ./mpi 1000 100
# mpirun -np 20 ./mpi 1000 100
# mpirun -np 40 ./mpi 1000 100