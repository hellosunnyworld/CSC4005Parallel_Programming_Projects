#!/bin/bash
#SBATCH --job-name=your_job_name # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # number of processes = 200 
#SBATCH --cpus-per-task=40      # Number of CPU cores allocated to each process
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)
cd ..
for num in 200 400 800 1000 2000
do
    for cores in 1 4 20 40
    do  
    ./openmp $num 100 $cores
    done
done
# cd /nfsmnt/119010355/CSC4005_2022Fall_Demo/project3_template/
# ./openmp 1000 100 4
# ./openmp 1000 100 20
# ./openmp 1000 100 40
# ./openmp 1000 100 80
# ./openmp 1000 100 120
# ./openmp 1000 100 200