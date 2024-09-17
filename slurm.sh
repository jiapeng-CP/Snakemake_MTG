#!/bin/bash

#SBATCH --nodes=2
#SBATCH --mail-user=jiapeng_chen@colpal.com
#SBATCH --nodelist=alustt132,alustt133
#SBATCH --cpus-per-task=96
#SBATCH --distribution=cyclic:cyclic    # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --mail-type=BEGIN,END,FAIL      # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --partition=defq
#SBATCH --time=180-00:00:00 # 6 months

Username=`whoami`

FreeMem=`scontrol show nodes | grep Free | awk '{print $3}' | awk -F'=' '{print $2}' | awk '{ sum+=$1 } END {print sum/1024}'`

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo "Free Memory on Cluster         = $FreeMem Gb"

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

mkdir -p /home/$Username/slurm/

conda activate /home/jiapengc/bin/biobakery4


/home/jiapengc/mambaforge/envs/snakemake/bin/snakemake --core 96
