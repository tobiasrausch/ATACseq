#!/bin/bash
#SBATCH -A korbel                   # group to which you belong
#SBATCH -p 1month                   # partition (queue)
#SBATCH -N 1                        # number of nodes
#SBATCH -n 4                        # number of cores
#SBATCH --no-requeue                # never requeue this job
#SBATCH -J atac                     # job name
#SBATCH --mem 8000M                 # memory pool for all cores
#SBATCH -t 2-23:55:00               # time
#SBATCH -o atac.%N.%j.out           # STDOUT
#SBATCH -e atac.%N.%j.err           # STDERR
#SBATCH --mail-type=FAIL,END        # notifications for job done & fail
#SBATCH --mail-user=rausch@embl.de  # send-to address

# Fetch ATAC-Seq script
ATACSCRIPT=${1}
shift

# Run analysis pipeline
${ATACSCRIPT} $@
