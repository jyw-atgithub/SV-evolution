#!/bin/bash

#SBATCH --job-name=np    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p highmem        ## partition/queue name
#SBATCH --array=18      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=24   ## number of cores the job needs
#SBATCH --mem-per-cpu=10G     # requesting memory per CPU