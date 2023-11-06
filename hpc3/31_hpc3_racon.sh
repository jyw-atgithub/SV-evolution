#!/bin/bash

#SBATCH --job-name=polish    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard        ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=36   ## number of cores the job needs
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU

