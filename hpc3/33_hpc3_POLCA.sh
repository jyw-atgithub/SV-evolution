#!/bin/bash

#SBATCH --job-name=POLCA    ## Name of the job.
#SBATCH -A jje_lab       ## account to charge
#SBATCH -p standard       ## partition/queue name
#SBATCH --array=1      ## number of tasks to launch (wc -l prefixes.txt)
#SBATCH --cpus-per-task=16   ## number of cores the job needs-x map-ont
#SBATCH --mem-per-cpu=6G     # requesting memory per CPU