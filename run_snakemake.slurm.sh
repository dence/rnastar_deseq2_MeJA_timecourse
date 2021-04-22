#!/bin/sh
#SBATCH --job-name=timecourse # Job name
#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mallory.morgan@ufl.edu  # Where to send mail
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=15gb                  # Total memory limit
#SBATCH --time=192:00:00              # Time limit hrs:min:sec
#SBATCH --output=master_%j.out     # Standard output and error log
#SBATCH --qos=peter
#SBATCH --account=peter


unset TMPDIR
module load python3

snakemake --rerun-incomplete --configfile config.yaml --snakefile Snakefile -c 20 --jobs 10 --directory . --cluster-config hipergator.cluster.json --cluster "sbatch --qos={cluster.qos} -c {cluster.c} -n {cluster.N} --mail-type=FAIL --mail-user=mallory.morgan@ufl.edu -t {cluster.time} --mem={cluster.mem} -J "timecourse" -o timecourse_%j.out -D /blue/peter/mallory.morgan/rnastar_deseq2_MeJA_timecourse"
#snakemake --configfile config/config.yaml --snakefile ./workflow/Snakefile -c 1 --jobs 1 --directory . --cluster-config ../hipergator.cluster.json --cluster "sbatch --qos={cluster.qos} -p {cluster.partition} -c {cluster.c} -n {cluster.N} --mail-type=FAIL --mail-user=d.ence@ufl.edu -t {cluster.time} --mem={cluster.mem} -J "fr1_align" -o fr1_align_%j.out -D /home/d.ence/projects/pinus_taeda_L/Fr1_project/test_pipelines/snakemake_pipelines/pipelines_for_hipergator/fr1_project_snakefiles"
