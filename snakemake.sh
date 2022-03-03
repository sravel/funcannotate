#!/bin/bash 
#
#SBATCH -J functannotate
#SBATCH -o functannotate."%j".out
#SBATCH -e functannotate."%j".err 

#SBATCH -p agap_long

module purge
module load snakemake/5.13.0 
 
mkdir -p logs/slurm/
 


# snakemake  --cluster-config cluster_config.yml  --configfile config.yaml --jobs 200 --rerun-incomplete
#snakemake  --cluster-config cluster_config.yml  --configfile config.yaml --jobs 200 --printshellcmds --dryrun --use-envmodules  --cluster "sbatch  -p {cluster.partition} --ntasks {cluster.ntasks} --mem-per-cpu={cluster.mem-per-cpu} -c {cluster.cpus-per-task} -e {cluster.error} -o {cluster.output}"
# snakemake  --cluster-config cluster_config.yml  --configfile config.yaml --jobs 200 --unlock  --use-envmodules
snakemake  --cluster-config cluster_config.yml  --configfile config.yaml --jobs 2 --dag  --use-envmodules | dot -Tpdf > dag.pdf
#snakemake --cluster-config cluster_config.yml  --configfile config.yaml --jobs 200  --use-envmodules  --cluster "sbatch  -p {cluster.partition} --ntasks {cluster.ntasks} --mem-per-cpu={cluster.memory} -c {cluster.cpus-per-task} -e {cluster.error} -o {cluster.output}"