#!/bin/csh

#SBATCH -A hyperion
#SBATCH -t 3:00:00

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

#SBATCH -J stats-full
#SBATCH -o /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/slurm_jobs.logs/out.AR-SST.extract_full_ocean_stats
#SBATCH -e /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/slurm_jobs.logs/error.AR-SST.extract_full_ocean_stats

#module load python/anaconda3.6


date

set ARtag='abs'
set flag_USstate=1
set flag_post_adj=1

srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step7.extract_IVT_SST.full_ocean.py >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.step7.txt


date
exit 0
