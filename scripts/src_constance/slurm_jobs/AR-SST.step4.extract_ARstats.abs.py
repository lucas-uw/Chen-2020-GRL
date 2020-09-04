#!/bin/csh

#SBATCH -A hyperion
#SBATCH -t 3:00:00

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

#SBATCH -J stats-abs
#SBATCH -o /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/slurm_jobs.logs/out.AR-SST.extract_ARstats_abs
#SBATCH -e /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/slurm_jobs.logs/error.AR-SST.extract_ARstats_abs

module load python/anaconda3.6


date

set ARtag='abs'
set flag_USstate=1
set flag_post_adj=1

srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step4.extract_ARstats.py $ARtag 100 $flag_USstate $flag_post_adj 500 >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.step4.$ARtag.100.$flag_USstate.$flag_post_adj.500.txt &

srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step4.extract_ARstats.py $ARtag 100 $flag_USstate $flag_post_adj 1000 >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.step4.$ARtag.100.$flag_USstate.$flag_post_adj.1000.txt &

srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step4.extract_ARstats.py $ARtag 100 $flag_USstate $flag_post_adj 2000 >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.step4.$ARtag.100.$flag_USstate.$flag_post_adj.2000.txt &

srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step4.extract_ARstats.py $ARtag 500 $flag_USstate $flag_post_adj 500 >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.step4.$ARtag.500.$flag_USstate.$flag_post_adj.500.txt &

srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step4.extract_ARstats.py $ARtag 500 $flag_USstate $flag_post_adj 1000 >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.step4.$ARtag.500.$flag_USstate.$flag_post_adj.1000.txt &

srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step4.extract_ARstats.py $ARtag 500 $flag_USstate $flag_post_adj 2000 >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.step4.$ARtag.500.$flag_USstate.$flag_post_adj.2000.txt &

srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step4.extract_ARstats.py $ARtag 1000 $flag_USstate $flag_post_adj 1000 >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.step4.$ARtag.1000.$flag_USstate.$flag_post_adj.1000.txt &

srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step4.extract_ARstats.py $ARtag 1000 $flag_USstate $flag_post_adj 2000 >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.step4.$ARtag.1000.$flag_USstate.$flag_post_adj.2000.txt &
wait

date
exit 0
