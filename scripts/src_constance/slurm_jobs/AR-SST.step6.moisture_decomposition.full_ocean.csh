#!/bin/csh

#SBATCH -A hyperion
#SBATCH -t 3:00:00

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12

#SBATCH -J moisture
#SBATCH -o /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/slurm_jobs.logs/out.AR-SST.step6.moisture
#SBATCH -e /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/slurm_jobs.logs/error.AR-SST.step6.moisture

date

set year=2003
foreach month(`seq 10 12`)
	srun -n1 /people/chen423/sw/anaconda3/bin/python /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step6.moisture_decomposition.no_AR_boundary.py $year $month >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.step6.full_ocean.HIST.$year.$month.txt &
end
wait

foreach year(`seq 2004 2014`)
	foreach month(`seq 1 12`)
		srun -n1 /people/chen423/sw/anaconda3/bin/python /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step6.moisture_decomposition.no_AR_boundary.py $year $month >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.step6.full_ocean.HIST.$year.$month.txt &
	end
	wait
end
		
set year=2015
foreach month(`seq 1 9`)
	srun -n1 /people/chen423/sw/anaconda3/bin/python /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step6.moisture_decomposition.no_AR_boundary.py $year $month >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.step6.full_ocean.HIST.$year.$month.txt &
end
wait


date
exit 0
