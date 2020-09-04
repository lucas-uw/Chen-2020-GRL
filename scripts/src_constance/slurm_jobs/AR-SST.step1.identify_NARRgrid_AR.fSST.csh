#!/bin/csh

#SBATCH -A hyperion
#SBATCH -t 3:00:00

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

#SBATCH -J AR-fSST
#SBATCH -o /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/slurm_jobs.logs/out.AR-SST.idenfity_AR.fSST
#SBATCH -e /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/slurm_jobs.logs/error.AR-SST.identify_AR.fSST


date

set model='fSST'
set ARmethod='abs'

set year=2003
foreach month(`seq 10 12`)
	srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step1.identify_AR.argv.py $model $year $month $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.$model.$year.$month.txt &
end
wait

foreach year(`seq 2004 2015`)
	foreach month(`seq 1 12`)
		srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step1.identify_AR.argv.py $model $year $month $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.$model.$year.$month.txt &
	end
	wait
end

#srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step1.identify_AR.py $model $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.$model.full.txt

date
exit 0
