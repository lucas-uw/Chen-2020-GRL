#!/bin/csh

#SBATCH -A hyperion
#SBATCH -t 48:00:00

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6

#SBATCH -J AR-p85
#SBATCH -o /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/slurm_jobs.logs/out.AR-SST.idenfity_AR.NARR.p85
#SBATCH -e /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/slurm_jobs.logs/error.AR-SST.identify_AR.NARR.p85


date

set ARmethod='p85'

set year=2003
foreach month(`seq 1 6`)
	srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step3.identify_NARR_AR.argv.py $year $month $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.NARR.$year.$month.txt &
end
wait
foreach month(`seq 7 9`)
	srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step3.identify_NARR_AR.argv.py $year $month $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.NARR.$year.$month.txt &
end
wait

foreach year(`seq 1981 2002`)
	foreach month(`seq 1 6`)
		srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step3.identify_NARR_AR.argv.py $year $month $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.NARR.$year.$month.txt &
	end
	wait

	foreach month(`seq 7 12`)
		srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step3.identify_NARR_AR.argv.py $year $month $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.NARR.$year.$month.txt &
	end
	wait
end

date
exit 0
