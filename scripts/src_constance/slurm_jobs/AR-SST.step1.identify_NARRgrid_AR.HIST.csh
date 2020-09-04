#!/bin/csh

#SBATCH -A hyperion
#SBATCH -t 1:00:00

#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5

#SBATCH -J AR-HIST
#SBATCH -o /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/slurm_jobs.logs/out.AR-SST.idenfity_AR.HIST
#SBATCH -e /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/slurm_jobs.logs/error.AR-SST.identify_AR.HIST


date

set model='HIST'
set ARmethod='abs'

#set year=2003
#foreach month(`seq 10 12`)
#	srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step1.identify_AR.argv.py $model $year $month $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.$model.$year.$month.txt &
#end
#wait

#foreach year(`seq 2004 2015`)
#	foreach month(`seq 1 12`)
#		srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step1.identify_AR.argv.py $model $year $month $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.$model.$year.$month.txt &
#	end
#	wait
#end


# extended

#foreach year(`seq 1981 2002`)
#	foreach month(`seq 1 12`)
#		srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step1.identify_AR.argv.py $model $year $month $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.$model.$year.$month.txt &
#	end
#	wait
#end

#set year=2003
#foreach month(`seq 1 9`)
#	srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step1.identify_AR.argv.py $model $year $month $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.$model.$year.$month.txt &
#end
#wait

srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step1.identify_AR.argv.ext.py $model 1981 1 $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.$model.1981.1.txt &
srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step1.identify_AR.argv.ext.py $model 1989 11 $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.$model.1989.11.txt &
srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step1.identify_AR.argv.ext.py $model 1995 4 $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.$model.1985.4.txt &
srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step1.identify_AR.argv.ext.py $model 1997 12 $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.$model.1997.12.txt &
srun -n1 /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/step1.identify_AR.argv.ext.py $model 1999 3 $ARmethod >& /pic/projects/hyperion/chen423/tools/paper_tools/AR-SST/logs/log.identify_AR.$model.1999.3.txt &
wait

date
exit 0
