#(run from current directory)
#$ -cwd
#(stderr and stdout are merged together to stdout)
#$ -j y
#(use bash as shell)
#$ -S /bin/bash
#(Tell the SGE that this is an array job, with "tasks" to be numbered 1 to 10000)
#$ -t 201-281
# -l hio=1
#User job limit -> 1000/nJobs = ujl310
#$ -l ujl=1
#(the cpu time for this job)
#$ -l h_cpu=100:59:59
#$-l os=slc5 
#(the maximum memory usage of this job)
#$-l h_vmem=9000M
#$ -V

#(setup the environment)
. /local/AZ/envfor.sh root 5.34.10;

#$ -i  ./model_2.txt
#$ -o ./model_2/out_$TASK_ID.log

./ampFit $SGE_TASK_ID