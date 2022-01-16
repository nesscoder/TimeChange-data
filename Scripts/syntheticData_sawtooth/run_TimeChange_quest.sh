#!/bin/bash

while IFS=$'\t' read fileName fileIn fileOut
do
	here=`pwd`
	JOB=`sbatch << EOJ
#!/bin/bash
#SBATCH -A p30673
#SBATCH --partition=normal
#SBATCH --time=48:00:00
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=elanness-cohn2017@u.northwestern.edu
#SBATCH -J $fileName
#SBATCH --nodes=1
#SBATCH -n 7
#SBATCH --mem=5G

# unload modules that may have been loaded when job was submitted
module purge all

# load the module you want to use 
module load R/4.0.3

# By default all file paths are relative to the directory where you submitted the job.
cd /home/emn6548/TimeChange-data/Scripts/syntheticData_sawtooth/

Rscript processSyntheticData_sawtooth_quest.R $fileName $fileIn $fileOut
EOJ
`
# print out the job id for reference later
echo "JobID = ${JOB} for parameters $fileName submitted on `date`"
done < masterFilePath_timeChange_sawtooth.txt
exit

# make this file executable and then run from the command line
# chmod u+x run_TimeChange_quest.sh
# ./run_TimeChange_quest.sh
