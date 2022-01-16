#!/bin/bash
while IFS=$'\t' read fileName fileIn fileOut
do
	Rscript processSyntheticData_sawtooth_limoRhyde.R $fileName $fileIn $fileOut
	# print out the job id for reference later
	echo "JobID = ${JOB} for parameters $fileName submitted on `date`"
done < masterFilePath_limoRhyde_sawtooth.txt
exit

# make this file executable and then run from the command line
# chmod u+x run_limoRhyde_local.sh
# ./run_limoRhyde_local.sh
