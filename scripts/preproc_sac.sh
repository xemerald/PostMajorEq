#!/bin/bash
#
#
if [ $# -ne 2 ]
then
	echo "Usage: concat_cwb_sac.sh <Directory Path> <SCNL Info>"
	exit 0
fi
echo "Creating the folder to store the coverted files..."
OUTPUT_DIR="_preproced"
mkdir -p ${1}/${OUTPUT_DIR}
echo "Listing all the archived SAC files..."
FILELIST=`cd ${1}; ls | grep TW`
for file in ${FILELIST}
do
	if [ -f ${1}/${file} ]
	then
		echo "Finding the SAC files..."
		STATION=`echo ${file} | cut -d. -f1`
		NETWORK=`echo ${file} | cut -d. -f3`
		LOCATION=`echo ${file} | cut -d. -f4`
		CHANNEL=`echo ${file} | cut -d. -f2`
		dirc=`echo ${CHANNEL} | cut -b 3`
		if [ ${dirc} == 'Z' ]
		then
			dirc=8
		elif [ ${dirc} == 'N' ]
		then
			dirc=10
		elif [ ${dirc} == 'E' ]
		then
			dirc=12
		fi
		factor=(`cat ${2} | grep ${STATION} | cut -d' ' -f${dirc}`)
		./preproc_sac -g ${factor} ${1}/${file} ${1}/${OUTPUT_DIR}/${file}
	fi
done
#
exit
