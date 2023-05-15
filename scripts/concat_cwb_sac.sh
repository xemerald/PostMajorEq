#!/bin/bash
#
#
if [ $# -ne 1 ]
then
	echo "Usage: concat_cwb_sac.sh <Directory Path>"
	exit 0
fi

echo "Listing all the archived SAC files..."
FILELIST=`cd ${1}; ls | grep .SAC`
for file in ${FILELIST}
do
	if [ -f ${1}/${file} ]
	then
		echo "Finding the SAC files for the same SCNL of ${file}..."
		STATION=`echo ${file} | cut -d. -f8`
		NETWORK=`echo ${file} | cut -d. -f7`
		LOCATION=`echo ${file} | cut -d. -f9`
		CHANNEL=`echo ${file} | cut -d. -f10`
		NSLC=${NETWORK}.${STATION}.${LOCATION}.${CHANNEL}
		SAMELIST=(`cd ${1}; ls | grep ${NSLC}`)
		if [ ${#SAMELIST[@]} -gt 1 ]
		then
			for samefile in ${SAMELIST[@]}
			do
				if [ ${file} != ${samefile} ]
				then
					./concat_sac ${1}/${file} ${1}/${samefile} ${1}
					rm -f ${1}/${samefile}
					#echo ${samefile}
				fi
			done
		fi
	fi
done
#
exit
