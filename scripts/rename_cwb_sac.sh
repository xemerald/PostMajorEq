#!/bin/bash
#
#
if [ $# -ne 1 ]
then
	echo "Usage: rename_cwb_sac.sh <Directory Path>"
	exit 0
fi

echo "Listing all the archived SAC files..."
FILELIST=`cd ${1}; ls | grep .SAC`
for file in ${FILELIST}
do
	echo "Renaming the SAC files ${file}..."
	STATION=`echo ${file} | cut -d. -f8`
	NETWORK=`echo ${file} | cut -d. -f7`
	LOCATION=`echo ${file} | cut -d. -f9`
	CHANNEL=`echo ${file} | cut -d. -f10`
	mv ${1}/${file} ${1}/${STATION}.${CHANNEL}.${NETWORK}.${LOCATION}
done
#
echo "Listing all the new archived SAC files..."
cd ${1}; echo "# SAC files list" > saclist; ls *.*.*.* >> saclist; cd -
exit
