#!/bin/bash
#
#
if [ $# -ne 1 ]
then
	echo "Usage: chan_mod_dir.sh <Directory Path>"
	exit 0
fi

echo "Listing all the P-Alert archived SAC files..."
FILELIST=`cd ${1}; ls | grep "^[DHLPSTW][0-9]\{2,3\}[A-F]\{0,1\}\.HH[ZNE]\.TW\.--"`
for file in ${FILELIST}
do
	NEW_CHAN=`echo ${file} | cut -d. -f2 | sed 's/HH/HL/'`
	chan_mod ${1}/${file} ${NEW_CHAN} ${1}
	rm -f ${1}/${file}
done
#
echo "Listing all the new archived SAC files..."
cd ${1}; echo "# SAC files list" > saclist; ls *.*.*.* >> saclist; cd -
exit
