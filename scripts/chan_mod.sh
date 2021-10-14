#!/bin/bash
#
#
WORKING_DIR=/root/working
TMP_DIR=${WORKING_DIR}/tmp
SRC_DIR=${WORKING_DIR}/data_src
RES_DIR=${WORKING_DIR}/data_res

#
if [ $# -ne 6 ]
then
	echo "Usage: chan_mod.sh <Year> <Month> <Day> <Hour> <Minute> <Second>"
	exit 0
fi

#
EV_YMD=`printf "%d%02d%02d" $1 $2 $3`
EV_HM=`printf "%02d%02d %02d%02d %02d%02d" $4 $[$5 - 3] $4 $[$5 - 2] $4 $[$5 - 1]`
FILEPATH=${SRC_DIR}/EVENTS_NTU/${EV_YMD:0:6}
SAC_DIR=''
#
for _EV_HM in ${EV_HM}
do
	SAC_DIR=`cd ${FILEPATH}; ls | grep ${EV_YMD}_${_EV_HM}*_MAN`
	SAC_DIR=${SAC_DIR//[[:space:]]/}
	if [ ${SAC_DIR} ] && [ -e "${FILEPATH}/${SAC_DIR}" ]
	then
		echo "Found the archived SAC path: '${SAC_DIR}'!"
		break
	fi
done
#
if [ ! ${SAC_DIR} ]
then
	echo "Can't find the SAC path for this event! Just exit."
	exit 0
fi

#
_SEC=`echo $6 | cut -d. -f1`
EVID=`printf "%s%02d%02d%02d" ${EV_YMD} $4 $5 ${_SEC}`
#
echo "Start to process Eq. ${EVID}..."
#
echo "Moving the archived files to local..."
cp -R ${FILEPATH}/${SAC_DIR} ${TMP_DIR}
#echo "Decompressing the archived SAC files..."
#tar -C ${TMP_DIR} -vjxf ${TMP_DIR}/${SAC_DIR}
#rm -f ${TMP_DIR}/${SAC_DIR}
echo "Listing all the P-Alert archived SAC files..."
FILELIST=`cd ${TMP_DIR}/${SAC_DIR}; ls | grep "[A-Z][0-9]\{2,3\}[A-F]\{0,1\}.HH[ZNE].TW.--"`
for file in ${FILELIST}
do
	NEW_CHAN=`echo ${file} | cut -d. -f2 | sed 's/HH/HL/'`
	chan_mod ${TMP_DIR}/${SAC_DIR}/${file} ${NEW_CHAN} ${TMP_DIR}/${SAC_DIR}
	rm -f ${TMP_DIR}/${SAC_DIR}/${file}
done
#
echo "Compressing all archived SAC files..."
cd ${TMP_DIR}; tar -vjcf ${SAC_DIR}.tar.bz2 ${SAC_DIR}; cd -
#
echo "Deleting the temporary folder..."
rm -rf ${TMP_DIR}/${SAC_DIR}
#
#
echo "Event processing complete!!"
