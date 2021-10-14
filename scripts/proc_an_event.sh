#!/bin/bash
#
#
WORKING_DIR=/root/working
TMP_DIR=${WORKING_DIR}/tmp
SRC_DIR=${WORKING_DIR}/data_src
RES_DIR=${WORKING_DIR}/data_res

#
if [ $# -ne 9 ]
then
	echo "Usage: proc_an_event.sh <Year> <Month> <Day> <Hour> <Minute> <Second> <Lat.> <Lon.> <Depth>"
	exit 0
fi

#
EV_YMD=`printf "%d%02d%02d" $1 $2 $3`
EV_HM=`printf "%02d%02d %02d%02d %02d%02d" $4 $[$5 - 3] $4 $[$5 - 2] $4 $[$5 - 1]`
FILEPATH=${SRC_DIR}/EVENTS_NTU/${EV_YMD:0:6}
FILENAME=''
#
for _EV_HM in ${EV_HM}
do
	FILENAME=`cd ${FILEPATH}; ls ${EV_YMD}_${_EV_HM}*_MAN.tar.bz2 2> /dev/null`
	FILENAME=${FILENAME//[[:space:]]/}
	if [ ${FILENAME} ] && [ -e "${FILEPATH}/${FILENAME}" ]
	then
		echo "Found the archived SAC file: '${FILENAME}'!"
		break
	fi
done
#
if [ ! ${FILENAME} ]
then
	echo "Can't find the SAC file for this event! Just exit."
	exit 0
fi

#
_SEC=`echo $6 | cut -d. -f1`
EVID=`printf "%s%02d%02d%02d" ${EV_YMD} $4 $5 ${_SEC}`
EQINFO=${TMP_DIR}/eqinfo
STALIST=${WORKING_DIR}/stalist
RESULT=${RES_DIR}/${EVID}_res.txt
SAC_DIR=${TMP_DIR}/${FILENAME:0:19}
#
echo "Start to process Eq. ${EVID}..."
#
echo "Moving the archived file to local..."
cp ${FILEPATH}/${FILENAME} ${TMP_DIR}
echo "Decompressing the archived SAC files..."
tar -C ${TMP_DIR} -vjxf ${TMP_DIR}/${FILENAME}
rm -f ${TMP_DIR}/${FILENAME}

#
echo "$1 $2 $3 $4 $5 $6 $7 $8 $9" > ${EQINFO}
postmajor ${EQINFO} ${STALIST} ${SAC_DIR} > ${RESULT}
#
rm -rf ${SAC_DIR}
rm -f ${EQINFO}
#
echo "Event processing complete!!"
