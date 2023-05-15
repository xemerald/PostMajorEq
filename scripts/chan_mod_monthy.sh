#!/bin/bash
#
#
WORKING_DIR=/root/working
TMP_DIR=${WORKING_DIR}/tmp
SRC_DIR=${WORKING_DIR}/data_src
RES_DIR=${WORKING_DIR}/data_res

#
if [ $# -ne 2 ]
then
	echo "Usage: chan_mod_monthy.sh <Year> <Month>"
	exit 0
fi

#
EV_YM=`printf "%d%02d" $1 $2`
FILEPATH=${SRC_DIR}/unarranged/${EV_YM}
SAC_DIRS=`cd ${FILEPATH}; ls | grep _MAN`
#
for dir in ${SAC_DIRS}
do
	if [ ${dir} ] && [ -d "${FILEPATH}/${dir}" ]
	then
		echo "Found the archived SAC path: '${dir}'!"
		echo "Start to process the SAC path: ${dir}..."
		#
		echo "Moving the archived files to local..."
		cp -R ${FILEPATH}/${dir} ${TMP_DIR}
		#echo "Decompressing the archived SAC files..."
		#tar -C ${TMP_DIR} -vjxf ${TMP_DIR}/${SAC_DIR}
		#rm -f ${TMP_DIR}/${SAC_DIR}
		chan_mod_dir.sh ${TMP_DIR}/${dir}
		#
		echo "Compressing all archived SAC files..."
		cd ${TMP_DIR}; tar -vjcf ${dir}.tar.bz2 ${dir}; cd -
		mv ${TMP_DIR}/${dir}.tar.bz2 ${RES_DIR}
		#
		echo "Deleting the temporary folder..."
		rm -rf ${TMP_DIR}/${dir}
	fi
done
#
echo "All the processing complete!!"
exit
