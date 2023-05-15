#!/bin/bash
#
#
if [ $# -ne 1 ]
then
	echo "Usage: convert_cwb_stainfo.sh <Input station info>"
	exit 0
fi

datas=`cat ${1}`
#reg='^[+-]?[0-9]+([.][0-9]+)?$'
IFS=$'\n'
sta=" "
chan_z=" "
chan_n=" "
chan_e=" "
result=" "
fin=0
for data in ${datas}
do
#echo "Working on ${_file}..."
	#echo ${data}
	_sta=`echo ${data} | cut -d' ' -f1`
	if [ ${_sta} != ${sta} ]
	then
		sta=${_sta}
		net="TW"
		loc=`echo ${data} | cut -d' ' -f3`
		coor=`echo ${data} | cut -d' ' -f4-6`
		result="${sta} ${net} ${loc} ${coor}"
		fin=0
	fi
#
	_chan=`echo ${data} | cut -d' ' -f2`
	_fac=`echo ${data} | cut -d' ' -f11`
	_fac=`printf "%.10f" ${_fac}`
	_fac=`echo "${_fac} * 100.0" | bc`
	dirc=`echo ${_chan} | cut -b 3`
	if [ ${dirc} == 'Z' ]
	then
		chan_z="${_chan} ${_fac}"
		fin=$((${fin} + 1))
	elif [ ${dirc} == 'N' ]
	then
		chan_n="${_chan} ${_fac}"
		fin=$((${fin} + 1))
	elif [ ${dirc} == 'E' ]
	then
		chan_e="${_chan} ${_fac}"
		fin=$((${fin} + 1))
	fi
#
	if [ ${fin} -eq 3 ]
	then
		echo "${result} ${chan_z} ${chan_n} ${chan_e}"
	fi
done
exit
