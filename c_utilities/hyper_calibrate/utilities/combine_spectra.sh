#/bin/bash

file_1=$1
file_2=$2

if [ -z $file_1 -o -z $file_2 ]; then
	echo "Not enough input arguments"
	exit
fi

#strip wavelength column
awk -v f=1 -v t=1 '{for(i=1;i<=NF;i++)if(i>=f&&i<=t)continue;else printf("%s%s",$i,(i!=NF)?OFS:ORS)}' $2 | paste -d' ' $1 -
