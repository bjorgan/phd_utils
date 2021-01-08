#/bin/bash

file=$1
num_times=$2

tmpfile=$(mktemp /tmp/dup-spectrum-wlens.XXXXXX)

#output wavelength column to temporary file
awk '{ print $1 " "  }' $file > $tmpfile

#print the rest of the columns the number of times as specified and combine with wavelengths
awk '{ print $2 " "  }' $file | perl -lne "print \$_ x $num_times" | paste -d' ' $tmpfile -

rm $tmpfile

