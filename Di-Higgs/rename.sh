#!/bin/bash

string1=$1
string2=$2

obj_file=`find -maxdepth 1 -name "*"$string1"*"`
for files in $obj_file
do
    with_label=$(echo $files | grep $string2)
    if [ 0 != $? ]
    then
	    filename=${files##./}
        #echo $filename
        mv $filename $string2$filename
    fi
done