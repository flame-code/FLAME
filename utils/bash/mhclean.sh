#!/bin/sh

echo "********** backup/prepare for restarting minhopp **********"

ZIP=gzip
scriptname=`basename "$0"`

if [ -d $1 ]; then
    echo "ERROR: directory $1 exists, $scriptname stops."
    exit
else
    mkdir $1 -m 0750
fi

filename=posinp.acf
if [ -s "$filename" ]; then
    echo "mv $filename $1 and $ZIP $1/$filename"
    mv "$filename" $1 ; "$ZIP" $1/"$filename"
fi

filename=posout.acf
if [ -s "$filename" ]; then
    echo "cp $filename $1 and $ZIP $1/$filename"
    cp "$filename" $1 ; "$ZIP" $1/"$filename"
fi

filename=poslow.acf
if [ -s "$filename" ]; then
    echo "cp $filename $1 and $ZIP $1/$filename"
    cp "$filename" $1 ; "$ZIP" $1/"$filename"
fi

filename=earr.dat
if [ -s "$filename" ]; then
    echo "cp $filename $1 and $ZIP $1/$filename"
    cp "$filename" $1 ; "$ZIP" $1/"$filename"
fi

filename=input.minhopp
if [ -s "$filename" ]; then
    echo "cp $filename $1 and $ZIP $1/$filename"
    cp "$filename" $1 ; "$ZIP" $1/"$filename"
fi

echo "mv posout.acf posinp.acf"
mv posout.acf posinp.acf

echo "********** backup/prepare done **********"
