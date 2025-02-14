#!/bin/bash
first=$1
second=$2
while [ "$first" -le "$second" ]
do
    ./fy3g_mwri_retrieve_emiss_swath.exe $first		
    echo $first finished
    let first=`date -d "-1 days ago ${first}" +%Y%m%d`
done
