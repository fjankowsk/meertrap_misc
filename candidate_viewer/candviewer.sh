#!/bin/bash
#
#   2021 Fabian Jankowski
#   A simple viewer for single-pulse candidate plots.
#

files=(`ls -1 *.jpg | awk -F '_' '{print $0, $3}' | sort -rnk 2 | awk '{print $1}'`)
nfiles=${#files[@]}

date=`date --utc`
statefile=~/fj_viewer.state

echo "Total number of files: ${nfiles}"
echo "Current date: ${date}"

echo "# ${date} -- new run" >> ${statefile}

i=0

for file in ${files[@]}
do
    percent=$(( 100 * ${i} / ${nfiles} ))
    echo "Viewing: $i, ${file} (${percent} %)"
    display ${file}
    echo ${file} >> ${statefile}
    i=$((i + 1))
done

echo "All candidates viewed."
