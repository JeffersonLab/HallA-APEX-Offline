#!/bin/bash

rm *SieveCut.cut
rm *SieveCut.ncut
#LIST="4647 4648 4650"
LIST="4647"
N=2000000

for i in $LIST
do
    echo "$i"
    if [[ "$i" == "4647" ]];then
	j=1
    fi
    if [[ "$i" == "4648" ]];then
	j=0
    fi
    if [[ "$i" == "4650" ]];then
	j=2
    fi
    for k in {0..26}
    do
    f=apex_$i.root.SieveCut.${j}_${k}.cut
    echo $f
    while read l
    do
	if [[ "$l" == "R.tr.n > 1000" ]];then
	    echo -e "0" >> apex_$i.root.SieveCut.ncut 
	else
	    echo -e "$N" >> apex_$i.root.SieveCut.ncut
	fi
	echo -e "$l" >> apex_$i.root.SieveCut.cut
    done < $f
    done
done
