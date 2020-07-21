FOILS="3 4 5"

for foil in $FOILS
do
    for col in {0..26}
    do
	FILE=apex_4653_opt_5th_xfp_full_V_wires.root.SieveCut.${foil}_$col.cut
	if test -f "$FILE"; then
	    continue
	fi
	for row in {1..17}
	do
	    echo "R.tr.n > 1000" >> $FILE
	done	    
    done
done
    
