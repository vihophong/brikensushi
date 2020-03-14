datadir=beam
outputdir=outrootfiles
cat /dev/null >list/list_wasabi$1to$2.txt
for i in `seq $1 $2`;
do
    printf -v input "%05d" $i
    ls -tr $datadir/beamrun$input.root >> list/list_wasabi$1to$2.txt
done

root -b -q 'dosort.C("list/list_wasabi'$1'to'$2'.txt","'$outputdir'/beamrun'$3'_wrun'$1'to'$2'.root")'
