teststr=$'hadd\n'
teststr+=beamrootfiles/beam_aida$1_$2to$3.root
teststr+=$'\n'
for i in `seq $2 $3`;
do
    #./runbash.sh $i
    teststr+=beamrootfiles/wasabisorted$i.root
    teststr+=$'\n'
    #echo $i
done 

$teststr
