teststr=$'hadd -k\n'
teststr+=beamrootfiles/beam_aida$1_$2to$3.root
teststr+=$'\n'
for i in `seq $2 $3`;
do
    ./runbash.sh $i
    scp -i mykey beamdata/beamrun0$i.root daq@10.32.6.194:/home/daq/Desktop/dummy/briken1/wasabi
    teststr+=beamrootfiles/wasabisorted$i.root
    teststr+=$'\n'
    #echo $i
done 

#$teststr
