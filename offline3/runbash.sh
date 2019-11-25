inpfile=$(printf "beamrun%05d" $1)
./convert.sh beamdata/$inpfile.evt beamdata/$inpfile.root
./dosort.sh beamdata/$inpfile.root beamrootfiles/wasabisorted$1.root