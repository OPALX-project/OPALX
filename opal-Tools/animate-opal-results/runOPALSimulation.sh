#!/bin/bash


if [[ "$1" == "" || ! -f $1 ]]; then
    echo "usage: runRingCyclotron.sh input"
    echo "       where input is the OPAL input file without the .in"
    exit 0
fi

INPUT=`echo ${1##*/}`;
PATHTOINPUT=`echo $1 |sed "s/$INPUT//"`;
BASENAME=`echo ${INPUT} | sed 's/\.in//'`

cd ./${PATHTOINPUT};
pwd
rm -rf ${BASENAME}.mpg &> /dev/null

# run the opal
mpirun -np 4 opal ${INPUT} --commlib mpi &
MPIRUN_PID=$!
echo ${MPIRUN_PID} > mpirun_pid

# wait a few seconds, check if opal still runs and if so 
# then start AnimateSimulation
sleep 8
if ! kill -0 ${MPIRUN_PID} &> /dev/null; then
    exit
fi

AnimateSimulation ${BASENAME}.h5 &> animate_simulation_output.txt &
AS_PID=$!
echo ${AS_PID} > as_pid

# wait a few seconds, then start mplayer
sleep 8
while [ ! -f ${BASENAME}.mpg ]; do
    sleep 1
done

mplayer -fs -loop 0 ${BASENAME}.mpg &> mplayer_output &

# check if mpirun has finished
while true; do
    if ! kill -0 ${MPIRUN_PID} &> /dev/null; then
	break;
    fi

    sleep 1;
done

#wait until AnimateSimulation has finished, then kill it
TAIL=`tail -1 animate_simulation_output.txt`;
LIAT="Reloading file ${BASENAME}.h5"
while [[ ! "$TAIL" =~ "$LIAT" ]]; do
    sleep 2;
    TAIL=`tail -1 animate_simulation_output.txt`;
done

kill ${AS_PID};
