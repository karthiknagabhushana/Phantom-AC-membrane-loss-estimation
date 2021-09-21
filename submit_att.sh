#!/bin/bash
source /etc/bashrc

#Setup files
INPUTFILE=Txdc_10M_F4_3D_att.h5
OUTPUTFILE=Txdc_10M_F4_3D_att_out.h5
CHECKPOINT_FILE=Txdc_10M_F4_3D_att_bkp.h5
CPP_BIN_ADD_PARAM='-p'

#one day in seconds
CHECKPOINT_TIME=43200
loop_rerun= true

while $loop_rerun; do

#Running the Simulation with the c++ binaries
./kspaceFirstOrder3D-OMP -i $INPUTFILE -o $OUTPUTFILE $CPP_BIN_ADD_PARAM --checkpoint_interval $CHECKPOINT_TIME --checkpoint_file $CHECKPOINT_FILE

##Running the Simulation with the c++ binaries (without periodically backing up the output)
#./kspaceFirstOrder3D-OMP -i $INPUTFILE -o $OUTPUTFILE $CPP_BIN_ADD_PARAM


TINDEX=$(h5dump -d 't_index' $OUTPUTFILE | grep '(0,0,0)' | cut -c13-)
NT=$(h5dump -d 'Nt' $OUTPUTFILE | grep '(0,0,0)' | cut -c13-)

if [ "$NT" -le "$TINDEX" ]; then
loop_rerun=false
echo "End of Checkpoint simulation, finishing"
fi

echo "$TINDEX $NT $loop_rerun"

done
