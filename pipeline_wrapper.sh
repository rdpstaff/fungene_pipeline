#!/bin/bash

. /home/web/gridware_scripts/python-virtenv/bin/activate
source /home/sgeadmin/ge6.2u3/rack_cell/common/settings.sh
export PATH=$PATH:/home/web/bin

echo $*
/scratch/fishjord/fungene_pipeline_scripts/fgp_wrapper.py $*
