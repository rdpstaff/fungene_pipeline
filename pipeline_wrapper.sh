#!/bin/bash

export PATH=$PATH:/bin:/home/web/bin
. /home/web/gridware_scripts/python-virtenv/bin/activate
#. /home/web/gridware_scripts/fungene_pipeline_staging/activate
source /home/sgeadmin/ge6.2u3/rack_cell/common/settings.sh

echo $*
/home/web/gridware_scripts/fungene_pipeline/fgp_wrapper.py $*
#testing /home/web/gridware_scripts/fungene_pipeline_staging/fgp_wrapper.py $*
