#!/bin/bash

WD=$(pwd)
FOLDER="CERENA"
URL="https://github.com/CERENADevelopers/CERENA.git"

if [ ! -d "$FOLDER" ] ; then
    git clone $URL $FOLDER
else
    cd "$FOLDER"
    git pull $URL
fi

cd $WD
matlab -nodisplay -nosplash -nodesktop -r "run('CERENA/CERENA/install_cerena.m');cd 'running_example'; run('modelDef_running.m'); run('simulations_running.m');exit;"
