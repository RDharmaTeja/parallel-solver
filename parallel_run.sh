#!/bin/bash

#output_dir='results'
ex='parallel_solver.out'
runlog='out'
#vislog='vislog'
#visdaemon='visdaemon.py'
echo | tee -a $runlog
echo "Run: `date +%Y/%m/%d-%H/%M/%S`" | tee -a $runlog


#if ls output*.*vtk 1> /dev/null 2>&1; then
    ## Output files exist
    ## Back them up
    #echo 'Backing up existing data' | tee -a $runlog
    #. backupcurrent.sh
    #echo "Previous output files moved to $backup_dir""." | tee -a $runlog
#fi

if [ -f $ex ]; then
    #echo 'Clearing old log files' | tee -a $runlog
    #. clearlogs.sh
    #echo 'Running solver.' | tee -a $runlog
    #python $visdaemon start './' ./$vislog  
    #echo >> $runlog
    #echo 'Solver output:' >> $runlog
    #echo >> $runlog
    #espeak -a10 'Running solver.'
    mpirun -np 2 ./$ex >> $runlog
fi

echo | tee -a $runlog
echo 'End of run script.' | tee -a $runlog
