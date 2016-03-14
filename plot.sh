. clearlogs.sh
cp buffer/process_0$1/output*.fvtk . 
python $visdaemon start './' ./$vislog
