#!/bin/bash

echo "************************************************************"
echo ""
echo "          SUBMITTING cfNIPT PIPELINE TO QUEUE"
echo ""
echo "************************************************************"
sbatch ./cfnipt
echo "************************************************************"
sleep 3

echo "     Contact machan@rm.dk if you experience problems "
echo "------------------------------------------------------------"

sleep 3
echo ""
echo " The analyses will begin shortly if resources are available "
sleep 5
echo " Jobs running or scheduled on GenomeDK for current user...  "

echo "------------------------------------------------------------"
squeue --me
echo "------------------------------------------------------------"
sleep 5
echo "************************************************************"
echo " You may close window (closing automatically in 30 seconds) "
echo "************************************************************"

sleep 30
