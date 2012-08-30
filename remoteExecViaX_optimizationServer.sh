#!/bin/bash
rm -f /scratch/popsampling_log.txt
time ./populationSampleOpt/plumeOptimization --optfile=/scratch/optTest8D_stage3.opt >& /scratch/popsampling_log.txt
