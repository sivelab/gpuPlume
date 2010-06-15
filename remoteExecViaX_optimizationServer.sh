#!/bin/bash
rm -f /scratch/popsampling_log.txt
./populationSampleOpt/plumeOptimization --usedb >& /scratch/popsampling_log.txt
