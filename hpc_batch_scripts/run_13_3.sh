#!/bin/bash
#BSUB -n 1
#BSUB -W 24:00
#BSUB -J PyTreeMotif
#BSUB -o outputs/stdout.%J
#BSUB -e outputs/stderr.%J

source ~/.bashrc
conda activate /usr/local/usrapps/csc522s24/lrwilli7/env_synth_2
python ../13_3_test.py
conda deactivate