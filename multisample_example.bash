#!/bin/bash
cd /srv/qgen/qiaseq-dna
source /srv/conda/bin/activate base
time python run_qiaseq_dna.py smcounter2.multisample_example.params.txt v2 single out2v6_{0} NEB_S2 NEB_S2_bis &> run_v6.2.2.log &