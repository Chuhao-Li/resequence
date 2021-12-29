#!/usr/bin/env python3
import sys
import pandas as pd

qc_dir = sys.argv[1].rstrip('/') + '/'
qc_file = qc_dir + 'multiqc_data/multiqc_fastqc.txt'
qc_file = pd.read_csv(qc_file, sep='\t')

ok = True
for i, row in qc_file.iterrows():
    if row.values.count('fail') >= 1:
        ok = False
        break

if ok == True:
    with open(qc_dir + 'qc.ok', 'w') as f:
        f.write('')

