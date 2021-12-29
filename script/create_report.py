#!/usr/bin/env python3

# 使用pandas把SV输出结果转化为html格式。
# 在每行末尾增加web_igv的链接。

import sys
import os
import json
import glob
import pandas as pd
pd.set_option('display.max_colwidth', -1)

genome = 'reference/reference.fasta'
gff = 'reference/reference.gff'

# 总是展示所用样品的比对信息。
bams = glob.glob("BWA/*/*.bam")
# /tool/igv/?genome=public/assembly/Dickeya_zeae_EC1_genomic.fna&gff=public/assembly/Dickeya_zeae_EC1_genomic.gff&bams=noodle/2020-08-13_yufan_resequence/BWA/S-15D/S-15D.bam&range=NZ_CP006929.1:1000-1100
prefix = '/tool/igv/?genome={genome}&gff={gff}&{bam}&range='.format(genome=genome, gff=gff, bam='&'.join(['bams=' + i for i in bams]))

var = sys.argv[1] # SV/SNV
infile = 'report/{}/clustered.{}.txt'.format(var, var)
if not os.path.exists(infile):
    sys.stderr.write("input file not found! \n")
    exit()

a = pd.read_csv(infile, sep='\t')

pos_list = []
for i,j in a.iterrows():
    if var == "SV":
        if j['chromosome'].startswith('-'):
            pos_list.append('')
            continue
        pos = j['chromosome'] + ':' + str(int(j['start'])) + '-' + str(int(j['end']))
    elif var == "SNV":
        if j['sample'].startswith('-'):
            pos_list.append('')
            continue
        pos = str(j['chromosome']) + ':' + str(int(j['position']))
    pos_list.append('<a target="_blank" rel="noreferrer" href="{}">view in igv</a>'.format(prefix + pos))

a['igv_link'] = pd.Series(pos_list)
a = a.fillna('')
a.to_html('report/{}/clustered.{}.html'.format(var, var), float_format=lambda x: '{:.0f}'.format(x), escape=False)
