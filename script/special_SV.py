#!/usr/bin/env python3
import re
import os
from glob import glob

def from_vcf(vcfs):
    all = [] # {sample: {(chro, start, ref, alt)}}
    for f in vcfs:
        sample = f.split('/')[-1].split('.')[0]
        svs = open(f)
        for sv in svs:
            if sv.startswith('#'):
                continue
            sv = sv.split('\t')
            PE = re.search(r'(?<=PE=)\d+', sv[7]).group()
            end = re.search(r'(?<=END=)\d+', sv[7]).group()
            sv_type = re.search(r'(?<=SVTYPE=)[^;\n]+', sv[7]).group()
            end = int(end)
            start = int(sv[1])
            length = end-start
            sv_tuple = (sv[0], start, end, length, sv_type, 'sample={s};source=delly;PE={pe}'.format(s=sample, pe=PE))
            all.append(sv_tuple)
    return all

def from_cnv(cnvs):
    all = []
    for f in cnvs:
        sample = f.split('/')[-1].split('.')[0] # 这里用‘.’来分割，那么，样品名不能包含.
        if f.endswith('.del.xls'):
            sv_type = "LOW"
        elif f.endswith('.high.xls'):
            sv_type = "HIGH"
        f = open(f)
        next(f)
        for line in f:
            if not line:
                break
            pos, depth, length = re.split(r'\s+', line.rstrip())
            length = int(length)
            chro = pos.split(':')[0]
            start, end = pos.split(':')[1].split('-')
            start = int(start)
            end = int(end)
            all.append((chro, start, end, length, sv_type, 'sample={s};source=samtools;depth={d}'.format(s=sample, d=depth)))
    return all

def cluster_by_pos(all_svs):
    # all_svs = [] # [(chro, start, end, length, sv_type, others), ...]
    last = 0
    last_chro = ''
    if not os.path.exists('report/SV'):
        os.makedirs('report/SV')
    out = open('report/SV/clustered.SV.txt', 'w')
    print("chromosome\tstart\tend\tlength\tsv_type\tothers", file=out)
    for i in sorted(all_svs, key=lambda x: tuple(x[:4])): # 函数的关键是这个排序。
        this_chro = i[0]
        this = i[1]
        if this_chro != last_chro:
            print('-\n-{}-\n-'.format(this_chro), file=out)
        elif abs(this - last) >= 100:
            print('-', file=out)
        print('\t'.join([str(j) for j in i]), file=out)
        last_chro = this_chro
        last = this
    out.close()

vcfs = glob("delly/*/*.filtered.vcf")
cnvs = glob("cnv/*/*.xls")
all_svs = from_vcf(vcfs) + from_cnv(cnvs)
cluster_by_pos(all_svs)
