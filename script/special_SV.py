#!/usr/bin/env python3
import sys
import re
import os

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
            end = int(end)
            start = int(sv[1])
            length = end-start
            sv_tuple = (sv[0], start, end, length, 'sample={s};source=delly;PE={pe}'.format(s=sample, pe=PE))
            all.append(sv_tuple)
    return all

def from_cnv(cnvs):
    all = []
    for f in cnvs:
        sample = f.split('/')[-1].split('.')[0]
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
            all.append((chro, start, end, length, 'sample={s};source=samtools;depth={d}'.format(s=sample, d=depth)))
    return all

def cluster_by_pos(all_svs):
    # all_svs = [] # [(chro, start, end, length, others), ...]
    last = 0
    last_chro = ''
    # print("# SV(Structural Variant) 结果文件。",
    #         "# 每个样品都与对应的参考序列进行了比对，使用samtools得到低覆盖度区域(使用source=samtools标注的SV)，使用delly获得了SVs(使用source=delly标注的SV)。",
    #         "# 按照在参考基因组序列上的位置对SVs进行了排列，间隔大于100bp的用空行隔开形成多个小组。同一个SV应该落在同一个小组中。",
    #         "# 如果某个位置上，所有样品都有同样的SV，那可能是野生型本来就有的突变。", 
    #         "chromosome\tstart\tend\tlength\tothers", sep = '\n')
    if not os.path.exists('report/SV'):
        os.makedirs('report/SV')
    out = open('report/SV/clustered.SV.txt', 'w')
    print("chromosome\tstart\tend\tlength\tothers", file=out)
    for i in sorted(all_svs, key=lambda x: tuple(x[:4])):
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

vcfs = []
cnvs = []
for i in sys.argv[1:]:
    if i.endswith('.vcf'):
        vcfs.append(i)
    elif i.endswith('.xls'):
        cnvs.append(i)
all_svs = from_vcf(vcfs) + from_cnv(cnvs)
cluster_by_pos(all_svs)
