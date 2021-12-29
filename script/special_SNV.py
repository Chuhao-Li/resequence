#!/usr/bin/env python3
# description: compare paired samples to get special snvs. 

def collect_snvs(samples):
    snvs_dict = {} # {'sample':('sample', chro, start, end, alt), ...}
    for sample in samples:
        snvs = open('BWA/{s}/{s}.final.vcf'.format(s=sample))
        snvs_dict[sample] = set()
        for snv in snvs:
            if snv.startswith('#'):
                continue
            snv = snv.split('\t')
            snv_tuple = (sample, snv[0], snv[1], snv[3], snv[4])
            snvs_dict[sample].add(snv_tuple)
    return snvs_dict

def create_gff_db(gfffile):
    import os
    dbfn=gfffile + '.db'
    if os.path.exists(dbfn):
        db = gffutils.FeatureDB(dbfn, keep_order=True)
    else:
        db = gffutils.create_db(gfffile, dbfn=dbfn, force=True, keep_order=True,
                merge_strategy='merge', sort_attribute_values=True)
    return db

def cluster_by_pos(snvs_dict, samples, gff):
    import re
    db = create_gff_db(gff)
    allsets = snvs_dict.values()
    all_snvs = [] # [('sample', chro, start, ref, alt), ...]
    anno = get_anno(samples)
    for i in allsets:
        for j in i:
            all_snvs.append(j)
    
    if not os.path.exists('report/SNV'):
        os.makedirs('report/SNV')
    out = open('report/SNV/clustered.SNV.txt', 'w')
    # print("# SNV(Single Nucleotide Variant) 结果文件。", 
    #         "# 每个样品都与对应的参考序列进行了比对，使用bcftools获得了SNVs。",
    #         "# 按照在参考基因组序列上的位置对SNVs进行了排列，间隔大于100bp的用空行隔开形成多个小组。同一个SNVs应该落在同一个小组中。", 
    #         "# 如果该小组中只包含了非野生型的SNVs，更有可能与表型相关。", 
    #         "# sample\tchromosome\tposition\tREF\tALT\tannotation", sep = '\n')
    print("sample\tchromosome\tposition\tREF\tALT\tannotation", file=out)
    last = 0
    last_chro = ''
    for i in sorted(all_snvs, key=lambda x: (x[1], int(x[2]))):
        this_chro = i[1]
        this = int(i[2])
        if this_chro != last_chro:
            print('-\n-{}\n-'.format(this_chro), file=out)
        elif abs(this - last) >= 100:
            print('-', file=out)
        if i in anno:
            gene = re.match(r'[^:#]+', anno[i][2]).group()
            try:
                product = next(db.children(gene, featuretype='CDS')).attributes['product'][0]
            except(StopIteration):
                product = 'Unknow'
            a = "SNVType={};ChangeGene={};SNVScore={};readDepth={};product={}".format(anno[i][1], 
                    gene, 
                    anno[i][-2],
                    anno[i][-1], 
                    product)
        else:
            a = 'SNVtype=intergenic'
        print('\t'.join(i) + '\t' + a, file=out)
        last_chro = this_chro
        last = this

def get_special(snvs_dict):
    allsets = snvs_dict.values()
    counts = {}
    for i in allsets:
        for j in i:
            if j in counts:
                counts[j] += 1
            else:
                counts[j] = 1
    # for i in sorted(counts, key = lambda x: counts[x]):
    #     print(i, counts[i])
    special = {i for i in counts if counts[i]==1}
    return special

def get_anno(samples):
    m = {} # {("chro", "pos", "ref", "alt"):"annotation", ...}
    for sample in samples:
        for line in open('snv_anno/{s}/{s}_snp.exonic_variant_function'.format(s=sample)):
            line = line.rstrip().split('\t')
            m[(sample, line[3], line[4], line[6], line[7])] = line
    return m

def pair_special(snvs_dict, special):
    for line in mapping:
        line = line.rstrip().split('\t')
        print(line[0])
        s1,s2 = line[1:]
        s1_anno = get_anno(s1)
        s2_anno = get_anno(s2)
        print("{}: ".format(s1))
        for i in snvs_dict[s1] - snvs_dict[s2]:
            if i in special:
                if i in s1_anno:
                    print(s1_anno[i])
                else:
                    print('{}: intergenic'.format(i))
        print("{}: ".format(s2))
        for i in snvs_dict[s2] - snvs_dict[s1]:
            if i in special:
                if i in s2_anno:
                    print(s2_anno[i])
                else:
                    print('{}: intergenic'.format(i))
        print()

if __name__ == '__main__':
    import gffutils
    import os
    gff = "reference/reference.gff"
    samples = os.listdir('BWA')
    # samples = ['plsX', 'ycgR']
    print(samples)
    snvs_dict = collect_snvs(samples)
    cluster_by_pos(snvs_dict, samples, gff)
    # special = get_special(snvs_dict)
    # pair_special(snvs_dict, special) # 这种方法会忽略掉不同样品共有的突变。

