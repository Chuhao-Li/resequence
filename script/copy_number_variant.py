#!/usr/bin/env python

import sys
import numpy as np

samples = open('samples.list').read().rstrip().splitlines()

for sample in samples:
    depth_file = open("cnv/{s}/{s}.depth".format(s=sample))
    low_out = "cnv/{s}/{s}.del.xls".format(s=sample)
    high_out = "cnv/{s}/{s}.high.xls".format(s=sample)
    
    raw_data = []
    lower_limit = 20
    upper_limit_ratio = 1.8 # high depth的定义：中位数的1.7倍。平均值容易受到极端值影响，如果遇到高表达的质粒、大片段的缺失，就会影响鉴定结果。
    min_len = 20 # 只有深度异常的区域大于该值时才报告。
    upper_collector = []
    lower_collector = []
    genome_size = 0
    total_depth = 0
    
    # get low coverage region
    for line in depth_file:
        genome_size += 1
        chro, pos, depth = line.strip().split('\t')
        depth = int(depth)
        pos = int(pos)
        raw_data.append([chro, pos, depth])
        total_depth += depth
    
    average_depth = total_depth/genome_size
    depth_std = np.std([i[2] for i in raw_data])
    depth_median = np.median([i[2] for i in raw_data])
    upper_limit = depth_median*upper_limit_ratio
    
    for chro, pos, depth in raw_data:
        if depth < lower_limit:
            lower_collector.append((chro, pos, depth))
        if depth >upper_limit:
            upper_collector.append((chro, pos, depth))
    
    # join them together
    
    def join(raw_collector):
        collector = [] # [[region1], [region2], ...]
        current = [] # [chro, start, end, depth_sum]
        for c,s,d in raw_collector:
            if current == []:
                current = [c, s, s, d]
            else:
                if c == current[0] and s - current[2] <= 1:
                    current[2] = s
                    current[3] += d
                else:
                    collector.append(current)
                    current = []
        else:
            if current:
                collector.append(current)
        for i in collector:
            l = i[2] - i[1] + 1
            i[3] = round(i[3]/l, 1)
            i.append(l)
        return collector
    
    def print_result(c, f=None, min_len=min_len):
        print('#chro:start-end\taverage_depth\tlength', file=f)
        for c,s,e,a,l in c:
            if l >= min_len:
                print('{}:{}-{}\t{}\t{}'.format(c, s, e, a, l), file=f)
    
    from pprint import pprint
    print('genome size:', genome_size)
    print('average depth:', round(average_depth, 1))
    print('median depth:', round(depth_median, 1))
    print('standard error of depth: ', round(depth_std, 1))
    print('low coverage region(depth lower than {}): '.format(lower_limit))
    with open(low_out, 'w') as f:
        print_result(join(lower_collector), f=f)
    print()
    print('high coverage region(depth over {}): '.format(upper_limit))
    with open(high_out, 'w') as f:
        print_result(join(upper_collector), f=f)

