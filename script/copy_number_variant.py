#!/usr/bin/env python
usage = 'copy_number_variant.py [*.depth] [low_coverage_out_file] [high_coverage_out_file]'
import sys
import numpy as np
if len(sys.argv) == 1:
    print(usage)
    exit()

depth_file = open(sys.argv[1])
raw_data = []
lower_limit = 20
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
upper_limit = average_depth*2

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

def print_result(c, f=None):
    print('#chro:start-end\taverage_depth\tlength', file=f)
    for c,s,e,a,l in c:
        print('{}:{}-{}\t{}\t{}'.format(c, s, e, a, l), file=f)

from pprint import pprint
print('genome size:', genome_size)
print('average depth:', round(average_depth, 1))
print('standard error of depth: ', round(depth_std, 1))
print('low coverage region(depth lower than {}): '.format(lower_limit))
print_result(join(lower_collector))
with open(sys.argv[2], 'w') as f:
    print_result(join(lower_collector), f=f)
print()
print('high coverage region(depth over {}): '.format(upper_limit))
print_result(join(upper_collector))
with open(sys.argv[3], 'w') as f:
    print_result(join(upper_collector), f=f)

