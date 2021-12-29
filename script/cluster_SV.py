#!/usr/bin/env python3
import sys
import re

all = []
for f in sys.argv[1:]:
    sample = f.split('/')[-1].split('.')[0]
    f = open(f)
    next(f)
    for line in f:
        if not line:
            break
        pos, depth, length = re.split(r'\s+', line.rstrip())
        chro = pos.split(':')[0]
        start, end = pos.split(':')[1].split('-')
        start = int(start)
        end = int(end)
        all.append((chro, start, end, depth, length, sample))

all.sort()
last = 0
chro = ''
for i in all:
    if chro != i[0]:
        chro = i[0]
        last = 0
        print()
    if i[1] - last > 100:
        print()
    print('\t'.join([str(j) for j in i]))
    last = i[1]
