#!/usr/bin/python3
import sys
import re
reads_threshold = 100
FILTER_threshold = 'PASS'
print("#threshold: SUPPORT_READS>{} and FILTER=={}".format(str(reads_threshold), FILTER_threshold))

# print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS-7PDE_FDSW202397866-1r.markDuplicates")
# print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tLENTH\tSUPPORT_READS")
# for line in sys.stdin:
#     if line.startswith('#'):
#         continue
#     line = line.rstrip().split('\t')
#     END = re.search(r'END=(\d+?);', line[7]).group(1)
#     PE = re.search(r'PE=(\d+?);', line[7]).group(1)
#     line[7] = str(int(END) - int(line[1])) #length
#     line[8] = PE
#     if int(PE) > reads_threshold and line[6] == FILTER_threshold:
#         print('\t'.join(line[:9]))

print("#CHROM:START-END\tLENGTH\tNAME\tSUPPORTED_READS")
for line in open(sys.argv[1]):
    if line.startswith('#'):
        continue
    line = line.rstrip().split('\t')
    CHROM = line[0]
    START = line[1]
    END = re.search(r'END=(\d+?);', line[7]).group(1)
    NAME = line[2]
    PE = re.search(r'PE=(\d+?);', line[7]).group(1)
    if int(PE) > reads_threshold and line[6] == FILTER_threshold:
        print('{}:{}-{}\t{}\t{}\t{}'.format(CHROM, START, END, int(END)-int(START), NAME, PE))
