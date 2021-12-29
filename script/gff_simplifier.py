import re
infile = 'PE100_pass.bed.gff'
START = None
for line in open(infile):
    if line.startswith('#'):
        continue
    line = line.rstrip().split('\t')
    # print(line)
    ID = re.search(r'ID=([^;]*)', line[13]).group(1)
    product = re.search(r'product=([^;]*)', line[13]).group(1)
    new = line[:5]
    new.append(ID)
    new.append(product)
    if START != line[1]:
        print()
        START = line[1]
    print('\t'.join(new))
