
#################
# 本脚本未完成！#
#################

# 1. 统计各样品SV数量。
# 2. 统计groups数量。
# 3. 统计各菌株特有的groups数量。
# 4. 统计各菌株特有的变异数量。（人工检查）
# 5. 描述菌株特有变异。

def is_specific(group):
    '''判断该group是不是菌株特有的。
    '''
    return true

samples = []
from collections import Counter
counter = Counter() # SV数量
group_id = 0 # groups 数量
this_group = []
for line in open(""):
    sline = line.rstrip('\n').split('\t')
    if not is_blank():
        sample = re.search(r'?<=sample=[^;]+', sline[5]).group()
        counter[sample] += 1
    if is_blank() and not last_is_blank():
        group_id += 1
