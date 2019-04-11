def remove_adjacent(list):
    new_list=[]
    new_list.append(list[0])
    temp = 0
    for j in range(1, len(list)):
        if list[j] != new_list[temp]:
            new_list.append(list[j])
            temp += 1
    return new_list

import numpy
import re
for line in open('temp_paths'):
    line=line.strip().split()
    line = list(map(int, line))
    new_line = remove_adjacent(line)
    print(re.sub("[^\d\.]", "", str(new_line)))

