# compare 2 pipeline
# compare base on clusters_2


import csv
import sys
file1 = sys.argv[1]
file2 = sys.argv[2]

clusters_1 = []
for row in csv.reader(open(file1, 'r')):
    this_cluster = set()
    if row[0] == 'Gene':
        continue
    for cell in row[8:]:
        genes = cell.split('\t')
        for gene in genes:
            if gene == '':
                continue
            this_cluster.add(gene)
    clusters_1.append([row[3], this_cluster])  # row[3] is the number of sequences


clusters_2 = []
for row in csv.reader(open(file2, 'r')):
    this_cluster = set()
    if row[0] == 'Gene':
        continue
    for cell in row[8:]:
        genes = cell.split('\t')
        for gene in genes:
            if gene == '':
                continue
            this_cluster.add(gene)
    clusters_2.append([row[3], this_cluster])



same = 0
diff = 0
total = 0
for i in clusters_2:
    maxx = 0
    match = None
    
    for j in clusters_1:
        intersection = i[1].intersection(j[1])
        if len(intersection) > maxx:
            maxx = len(intersection)
            match = j

    if match == None:
        diff += 1
        continue    
    
    diff_1 = i[1].difference(match[1])
    diff_2 = match[1].difference(i[1])

    ls = [i[0]]

    if len(diff_1) == 0 and len(diff_2) == 0:
        same += 1
        ls.append('Same')
    else:
        diff += 1
        ls.append('Diff')

    total +=1
    print('\t'.join(ls))

print(str(same) + '\t' + str(diff) + '\t'+ str(total))