import json
from Bio import SeqIO
with open('data301_400.json') as data_file:
    memo = json.load(data_file); data_file.close()

def p(num):
    row = []
    logic = ['id','blast', 'pfam', 'prosite', 'kegg', 'go', 'comments']
    for ref in logic:
        if memo[str(num)][ref] == 'noHomology':
            row.append('noHomolog')
        else:
            row.append(str(memo[str(num)][ref]))
    return ','.join(row[:])

output=open('Annotations301_400.csv', 'w')
header = ','.join(['id','blast', 'pfam', 'prosite', 'kegg', 'go', 'comments']) + '\n'
output.write(header)
output.write('\n'.join([p(i) for i in range(301, 401)]))
output.close()

print('csv done')
