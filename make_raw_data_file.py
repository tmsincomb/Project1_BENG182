import json
from Bio import SeqIO

with open('data0_49.json') as data_file:
    data1 = json.load(data_file); data_file.close()

with open('data50_99.json') as data_file:
    data2 = json.load(data_file); data_file.close()

rows = open('Comments301_400.csv','r').read().split('\n')
comments = [v.split(',')[0]+' -> '+v.split(',')[1]+' -> '+','.join(v.split(',')[2:])[:-1] for v in rows[2:]]
queries = list(SeqIO.parse("project_seqs.fasta", format="fasta"))
ints_strings = list(map(str, list(range(100))))
memo = {}

for index, seq_num  in zip(ints_strings, list(range(301,351))):
    memo[seq_num] = {'id' : queries[int(index)].id.split(':')[-1],
                     'blast' : data1[index]['blast'],
                     'pfam' : data1[index]['pfam'],
                     'prosite' : data1[index]['prosite'],
                     'kegg' : data1[index]['kegg'],
                     'go' : data1[index]['go'],
                     'comments' : comments[int(index)]}

for index, seq_num  in zip(ints_strings, list(range(351,401))):
    memo[seq_num] = {'id' : queries[int(index)+50].id.split(':')[-1],
                     'blast' : data2[index]['blast'],
                     'pfam' : data2[index]['pfam'],
                     'prosite' : data2[index]['prosite'],
                     'kegg' : data2[index]['kegg'],
                     'go' : data2[index]['go'],
                     'comments' : comments[int(index)+50]}

output = open('raw_data301_400.txt', 'w')
for i in range(301, 401):
    target=memo[i]['id']+'.txt'
    logic = ['id','blast', 'pfam', 'prosite', 'kegg', 'go', 'comments']
    output.write('\n\n'.join([ref+'\n'+str(memo[i][ref]) for ref in logic]))
output.close()

print('raw data done')
