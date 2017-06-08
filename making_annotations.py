import json
from Bio.Blast import NCBIXML
from io import StringIO
from Bio import SeqIO

s, f = 301,401
with open('data0_49.json') as data_file:
    data1 = json.load(data_file); data_file.close()

with open('data50_99.json') as data_file:
    data2 = json.load(data_file); data_file.close()

rows = open('Comments301_400.csv','r').read().split('\n')
comments = [v.split(',')[0]+' -> '+v.split(',')[1]+' -> '+','.join(v.split(',')[2:])[:-1] for v in rows[2:]]

queries = list(SeqIO.parse("project_seqs.fasta", format="fasta"))
ints_strings = list(map(str, list(range(100))))
memo = {}

def blast_parser(blast):
    fields = []
    result_handle = StringIO(blast)
    blast_records = NCBIXML.parse(result_handle)
    for rec in blast_records:
        for alignment in rec.alignments:
            for hsp in alignment.hsps:
                #fields = [rec.query_id, str(rec.query_length), alignment.title, alignment.hit_id, alignment.hit_def, alignment.accession, str(hsp.expect)]
                if hsp.expect < 1e-20:
                    fields.append([alignment.hit_id, alignment.hit_def, alignment.accession, str(hsp.expect)])
            if len(fields) >= 5: break
    result_handle.close()
    return fields

memo = {}
for index, seq_num  in zip(ints_strings, list(range(301,351))):
    memo[seq_num] = {'id' : queries[int(index)].id.split(':')[-1],
                     'blast' : blast_parser(data1[index]['blast']),
                     'pfam' : data1[index]['pfam'],
                     'prosite' : data1[index]['prosite'],
                     'kegg' : data1[index]['kegg'],
                     'go' : data1[index]['go'],
                     'comments' : comments[int(index)]}

for index, seq_num  in zip(ints_strings, list(range(351,401))):
    memo[seq_num] = {'id' : queries[int(index)+50].id.split(':')[-1],
                     'blast' : blast_parser(data2[index]['blast']),
                     'pfam' : data2[index]['pfam'],
                     'prosite' : data2[index]['prosite'],
                     'kegg' : data2[index]['kegg'],
                     'go' : data2[index]['go'],
                     'comments' : comments[int(index)+50]}


from collections import defaultdict
pfam = open('pfam2go.txt', 'r').read().split('\n')
pfam2go = defaultdict(set)
for line in pfam[6:-2]:
    pfam2go[line.split(':')[1].split()[0]].add(str(line.split('>')[1][1:]))

prosite = open('prosite2go.txt', 'r').read().split('\n')
prosite2go = defaultdict(set)
for line in prosite[6:-2]:
    prosite2go[line.split(':')[1].split()[0]].add(str(line.split('>')[1][1:]))

for i in range(s, f):
    try:
        curr = list(memo[i]['pfam'].keys())[0]
        if list(pfam2go.get(curr)):
            memo[i]['go'] = list(pfam2go.get(curr))
        else:
            memo[i]['go'] = 'noHomolog'
    except:
        memo[i]['go'] = 'noHomolog'

#just in case prosite has something pfam doesnt. So far no
for i in range(s, f):
    try:
        curr = list(memo[i]['prosite'].keys())[0]
        length = len(memo[i]['go'])
        if list(prosite2go.get(curr)):
            memo[i]['go'] = list(set(memo[i]['go']+list(prosite2go.get(curr))))
            #if length < len(memo[i]['go']):
            #    print(i)
    except:
        pass

#print(len(memo.keys()))
with open('parsed_data301_400.json', 'w') as outfile:
    json.dump(memo, outfile)

print('parsing annotations done')
