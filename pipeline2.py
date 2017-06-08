from prody import *
from matplotlib.pylab import *
from Bio.Blat import NCBIWWW
from Bio import SeqIO
from Bio.ExPASy import ScanProsite
from bioservices.kegg import KEGG
import json
k = KEGG()

"""
*How annotations are queried*

blast -> seq
pfam -> seq
prosite -> seq
kegg -> original gene annotation
go -> pfam and prosite ontology later
"""

DNA=open('UP000006737.fasta.txt', 'r').read().split('>')
queries = list(SeqIO.parse("project_seqs.fasta", format="fasta"))
dna = DNA[301:401] #first index is blank from split so index is "normal" for sequences
dna = dna[50:] #can be used to pick which sequences to query
queries = queries[50:] #can be used to pick which sequences to query
fasta = ['>'+dna[i] for i in range(len(dna))]
dna = [''.join(v.split('\n')[1:-1]) for v in dna]
gn = [''.join(v.split('\n')[0]).split('GN=')[1].split()[0] for v in fasta]

def blast(query):
    try:
        result_handle = NCBIWWW.qblast("blastp", "nr", query.seq)
        result_handle = result_handle.read()
        return result_handle
    except:
        return 'noHomology'

def pfam(seq):
    try:
        matches = searchPfam(seq)
        return matches
    except:
        return 'noHomology'

def prosite(seq):
    try:
        handle = ScanProsite.scan(seq=seq, lowscore=10) #add lowscore=# to get the next # 'possible' additional hits
        result = ScanProsite.read(handle)
        if result == []:
            return 'noHomology'
        else:
            return result
    except:
        return 'noHomology'

def kegg(gn):
    try:
        hit = k.get_pathway_by_gene(gn, "acb")
        if not hit:
            return 'noHomology'
        else:
            return hit
    except:
        return 'noHomology'

def go(go_id):
    "Will be inputed later"
    pass

memo = {num:{'blast':blast(queries[num]), 'pfam':pfam(seq), 'prosite':prosite(seq), 'kegg':kegg(gn[num]), 'go':'future go data'} for num, seq in list(enumerate(dna[:]))}

with open('data50_99.json', 'w') as outfile:
    json.dump(memo, outfile)

print('second half done')
