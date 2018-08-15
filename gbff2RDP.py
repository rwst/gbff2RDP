from sys import *
from subprocess import run
from Bio import Seq, SeqRecord, SeqIO, Entrez
from Bio.Emboss.Applications import SeqretCommandline
from anytree import Node, find, findall, PreOrderIter, RenderTree, AsciiStyle
from urllib.error import HTTPError
import xml.etree.ElementTree as ET

def disp_status(timelapse, timeout):
    if timelapse and timeout:
        percent = 60 * (float(timelapse)/float(timeout))
    stdout.write("Progress : ["+"*"*int(percent)+" "*(60-int(percent-1))+"]")
    stdout.flush()
    stdout.write("\r\r")

class TaxNode(Node):
    separator = ';'

def insert_data(rec, tax, taxtree):
    """
    Insert the taxonomic data in tax associated with the specific
    fasta record rec into rec's description and into our own taxtree.

    We assume there is no taxon with two parents.
    """
    terms = ''
    node = taxtree
    for line in tax:
        terms = terms + ' ' + line.lstrip(' ').rstrip()
    terms = list(filter(lambda x: len(x), terms.split(';')))
    for t in terms:
        term = t.lstrip().rstrip('.').replace(' ', '_')
        f = find(node, lambda n: n.name == term)
        if f:
            node = f
        else:
            node = TaxNode(term, parent=node)
        if rankmap.get(term) == 'genus':
            break
    if rankmap.get(term) != 'genus':
        term = term + '_genus'
        f = find(node, lambda n: n.name == term)
        if f:
            node = f
        else:
            node = TaxNode(term, parent=node)
        if not rankmap.get(term):
            rankmap[term] = 'genus'
    rec.description = str(node).split("'")[1][1:]

arg = argv[1]
Entrez.email = "ralf@ark.in-berlin.de"
Entrez.api_key = '788bafaf966e618652a65887a6314f458a08'
run(['seqret', arg, 'tmp.fa'])
rdict = SeqIO.to_dict(SeqIO.parse('tmp.fa', 'fasta'))
hdl = open(arg, 'r')
s = hdl.readlines()

seqid = None
org_seen = False
rankmap = {}
ranks = set(['domain', 'phylum', 'class', 'subclass', 'order',
             'suborder', 'family', 'genus'])
try:
    hdl = open('taxa_rank_db.txt', 'r')
except FileNotFoundError:
    pass
else:
    for line in hdl.readlines():
        p = line.rfind(' ')
        name = line[:p]
        rank = line[p+1:].rstrip()
        rankmap[name.replace(' ', '_')] = rank
        hdl.close()

taxtree = TaxNode("Root")
tax = []
count = 0
rcount = len(rdict)
for line in s:
    line = line.rstrip()
    if org_seen:
        if line.startswith('REFERENCE'):
            org_seen = False
            rec = rdict.get(seqid)
            if rec:
                insert_data(rec, tax, taxtree)
                count = count + 1
                disp_status(count, rcount)
            tax = []
            continue
        tax.append(line)
    if line.startswith('ACCESSION'):
        seqid = line.split()[1].lstrip().rstrip()
        continue
    if line.startswith('  ORGANISM'):
        org_seen = True
        continue

rlist = list(rdict.values())
print('\nNumber of sequences: {}, taxonomic annotations added: {}'.
        format(len(rlist), count))
SeqIO.write(rlist, 'trainset.fa', 'fasta')

def fill_rankmap():
    children = set(node.name for node in PreOrderIter(taxtree))
    new_taxa = children.difference(rankmap.keys())
    if new_taxa is None or len(new_taxa) == 0:
        return
    print(new_taxa)
    return
    print('Querying Entrez for {} new taxa...'.format(len(new_taxa)))
    ids = []
    for t in new_taxa:
        hdl = Entrez.esearch(db="taxonomy", term=t,
                field='Scientific Name', usehistory="y")
        entry = Entrez.read(hdl)
        l = entry.get('IdList')
        if len(l) == 0:
            print('{} not found'.format(t))
        else:
            ids = ids + [l[0]]
        if len(l) > 1:
            print('found {} ids for {}'.format(len(l), t))
        disp_status(len(ids), len(new_taxa))
    search_results = Entrez.read(Entrez.epost("taxonomy", id=",".join(ids)))
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    count = len(ids)
    batch_size = 100
    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)
        print("Going to download record %i to %i" % (start+1, end))
        attempt = 0
        while attempt < 3:
            attempt += 1
            try:
                fetch_handle = Entrez.esummary(db="taxonomy",
                                   retstart=start, retmax=batch_size,
                                   webenv=webenv, query_key=query_key)
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt)
                    time.sleep(15)
                else:
                    raise
        data = ET.fromstring(fetch_handle.read())
        fetch_handle.close()
        for entry in data.findall('./DocSum'):
            for item in entry.findall('./Item'):
                if item.get('Name') == 'ScientificName':
                    name = item.text
                if item.get('Name') == 'Rank':
                    rank = item.text
            if rank is None or len(rank) == 0:
                rank = 'unknown'
            rankmap[name] = rank
    hdl = open('taxa_rank_db.txt', 'w')
    for k,v in rankmap.items():
        hdl.write('{} {}\n'.format(k, v))

def prune_taxtree():
    for node in PreOrderIter(taxtree,
            filter_=lambda n: rankmap.get(n.name) == 'genus'
            and len(n.children) > 0):
        node.children = []
#    for node in PreOrderIter(taxtree,
#            filter_=lambda n: len(n.children) == 0
#            and rankmap.get(n.name) != 'genus'):

fill_rankmap()
prune_taxtree()
t = RenderTree(taxtree, style=AsciiStyle()).by_attr()
hdl = open('taxtree.txt', 'w')
hdl.write(t+'\n')

print('Writing taxid db')
count = -1
m = {}
tfile = open('trainset_db_taxid.txt', 'w')
rankmap['Root'] = 'rootrank'
for node in PreOrderIter(taxtree):
    count = count + 1
    m[node.name] = count
    if count == 0:
        parentno = -1
    else:
        parentno = m.get(node.parent.name)
    tfile.write('{}*{}*{}*{}*{}\n'.format(count,
        node.name, parentno, node.depth, rankmap.get(node.name)))
