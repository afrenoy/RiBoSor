#!/opt/local/bin/python2.7
# -*- coding:UTF-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools
import csv
import os
import re
import argparse
import sys
sys.path.append('/Users/antoine/code/bioinfo')
import tunestopfs

# Parse the arguments
#parser = argparse.ArgumentParser(description='Turn DO(t) and Fluo(t) csv files into one Fluo(DO) csv file')
#parser.add_argument("--do", type=str, help="input do file")
#parser.add_argument("--fluo", type=str, help="input fluo file")
#parser.add_argument("--output", type=str, help="output file")
#args = parser.parse_args()

######################    Global variables   ########################

#****** Fonction codon synonymes ******
allsyn=tunestopfs.SynonymousCodons

#****** Fonction codon similaires ******
allnonsyn=dict()

allnonsyn['TTC']=['TAT','TAC','TGG','TGG','ATC','ATA','ATT','TTG','TTA','CTC','CTG', 'CTA','CTT']
allnonsyn['TTT']=['TAT','TAC','TGG','TGG','ATC','ATA','ATT','TTG','TTA','CTC','CTG', 'CTA','CTT']
allnonsyn['TTG']=['ATC','ATA','ATT','ATG','GTC','GTA','GTG', 'GTT','TTT','TTC']
allnonsyn['TTA']=['ATC','ATA','ATT','ATG','GTC','GTA','GTG', 'GTT','TTT','TTC']
allnonsyn['CTT']=['ATC','ATA','ATT','ATG','GTC','GTA','GTG', 'GTT','TTT','TTC']
allnonsyn['CTC']=['ATC','ATA','ATT','ATG','GTC','GTA','GTG', 'GTT','TTT','TTC']
allnonsyn['CTA']=['ATC','ATA','ATT','ATG','GTC','GTA','GTG', 'GTT','TTT','TTC']
allnonsyn['CTG']=['ATC','ATA','ATT','ATG','GTC','GTA','GTG', 'GTT','TTT','TTC']
allnonsyn['ATT']= ['GTC','GTA','GTG', 'GTT','TTG','TTA','CTC','CTG', 'CTA','CTT','TTT','TTC','TGT','TGC']
allnonsyn['ATC']= ['GTC','GTA','GTG', 'GTT','TTG','TTA','CTC','CTG', 'CTA','CTT','TTT','TTC','TGT','TGC']
allnonsyn['ATA']= ['GTC','GTA','GTG', 'GTT','TTG','TTA','CTC','CTG', 'CTA','CTT','TTT','TTC','TGT','TGC']
allnonsyn['ATG']= ['GTC','GTA','GTG', 'GTT','TTG','TTA','CTC','CTG', 'CTA','CTT','TTT','TTC','TGT','TGC']
allnonsyn['GTT']= ['ATC','ATA','ATT','TTG','TTA','CTC','CTG', 'CTA','CTT','ATG','GCT','GCC','GCA','GCG']
allnonsyn['GTC']= ['ATC','ATA','ATT','TTG','TTA','CTC','CTG', 'CTA','CTT','ATG','GCT','GCC','GCA','GCG']
allnonsyn['GTA']= ['ATC','ATA','ATT','TTG','TTA','CTC','CTG', 'CTA','CTT','ATG','GCT','GCC','GCA','GCG']
allnonsyn['GTG']= ['ATC','ATA','ATT','TTG','TTA','CTC','CTG', 'CTA','CTT','ATG','GCT','GCC','GCA','GCG']
allnonsyn['TCT']= ['ACC','ACT','ACA','ACG', 'GCT','GCC','GCA','GCG','AAC','AAT','GAC','GAT']
allnonsyn['TCC']= ['ACC','ACT','ACA','ACG', 'GCT','GCC','GCA','GCG','AAC','AAT','GAC','GAT']
allnonsyn['TCA']= ['ACC','ACT','ACA','ACG', 'GCT','GCC','GCA','GCG','AAC','AAT','GAC','GAT']
allnonsyn['TCG']= ['ACC','ACT','ACA','ACG', 'GCT','GCC','GCA','GCG','AAC','AAT','GAC','GAT']
allnonsyn['CCT']= ['GCT','GCC','GCA','GCG','GAC','GAT','CAG','CAA','GAA','GAG']
allnonsyn['CCC']= ['GCT','GCC','GCA','GCG','GAC','GAT','CAG','CAA','GAA','GAG']
allnonsyn['CCA']= ['GCT','GCC','GCA','GCG','GAC','GAT','CAG','CAA','GAA','GAG']
allnonsyn['CCG']= ['GCT','GCC','GCA','GCG','GAC','GAT','CAG','CAA','GAA','GAG']
allnonsyn['ACT']= ['TCT','TCA','TCG', 'TCC', 'AGT', 'AGG','GTC','GTA','GTG', 'GTT','GAA','GAG','GGG','GGC','GGA','GGT']
allnonsyn['ACC']= ['TCT','TCA','TCG', 'TCC', 'AGT', 'AGG','GTC','GTA','GTG', 'GTT','GAA','GAG','GGG','GGC','GGA','GGT']
allnonsyn['ACA']= ['TCT','TCA','TCG', 'TCC', 'AGT', 'AGG','GTC','GTA','GTG', 'GTT','GAA','GAG','GGG','GGC','GGA','GGT']
allnonsyn['ACG']= ['TCT','TCA','TCG', 'TCC', 'AGT', 'AGG','GTC','GTA','GTG', 'GTT','GAA','GAG','GGG','GGC','GGA','GGT']
allnonsyn['GCT']= ['GTC','GTA','GTG', 'GTT','GGG','GGC','GGA','GGT','TGT','TGC','ACC','ACT','ACA','ACG']
allnonsyn['GCC']= ['GTC','GTA','GTG', 'GTT','GGG','GGC','GGA','GGT','TGT','TGC','ACC','ACT','ACA','ACG']
allnonsyn['GCA']= ['GTC','GTA','GTG', 'GTT','GGG','GGC','GGA','GGT','TGT','TGC','ACC','ACT','ACA','ACG']
allnonsyn['GCG']= ['GTC','GTA','GTG', 'GTT','GGG','GGC','GGA','GGT','TGT','TGC','ACC','ACT','ACA','ACG']
allnonsyn['TAT']= ['TGG','TGG','ATG','TTG','TTA','CTC','CTG', 'CTA','CTT','ATC','ATA','ATT']
allnonsyn['TAC']= ['TGG','TGG','ATG','TTG','TTA','CTC','CTG', 'CTA','CTT','ATC','ATA','ATT']
allnonsyn['CAT']= ['TAT','TAC','AGG','CGC','CGA','CGT','AGA','CGC','AAC','AAT','GAA','GAG']
allnonsyn['CAC']= ['TAT','TAC','AGG','CGC','CGA','CGT','AGA','CGC','AAC','AAT','GAA','GAG']
allnonsyn['CAA']= ['AGG','CGC','CGA','CGT','AGA','CGC','AAC','AAT','GAC','GAT','CAT','CAC']
allnonsyn['CAG']= ['AGG','CGC','CGA','CGT','AGA','CGC','AAC','AAT','GAC','GAT','CAT','CAC']
allnonsyn['AAT']= ['TCT','TCA','TCG', 'TCC', 'AGT', 'AGG','CAT','CAC','GAA','GAG','GGG','GGC','GGA','GGT']
allnonsyn['AAC']= ['TCT','TCA','TCG', 'TCC', 'AGT', 'AGG','CAT','CAC','GAA','GAG','GGG','GGC','GGA','GGT']
allnonsyn['AAA']= ['AGG','CGC','CGA','CGT','AGA','CGC','GAA','GAG','CAG','CAA','TCT','TCA','TCG', 'TCC', 'AGT', 'AGG']
allnonsyn['AAG']= ['AGG','CGC','CGA','CGT','AGA','CGC','GAA','GAG','CAG','CAA','TCT','TCA','TCG', 'TCC', 'AGT', 'AGG']
allnonsyn['GAT']= ['AAC','AAT','GAA','GAG','TCT','TCA','TCG', 'TCC', 'AGT', 'AGG','CAG','CAA']
allnonsyn['GAC']= ['AAC','AAT','GAA','GAG','TCT','TCA','TCG', 'TCC', 'AGT', 'AGG','CAG','CAA']
allnonsyn['GAA']= ['AAA','AAG','GCT','GCC','GCA','GCG','AAC','AAT','TCT','TCA','TCG', 'TCC', 'AGT', 'AGG']
allnonsyn['GAG']= ['AAA','AAG','GCT','GCC','GCA','GCG','AAC','AAT','TCT','TCA','TCG', 'TCC', 'AGT', 'AGG']
allnonsyn['TGT']= ['GCT','GCC','GCA','GCG','ATC','ATA','ATT','TTG','TTA','CTC','CTG', 'CTA','CTT','TTT','TTC']
allnonsyn['TGC']= ['GCT','GCC','GCA','GCG','ATC','ATA','ATT','TTG','TTA','CTC','CTG', 'CTA','CTT','TTT','TTC']
allnonsyn['TGG']= ['TAT','TAC','TTT','TTC','ACC','ACT','ACA','ACG','CAT','CAC']
allnonsyn['CGT']= ['AAA','AAG','CAG','CAA','CAT','CAC','GAA','GAG']
allnonsyn['CGC']= ['AAA','AAG','CAG','CAA','CAT','CAC','GAA','GAG']
allnonsyn['CGA']= ['AAA','AAG','CAG','CAA','CAT','CAC','GAA','GAG']
allnonsyn['CGG']= ['AAA','AAG','CAG','CAA','CAT','CAC','GAA','GAG']
allnonsyn['AGA']= ['AAA','AAG','CAG','CAA','CAT','CAC','GAA','GAG']
allnonsyn['AGG']= ['AAA','AAG','CAG','CAA','CAT','CAC','GAA','GAG']
allnonsyn['AGT']= ['ACC','ACT','ACA','ACG', 'GCT','GCC','GCA','GCG','AAC','AAT','GAC','GAT']
allnonsyn['AGG']= ['ACC','ACT','ACA','ACG', 'GCT','GCC','GCA','GCG','AAC','AAT','GAC','GAT']
allnonsyn['GGT']= ['GCT','GCC','GCA','GCG','AAC','AAT','TCT','TCA','TCG', 'TCC', 'AGT', 'AGG','GAC','GAT']
allnonsyn['GGC']= ['GCT','GCC','GCA','GCG','AAC','AAT','TCT','TCA','TCG', 'TCC', 'AGT', 'AGG','GAC','GAT']
allnonsyn['GGA']= ['GCT','GCC','GCA','GCG','AAC','AAT','TCT','TCA','TCG', 'TCC', 'AGT', 'AGG','GAC','GAT']
allnonsyn['GGG']= ['GCT','GCC','GCA','GCG','AAC','AAT','TCT','TCA','TCG', 'TCC', 'AGT', 'AGG','GAC','GAT']
allnonsyn['ATG']= ['TTG','TTA','CTC','CTG', 'CTA','CTT','GTC','GTA','GTG', 'GTT','ATC','ATA','ATT','CAG','CAA']

minspace=3
maxspace=7
rbs='AGGAGG'
lrbs=len(rbs)
atg='ATG'

numberofAAchangedautorized=0 # For now

def evaluateRBS(sequence,maxdist):
    global minspace
    global maxspace
    global rbs
    global lrbs
    global atg
    disttable=dict()
    for lspace in range(minspace,maxspace+1):
        drbs=[sequence[i]==x for (i,x) in enumerate(rbs)].count(False)
        datg=[sequence[i+lrbs+lspace]==x for (i,x) in enumerate(atg)].count(False)
        disttable[lspace]=(drbs+datg)
    return [(lspace,dist) for (lspace,dist) in disttable.items() if dist<=maxdist]

def treatRBS(genesequence,pos,rbssequence,lspace,distAA,save):
    global rbs
    global lrbs
    posfirstcds=pos+lrbs+lspace+3
    shift=posfirstcds%3
    if shift==0:
        # Then this RBS is not usefull because in frame
        return
    fr=posfirstcds-shift+3 # The position in the upstream gene where we will start considering making synonymous change to eliminate STOP codons in our downstream new frame
    # We do not start at posfirtscds-shift otherwise we would risk to modify the RBS we just created.
    assert len(genesequence[fr:])%3 == 0
    (modified_end_genesequence,remaining_stops)=tunestopfs.removestopinframepx(genesequence[fr:],shift)
    modified_genesequence=genesequence[:fr]+modified_end_genesequence
    assert( str(modified_genesequence.translate()) == str(genesequence.translate()) )
    sizeoverlap=(len(genesequence)-(pos+lrbs+lspace))
    protection=float(sizeoverlap)/len(genesequence)
    distBP=[modified_genesequence[i] == x for (i,x) in enumerate(genesequence)].count(False)
    save.writerow([pos,sizeoverlap,protection,shift,distBP,distAA,remaining_stops,str(modified_genesequence),'Antoine'])
    
def findRBS(genesequence,save):
    l=len(genesequence)
    local=18 # RBS (6) + spacer (<=7) + ATG (3) arrondi au codon près
    for pos in range(0,l-local): # For each position in the gene
        if pos<700:
            continue
        print pos
        candidate=genesequence[pos:pos+18] # We look at the 18 next nucleotides
        assert (len(candidate)==18)
        codonliste=[candidate[k:k+3].tostring() for k in range(0,local,3)] # We transform the list of 18 nucleotides into a list of 6 codons
        allpossible=[allsyn[k]+[k] for k in codonliste]
        for synonymous_combination in itertools.product(*allpossible): # We look at each possible combination of synonymous changes
            syn_candidate=Seq(reduce(lambda x,y: x+y,synonymous_combination))
            assert len(syn_candidate)==18
            assert( str(syn_candidate.translate()) == str(candidate.translate()) )
            # We calculate the score of this candidate. We have to be carreful that it can happens it is a good candidate for several different RBS configurations (different spacer sizes), so we want to find all the possible configurations
            allconfigs=evaluateRBS(str(syn_candidate),0)
            for (lspace,dist) in allconfigs:
                treatRBS(genesequence,pos,syn_candidate,lspace,dist,save)

save = csv.writer(open("findoverlap.csv", "wb"))
save.writerow(['Start position','Length of the Overlap','% protected','Frame','Number of bpchanged','Number of AA changed','Number of remaining STOPS''New seq'])

fasta_record = SeqIO.read(open("seqtest.fasta","r"),"fasta")
bioseq=fasta_record.seq
strseq=bioseq.tostring()

findRBS(bioseq,save)
