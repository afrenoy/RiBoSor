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
len_rbs=len(rbs)
atg='ATG'

def evaluateRBS(sequence,maxdist):
    global minspace
    global maxspace
    global rbs
    global len_rbs
    global atg
    disttable=dict()
    relativefirstpostochange=dict()
    newsequence=dict()
    for len_spacer in range(minspace,maxspace+1):
        matchrbs=[sequence[i]==x for (i,x) in enumerate(rbs)]
        drbs=matchrbs.count(False)
        matchatg=[sequence[i+len_rbs+len_spacer]==x for (i,x) in enumerate(atg)]
        datg=matchatg.count(False)
        disttable[len_spacer]=(drbs+datg)
        if drbs>0:
            relativefirstpostochange[len_spacer]=matchrbs.index(False)
            newsequence[len_spacer]=list(sequence[:len_rbs+len_spacer+3])
            newsequence[len_spacer][relativefirstpostochange[len_spacer]]=rbs[relativefirstpostochange[len_spacer]]
        elif datg>0:
            relativefirstpostochange[len_spacer]=len_rbs+len_spacer+matchatg.index(False)
            newsequence[len_spacer]=list(sequence[:len_rbs+len_spacer+3])
            newsequence[len_spacer][relativefirstpostochange[len_spacer]]=atg[matchatg.index(False)]
        else:
            relativefirstpostochange[len_spacer]=-1
            newsequence[len_spacer]=list(sequence[:len_rbs+len_spacer+3])
    return [(len_spacer,"".join(newsequence[len_spacer]),disttable[len_spacer],relativefirstpostochange[len_spacer]) for len_spacer in disttable.keys() if disttable[len_spacer]<=maxdist]


def treatRBS(origgenesequence,pos,rbssequence,len_spacer,nb_nonsyn,pos_nonsyn,differ,save): # TODO: modify so it reports every change.
    # Compute the new sequence from the original sequence and the RBS-START
    global len_rbs
    posfirstcds=pos+len_rbs+len_spacer+3
    assert len(rbssequence)==len_rbs+len_spacer+3 # posfirstcds-pos
    genesequence=origgenesequence[:pos]+rbssequence+origgenesequence[posfirstcds:] # We do not use rbssequence entirely because we do not need the end
    shift=posfirstcds%3
    if shift==0:
        # Then this RBS is not usefull because in frame
        return
    fr=posfirstcds-shift+3 # The position in the upstream gene where we will start considering making synonymous change to eliminate STOP codons in our downstream new frame
    # We do not start at posfirtscds-shift otherwise we would risk to modify the RBS we just created.
    assert len(genesequence[fr:])%3 == 0
    
    # Remove the stop codon that could exist in our new frame (APH) without modifying what is coded by the existing gene in main frame (GalK)
    (modified_end_genesequence,changedposition,remaining_stops)=tunestopfs.removestopinframepx(genesequence[fr:],shift,False)
    modified_genesequence=genesequence[:fr]+modified_end_genesequence
    
    # In one text file per RBS, output the detail of every changed nucleotide
    file=open("findoverlap/"+str(pos)+".txt","w")
    
    # Starting with the possible non synonymous change made to create a perfect RBS-START
    if nb_nonsyn==1:
        posnonsyn=pos+pos_nonsyn
        origcodon=origgenesequence[posnonsyn:posnonsyn+3].translate()
        newcodon=genesequence[posnonsyn:posnonsyn+3].translate()
        file.write("At position " + str(posnonsyn) + " replaced " + str(origgenesequence[posnonsyn]) + " by " + str(genesequence[posnonsyn]) + " (non synonymous, " + str(origcodon) + " -> " + str(newcodon) +  ") to create a perfect RBS-START\n")
    
    # Output the synonymous changes we made to create this RBS-START
    for i in differ:
        file.write("At position " + str(i) + " replaced " + str(origgenesequence[i]) + " by " + str(genesequence[i]) + " (synonymous) to create a perfect RBS-START\n")
    
    # Output information about the synonymous changes we made to remove the stop codons
    for p in changedposition:
        file.write("At position " + str(fr+p) + " replaced " + str(origgenesequence[fr+p]) + " by " + str(modified_genesequence[fr+p]) + " (synonymous) to eliminate a STOP codon\n")
    
    # Output information about potential remaining STOP codons we were not able to remove
    for r in remaining_stops:
        file.write("At position " + str(r) + " unable to remove a STOP codon\n")
    
    # Compute other information we want about this overlap (number of syn/nonsyn changes, absolute and relative size, ...)
    distAA=[str(modified_genesequence.translate())[i] == x for (i,x) in enumerate(str(origgenesequence.translate()))].count(False)
    distBP=[modified_genesequence[i] == x for (i,x) in enumerate(origgenesequence)].count(False)
    sizeoverlap=(len(genesequence)-(pos+len_rbs+len_spacer))
    protection=float(sizeoverlap)/len(genesequence)
    save.writerow([pos,sizeoverlap,protection,shift,len_spacer,distBP,distAA,len(remaining_stops),'Antoine',str(modified_genesequence)])
    file.close()
    
def findRBS(genesequence,save):
    l=len(genesequence)
    local=18 # RBS (6) + spacer (<=7) + ATG (3) arrondi au codon près
    for pos in range(0,l-local): # For each position in the gene, we try to start an RBS there
        totreat=dict()
        shift=pos%3
        candidate=genesequence[pos-shift:pos-shift+18] # We look at the 18 next nucleotides
        assert (len(candidate)==18)
        codonliste=[candidate[k:k+3].tostring() for k in range(0,local,3)] # We transform the list of 18 nucleotides into a list of 6 codons
        allpossible=[allsyn[k]+[k] for k in codonliste]
        for synonymous_combination in itertools.product(*allpossible): # We look at each possible combination of synonymous changes
            syn_candidate=Seq(reduce(lambda x,y: x+y,synonymous_combination))
            differ=set([i+pos-shift for (i,x) in enumerate(candidate) if syn_candidate[i]!=x])
            rbs_candidate=syn_candidate[shift:]
            nbps=[candidate[i]==x for (i,x) in enumerate(syn_candidate)].count(False)
            assert (len(differ)==nbps)
            # We calculate the score of this candidate. We have to be carreful that it can happens it is a good candidate for several different RBS configurations (different spacer sizes), so we want to find all the possible configurations
            allconfigs=evaluateRBS(str(rbs_candidate), 1) # We allow RBS+spacer+ATG that have one nucleotide different from the consensus (can be turned into consensus making one non synonymous change)
            for (len_spacer,mod_rbs_candidate,nb_nonsyn,pos_nonsyn) in allconfigs: # nb_nonsyn is 0 or 1, so pos_nonsyn is only one number
                if not totreat.has_key(mod_rbs_candidate) or totreat[mod_rbs_candidate][1]>nb_nonsyn or ( (totreat[mod_rbs_candidate][1]==nb_nonsyn) and (totreat[mod_rbs_candidate][3]>nbps) ):
                    totreat[mod_rbs_candidate]=(len_spacer,nb_nonsyn,pos_nonsyn,nbps,differ)
                    # TODO: maybe index by rbs_candidate and not mod_rbs_candidate
        # Treat the best RBS we found
        if totreat:
            besttotreat=sorted(totreat, cmp=lambda x,y: (totreat[x][1]-totreat[y][1])*100 + totreat[x][3]-totreat[y][3])[0]
            treatRBS(genesequence,pos,besttotreat,totreat[besttotreat][0],totreat[besttotreat][1],totreat[besttotreat][2],totreat[besttotreat][4],save)

save = csv.writer(open("findoverlap.csv", "wb"))
save.writerow(['Start position','Length of the Overlap','% protected','Frame','Spacer','Number of base pair changed','Number of amino acid changed','Number of remaining STOPS','Algorithm','New sequence'])

fasta_record = SeqIO.read(open("seqtest.fasta","r"),"fasta")
bioseq=fasta_record.seq

os.mkdir('findoverlap')

findRBS(bioseq,save)
