#!/usr/bin/env python2.7
# -*- coding:UTF-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools
import csv
import os
import re
import getopt
import sys
import startstop
import codons as codonsfun

allsyn=codonsfun.SynonymousCodons
allnonsyn=codonsfun.NonSynonymousCodons

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


def treatRBS(origgenesequence,pos,rbssequence,len_spacer,nb_nonsyn,pos_nonsyn,differ,removestart,removefs,save,detaildir): # TODO: modify so it reports every change.
    # In one text file per RBS, output the detail of every changed nucleotide
    rbsfile=open(detaildir+"/"+str(pos)+".txt","w")
    rbsfile.write('-- Original sequence\n'+str(origgenesequence)+'\n')

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

    # Output the synonymous changes we made to create this RBS-START
    rbsfile.write('-- Creation of the RBS:'+'\n')
    for i in differ:
        rbsfile.write("  At position " + str(i) + " replaced " + str(origgenesequence[i]) + " by " + str(genesequence[i]) + " (synonymous) to create a perfect RBS-START\n")
    
    # Report the possible non synonymous change made to create a perfect RBS-START
    if nb_nonsyn==1:
        posnonsyn=pos+pos_nonsyn
        origcodon=origgenesequence[posnonsyn:posnonsyn+3].translate()
        newcodon=genesequence[posnonsyn:posnonsyn+3].translate()
        rbsfile.write("  At position " + str(posnonsyn) + " replaced " + str(origgenesequence[posnonsyn]) + " by " + str(genesequence[posnonsyn]) + " (non synonymous, " + str(origcodon) + " -> " + str(newcodon) +  ") to create a perfect RBS-START\n")
    
    # Make synonymous changes in galK sequence to try to get rid of this fucking cryptic promoter
    if fr>25:
        optimized_beg_galK=startstop.optimize(str(genesequence[:fr-24]))
        optimized_genesequence=Seq(optimized_beg_galK)+genesequence[fr-24:]
    optimized_end_galK=startstop.optimize(str(genesequence[fr+6:]))
    optimized_genesequence=optimized_genesequence[:fr+6]+Seq(optimized_end_galK)
    rbsfile.write('-- Sequence after optimizing the frame of galK\n'+str(optimized_genesequence)+'\n')

    # Remove the stop codons that could exist in our new frame (aph, after the ATG) without modifying what is coded by the existing gene in main frame (GalK)
    (end_genesequence_wostop,changedpositionstop,remaining_stops)=startstop.removestopinframepx(optimized_genesequence[fr:],shift,False)

    # Remove the start codons that could exist in our new frame (aph) without modifying what is coded by the existing gene in main frame (GalK)
    if (removestart==1): # We only remove starts after the ATG of aph
        (end_genesequence_wostopstart,changedpositionstart,remaining_starts)=startstop.removestartinframepx(str(end_genesequence_wostop),shift,False)
        modified_genesequence=optimized_genesequence[:fr]+Seq(end_genesequence_wostopstart)
    elif (removestart==2): # We also remove starts before the ATG of aph
        (end_genesequence_wostopstart,rel_changedpositionstart,rel_remaining_starts)=startstop.removestartinframepx(str(end_genesequence_wostop),shift,False)
        (start_genesequence_wostart,abs_changedpositionstart,abs_remaining_starts)=startstop.removestartinframepx(str(optimized_genesequence[:pos-pos%3]),shift,False)
        modified_genesequence=Seq(start_genesequence_wostart) + optimized_genesequence[pos-pos%3:fr] + Seq(end_genesequence_wostopstart)
        changedpositionstart=[x-fr for x in abs_changedpositionstart]+rel_changedpositionstart
        remaining_starts=[x-fr for x in abs_remaining_starts]+rel_remaining_starts
    else: # We do not remove starts at all 
        modified_genesequence=optimized_genesequence[:fr]+end_genesequence_wostop
        changedpositionstart=[]
        remaining_starts=[]

    
    # Output information about the synonymous changes we made to remove the stop codons
    rbsfile.write('-- Remove stop codons in new frame\n')
    for p in changedpositionstop:
        rbsfile.write("  At position " + str(fr+p) + " replaced " + str(optimized_genesequence[fr+p]) + " by " + str(modified_genesequence[fr+p]) + " (synonymous) to eliminate a STOP codon\n")
    
    # Output information about potential remaining STOP codons we were not able to remove
    for r in remaining_stops:
        rbsfile.write("At position " + str(fr+r) + " unable to remove a STOP codon\n")
    
    # Output information about the synonymous changes we made to remove the start codons
    rbsfile.write('-- Remove start codons in new frame\n')
    for p in changedpositionstart:
        rbsfile.write("  At position " + str(fr+p) + " replaced " + str(optimized_genesequence[fr+p]) + " by " + str(modified_genesequence[fr+p]) + " (synonymous) to eliminate a START codon\n")
    
    # Output information about potential remaining START codons we were not able to remove
    for r in remaining_starts:
        rbsfile.write("At position " + str(fr+r) + " unable to remove a START codon\n")
    
    if removefs:
        rbsfile.write('-- Remove FS hotspots (single nucleotide repeats)\n')
        # Remove FS hotspots (everywhere, not only in the overlapping sequence)
        (nofs_sequence, remaininghotspots) = startstop.removeFShotspots2frame(str(modified_genesequence),shift,4,pos,pos+len_rbs+len_spacer+3)
        for pn in [i for (i,x) in enumerate(str(modified_genesequence)) if not x==nofs_sequence[i]]:
            rbsfile.write("  At position " + str(pn) + " replaced " + str(modified_genesequence[pn]) + " by " + str(nofs_sequence[pn]) + " (synonymous) to eliminate a FS hotspot\n")

        # Output information about potential FS hotspots we were not able to remove
        for pn in remaininghotspots:
            rbsfile.write("At position " + str(pn) + " unable to remove a FS hotspot\n")

        modified_genesequence = Seq(nofs_sequence)
    
    # Remove rare codons in alternative frame
    (end_genesequence_worare,changedpositionrare,remaining_rare)=startstop.removerarecodonsinframepx(str(modified_genesequence)[fr:],shift,4,8.0,False)
    modified_genesequence=modified_genesequence[:fr]+Seq(end_genesequence_worare)
    # Output information about changes we made to remove these codons
    rbsfile.write('-- Remove rare codons\n')
    for p in changedpositionrare:
        rbsfile.write("  At position " + str(fr+p) + " replaced " + str(optimized_genesequence[fr+p]) + " by " + str(modified_genesequence[fr+p]) + " (synonymous) to eliminate a rare codon\n")
    for r in remaining_rare:
        rbsfile.write("At position " + str(fr+r) + " unable to remove a rare codon\n")

    rbsfile.write('-- Final sequence\n'+str(modified_genesequence)+'\n')

    # Compute other information we want about this overlap (number of syn/nonsyn changes, absolute and relative size, ...)
    distAA=[str(modified_genesequence.translate())[i] == x for (i,x) in enumerate(str(origgenesequence.translate()))].count(False)
    distBP=[modified_genesequence[i] == x for (i,x) in enumerate(origgenesequence)].count(False)
    sizeoverlap=(len(genesequence)-(pos+len_rbs+len_spacer))
    proportion=float(sizeoverlap)/len(genesequence)
    save.writerow([pos,sizeoverlap,proportion,shift,len_spacer,distBP,distAA,len(remaining_stops),'Antoine',str(modified_genesequence)])
    rbsfile.close()
    
def findRBS(genesequence,save,detaildir):
    l=len(genesequence)
    local=18 # RBS (6) + spacer (<=7) + ATG (3) arrondi au codon près
    for pos in range(0,l-local): # For each position in the gene, we try to start an RBS there
        totreat=dict()
        shift=pos%3
        candidate=genesequence[pos-shift:pos-shift+18] # We look at the 18 next nucleotides
        assert (len(candidate)==18)
        codonliste=[str(candidate[k:k+3]) for k in range(0,local,3)] # We transform the list of 18 nucleotides into a list of 6 codons
        allpossible=[allsyn[k]+[k] for k in codonliste]
        for (rarity_score,synonymous_combination) in codonsfun.smartcodonproduct(*allpossible): # We look at each possible combination of synonymous changes
            syn_candidate=Seq(reduce(lambda x,y: x+y,synonymous_combination))
            differ=set([i+pos-shift for (i,x) in enumerate(candidate) if syn_candidate[i]!=x])
            rbs_candidate=syn_candidate[shift:]
            nbps=[candidate[i]==x for (i,x) in enumerate(syn_candidate)].count(False)
            assert (len(differ)==nbps)
            # We calculate the score of this candidate. We have to be carreful that it can happens it is a good candidate for several different RBS configurations (different spacer sizes), so we want to find all the possible configurations
            allconfigs=evaluateRBS(str(rbs_candidate), 1) # We allow RBS+spacer+ATG that have one nucleotide different from the consensus (can be turned into consensus making one non synonymous change)
            for (len_spacer,mod_rbs_candidate,nb_nonsyn,pos_nonsyn) in allconfigs: # nb_nonsyn is 0 or 1, so pos_nonsyn is only one number
                if not totreat.has_key(mod_rbs_candidate) or totreat[mod_rbs_candidate][1]>nb_nonsyn or ( (totreat[mod_rbs_candidate][1]==nb_nonsyn) and (totreat[mod_rbs_candidate][3]>nbps) ):
                    totreat[mod_rbs_candidate]=(len_spacer,nb_nonsyn,pos_nonsyn,nbps,differ,rarity_score)
                    # TODO: maybe index by rbs_candidate and not mod_rbs_candidate
        # Treat the best RBS we found
        if totreat:
            besttotreat=sorted(totreat, cmp=lambda x,y: (totreat[x][1]-totreat[y][1])*1000 + int((totreat[x][5]-totreat[y][5])*30) + totreat[x][3]-totreat[y][3])[0]
            treatRBS(genesequence,pos,besttotreat,totreat[besttotreat][0],totreat[besttotreat][1],totreat[besttotreat][2],totreat[besttotreat][4],2,True,save,detaildir)

opts, args = getopt.getopt(sys.argv[1:], "", [])
name=args[0]
inputfilename=name+".fasta"
outputfilename=name+".csv"
detaildir=name

save = csv.writer(open(outputfilename, "wb"))
save.writerow(['Start position','Length of the Overlap','% protected','Frame','Spacer','Number of base pair changed','Number of amino acid changed','Number of remaining STOPS','Algorithm','New sequence'])

fasta_record = SeqIO.read(open(inputfilename,"r"),"fasta")
bioseq=fasta_record.seq

os.mkdir(detaildir)

findRBS(bioseq,save,detaildir)

