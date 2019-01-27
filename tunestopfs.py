from Bio import SeqIO
from Bio.Seq import Seq
from functools import reduce

computed=False

def isstop(c):
  if c=='TGA' or c=='TAG' or c=='TAA':
    return True
  else:
    return False

"""Measuring and tuning the presence and accessibility of stop codons"""

SynonymousCodons = { 
'TGT': ['TGC'],
'TGC': ['TGT'], 
'GAC': ['GAT'],
'GAT': ['GAC'], 
'TCT': ['TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
'TCG': ['TCT', 'TCA', 'TCC', 'AGC', 'AGT'],
'TCA': ['TCT', 'TCG', 'TCC', 'AGC', 'AGT'],
'TCC': ['TCT', 'TCG', 'TCA', 'AGC', 'AGT'],
'AGC': ['TCT', 'TCG', 'TCA', 'TCC', 'AGT'],
'AGT': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC'], 
'CAA': ['CAG'], 
'CAG': ['CAA'], 
'ATG': [], 
'AAC': ['AAT'],
'AAT': ['AAC'], 
'CCT': ['CCG', 'CCA', 'CCC'],
'CCG': ['CCT', 'CCA', 'CCC'],
'CCA': ['CCT', 'CCG', 'CCC'],
'CCC': ['CCT', 'CCG', 'CCA'], 
'AAG': ['AAA'],
'AAA': ['AAG'], 
'TAG': ['TGA', 'TAA'],
'TGA': ['TAG', 'TAA'],
'TAA': ['TAG', 'TGA'], 
'ACC': ['ACA', 'ACG', 'ACT'],
'ACA': ['ACC', 'ACG', 'ACT'],
'ACG': ['ACC', 'ACA', 'ACT'],
'ACT': ['ACC', 'ACA', 'ACG'],
'TTC': ['TTT'],
'TTT': ['TTC'], 
'GCA': ['GCC', 'GCG', 'GCT'],
'GCC': ['GCA', 'GCG', 'GCT'],
'GCG': ['GCA', 'GCC', 'GCT'],
'GCT': ['GCA', 'GCC', 'GCG'],
'GGT': ['GGG', 'GGA', 'GGC'],
'GGG': ['GGT', 'GGA', 'GGC'],
'GGA': ['GGT', 'GGG', 'GGC'],
'GGC': ['GGT', 'GGG', 'GGA'],
'ATC': ['ATA', 'ATT'],
'ATA': ['ATC', 'ATT'],
'ATT': ['ATC', 'ATA'], 
'TTA': ['TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
'TTG': ['TTA', 'CTC', 'CTT', 'CTG', 'CTA'],
'CTC': ['TTA', 'TTG', 'CTT', 'CTG', 'CTA'],
'CTT': ['TTA', 'TTG', 'CTC', 'CTG', 'CTA'],
'CTG': ['TTA', 'TTG', 'CTC', 'CTT', 'CTA'],
'CTA': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG'], 
'CAC': ['CAT'],
'CAT': ['CAC'], 
'CGA': ['CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
'CGC': ['CGA', 'CGG', 'CGT', 'AGG', 'AGA'],
'CGG': ['CGA', 'CGC', 'CGT', 'AGG', 'AGA'],
'CGT': ['CGA', 'CGC', 'CGG', 'AGG', 'AGA'],
'AGG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGA'],
'AGA': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG'],
'TGG': [], 
'GTA': ['GTC', 'GTG', 'GTT'],
'GTC': ['GTA', 'GTG', 'GTT'],
'GTG': ['GTA', 'GTC', 'GTT'],
'GTT': ['GTA', 'GTC', 'GTG'],
'GAG': ['GAA'],
'GAA': ['GAG'], 
'TAC': ['TAT'],
'TAT': ['TAC'] 
}

stopabilitybycodon=dict()
lessstopablesynonymous=dict()

def hamming(seq1,seq2):
    d=0
    for s1,s2 in zip(seq1,seq2):
        if s1!=s2:
            d=d+1
    return d

def stopability1(codon):
    """First possible measure of stopability at codon level: number of BPS that turn the codon in STOP"""
    s=0
    assert(len(codon)==3)
    if (hamming(codon,'TAA')==1):
        s=s+1
    if (hamming(codon,'TAG')==1):
        s=s+1
    if (hamming(codon,'TGA')==1):
        s=s+1
    return s

def stopability2(codon):
    """Second possible measure. 1 if there exist at least one BPS turning the codon in STOP, 0 otherwise"""
    s=0
    assert(len(codon)==3)
    if (hamming(codon,'TAA')==1):
        s=s+1
    elif (hamming(codon,'TAG')==1):
        s=s+1
    elif (hamming(codon,'TGA')==1):
        s=s+1
    return s
     
def calculatestuff():
    """Calculate for each codon of the genetic code the accessibility of stop codon by BPS
     and searching the synonymous codon that decreases the most this accessibility"""
    global computed
    if computed:
        return
    global stopabilitybycodon
    for codon in list(SynonymousCodons.keys()):
        stopabilitybycodon[codon]=stopability1(codon)
    for codon in list(SynonymousCodons.keys()):
        syn=SynonymousCodons[codon]
        syn.append(codon)
        candidate=min(syn,key=lambda x: stopabilitybycodon[x])
        if stopabilitybycodon[candidate]<stopabilitybycodon[codon]:
            lessstopablesynonymous[codon]=candidate
        else:
            lessstopablesynonymous[codon]=codon
    computed=True
    return


def accstop(s0):
    """This function calculate the percentages of BPS that will create a STOP codon in every frame of a given sequence"""
    """The sequence is taken as a python string"""
    f=0
    global stopabilitybycodon
    import math
    calculatestuff()
    l0=len(s0)
    stop0=0
    stop1=0
    stop2=0
    nbcodons=len(list(range(3,l0-5,3)))
    for i in range(3,l0-5,3):
        c0=s0[i:i+3]
        if isstop(c0):
            f=f+1
            continue
        c1=s0[i+1:i+4]
        c2=s0[i+2:i+5]
        stop0=stop0+stopabilitybycodon[c0]
        stop1=stop1+stopabilitybycodon[c1]
        stop2=stop2+stopabilitybycodon[c2]
    return(nbcodons-f,(stop0,stop1,stop2))

def facultativeaccstop(s0):
    """This function calculate, for a given coding sequence in main frame, 
    the part of accessiblity of stop codon that could have been avoided
     by using more robust synonymous codons."""
    """The sequence is taken as a python string"""
    f=0
    calculatestuff()
    global stopabilitybycodon
    global lessstopablesynonymous
    l0=len(s0)
    stop0=0
    nbcodons=len(list(range(3,l0-5,3)))
    for i in range(3,l0-5,3):
        c0=s0[i:i+3]
        if isstop(c0):
            f=f+1
            continue
        bestsyn=lessstopablesynonymous[c0]
        stop0=stop0+stopabilitybycodon[c0]-stopabilitybycodon[bestsyn]
    return(nbcodons-f,stop0)
    

def removestopinframepx(s0,x,verbose=True):
    """Take a sequence and remove stop codons in alternative frame with only synonymous changes in main frame"""
    """The sequence is taken as a BioSeq object and not as a python string"""
    """Returns the modified sequence and the number of remaining stops"""
    import itertools
    
    assert (x>0) and (x<3)
    prot0=s0.translate()
    l0=len(s0)
    
    sx=s0[x:l0-3+x]
    tx=sx.translate()
    pt=tx.find('*')
    
    bummer=False
    
    changedposition=[]
    
    # While there remains some stop codons
    while(pt>-1):
        # Find the involved codons in s0
        p1=pt*3
        c1=s0[p1:p1+3]
        c2=s0[p1+3:p1+6]
        # Try to change them with synonymous
        candfound=False
        l1=list(SynonymousCodons[str(c1)])
        l2=list(SynonymousCodons[str(c2)])
        l1.append(str(c1))
        l2.append(str(c2))
        allcand=dict()
        for cand1,cand2 in itertools.product(l1,l2):
            new=list(s0)
            new[p1:p1+3]=cand1[0:3]
            new[p1+3:p1+6]=cand2[0:3]
            news0=Seq("".join(new))
            # Check whether the stop codon has been removed in new sequence
            newsx=news0[x:l0-3+x]
            newtx=newsx.translate()
            newpt=newtx.find('*',pt)
            if (newpt>pt) or (newpt==-1):
                candfound=True
                # Calculate the number of BPS
                allcand[(cand1,cand2)]=[(c1+c2)[i]==BP for (i,BP) in enumerate(cand1+cand2)].count(False)
        if candfound: # We found at least one candidate
            # Find the best one (less BPS) among all possible candidates
            (cand1,cand2)=sorted(allcand,key=lambda u: allcand[u])[0]
            changedposition=changedposition+[p1+i for (i,BP) in enumerate(cand1+cand2) if BP!=(c1+c2)[i]]
            new=list(s0)
            new[p1:p1+3]=cand1[0:3]
            new[p1+3:p1+6]=cand2[0:3]
            news0=Seq("".join(new))
            newsx=news0[x:l0-3+x]
            newtx=newsx.translate()
            newpt=newtx.find('*',pt)
            assert(len(s0)==len(news0))
            if verbose:
                print("We removed stop codon at position " + str(pt*3+x) + ". ")
            s0=news0
            sx=newsx
            tx=newtx
            pt=newpt
        else:
            bummer=True
            if verbose:
                print("Bummer ! We are unable to remove stop codon at position " + str(pt*3+x) + ". Let's skip this one and go on.")
            pt=tx.find('*',pt+1)
            continue
    # We have removed all stop codons in frame x
    assert (str(s0.translate())==str(prot0))
    assert (bummer or (not '*' in tx))
    return (s0,changedposition,[i*3+x for (i,A) in enumerate(tx) if A=='*'])

def frameshiftability(sequence):
    """Measuring the susceptibility for frameshifts due to runs of repeated nucleotides"""
    """The sequence is taken as a python string"""
    n=2
    nbrepeats=[0 for nucl in sequence]
    for (i,nucl) in enumerate(sequence):
        if (i==0):
            nbrepeats[i]=1
            continue
        if (sequence[i-1]==sequence[i]):
            nbrepeats[i]=nbrepeats[i-1]+1
        else:
            nbrepeats[i]=1
    return nbrepeats

def frameshiftability_score(sequence):
    nbrepeats=frameshiftability(sequence)
    return sum(nbrepeats)-len(sequence)

def removeFShotspots(sequence):
    """Try to remove runs of 3 base pair or more doing only synonymous changes."""
    import itertools
    initsequence=sequence
    nbrepeats = frameshiftability(sequence)
    startrepeats = [i-1 for (i,x) in enumerate(nbrepeats) if x==2]
    for firstpos in startrepeats:
        nucl=sequence[firstpos]
        i=firstpos+1
        while (i<=len(sequence)-1 and sequence[i]==nucl):
            i=i+1
        lastpos=i-1
        codons = [sequence[k:k+3] for k in range(firstpos-firstpos%3,lastpos-lastpos%3+1,3)]
        prefixe=''
        suffixe=''
        if firstpos>=3:
            prefixe=sequence[firstpos-firstpos%3-3:firstpos-firstpos%3]
        if lastpos<=len(sequence)-4:
            suffixe=sequence[lastpos-lastpos%3+3:lastpos-lastpos%3+6]
        allpossible = [SynonymousCodons[k]+[k] for k in codons]
        allcombination = itertools.product(*allpossible)
        bestcombination = min(allcombination, key = lambda combination: frameshiftability_score(prefixe+reduce(lambda x,y: x+y,combination)+suffixe))
        l=list(sequence)
        l[firstpos-firstpos%3:lastpos-lastpos%3+3]=bestcombination
        sequence=str('').join(l)
    assert(str(Seq(initsequence).translate())==str(Seq(sequence).translate()))
    return(sequence,frameshiftability_score(sequence))

def optimizeoverlapfusion(sequence, frame):
    """Optimize the overlapped part of our construct according to the following criteria:
    Constraints:
    1) Do not create STOP codon in the frame of the downstream gene
    2) Do not create ATG in the frame of the downstream gene
    3) Only make synonymous changes in the frame of the upstream gene
    Optimization goal:
    1) STOP codons are less accessible by BPS in the frame of the upstream gene
    2) The sequence is more prone to create frame shifts
    Inputs:
    sequence is given in the frame of the upstream gene.
    frame is the shift necessary to get the frame of the downstream gene.
    """
    """The sequence is taken as a BioSeq object and not as a python string"""
    protup=sequence.translate()
    lseq=len(sequence)
    sshift=sequence[frame:lseq]
    protdown=sshift.translate()
    
