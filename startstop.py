from Bio import SeqIO
from Bio.Seq import Seq

"""Functions to remove start/stop codons in alternative frames"""

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
                print "We removed stop codon at position " + str(pt*3+x) + ". "
            s0=news0
            sx=newsx
            tx=newtx
            pt=newpt
        else:
            bummer=True
            if verbose:
                print "Bummer ! We are unable to remove stop codon at position " + str(pt*3+x) + ". Let's skip this one and go on."
            pt=tx.find('*',pt+1)
            continue
    # We have removed all stop codons in frame x
    assert (str(s0.translate())==str(prot0))
    assert (bummer or (not '*' in tx))
    return (s0,changedposition,[i*3+x for (i,A) in enumerate(tx) if A=='*'])

def codonfold(sx):
   assert (len(sx)%3==0)
   if len(sx)==0:
       return []
   if len(sx)==3:
      return [sx]
   return [sx[0:3]]+codonfold(sx[3:])

def isstart(codon):
   return codon=='ATG' or codon=='GTG' or codon=='TTG'

def firststart(sfx,off=0):
   if len(sfx)==0:
      return -1
   elif isstart(sfx[0]):
      return off
   else:
      return firststart(sfx[1:],off+1)

def findfirststart(sx, m=0):
   sfx=codonfold(sx)
   assert(len(sfx)>=m)
   cand=firststart(sfx[m:])
   if cand>-1:
      return cand+m
   else:
      return -1

def countstart(sx):
    sfx=codonfold(sx)
    return len([c for c in sfx if isstart(c)])


def removestartinframepx(s0,x,verbose=True):
    """Take a sequence and remove start codons in alternative frame with only synonymous changes in main frame, without creating stop codons"""
    """The sequence is taken as a python string and NOT a BioSeq object"""
    """Returns the modified sequence and the number of remaining starts"""
    import itertools
   
    assert(type(s0)==str)
    assert (x>0) and (x<3)
    prot0=Seq(s0).translate()
    l0=len(s0)
    sx=s0[x:l0-3+x]
    pt=findfirststart(sx)
    
    nstop=str(Seq(sx).translate()).count('*') 
    
    bummer=False
    
    changedposition=[]
    
    # While there remains some start codons
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
            news0="".join(new)
            # Check whether the start codon has been removed in new sequence
            newsx=news0[x:l0-3+x]
            newpt=findfirststart(newsx,pt)
            newnstop=str(Seq(newsx).translate()).count('*')
            if ((newpt>pt) or (newpt==-1)) and (newnstop<=nstop): # We removed one start and did not create one stop
                candfound=True
                # Calculate the number of BPS
                allcand[(cand1,cand2)]=[(c1+c2)[i]==BP for (i,BP) in enumerate(cand1+cand2)].count(False)
                if newnstop<nstop:
                    print "Weird ! While our stop-removing algorithm did not manage to remove a stop around position " + str(pt*3+x) + ", our  start-removing algorithm somehow manage to remove this stop."
                    nstop=newnstop
        if candfound: # We found at least one candidate
            # Find the best one (less BPS) among all possible candidates
            (cand1,cand2)=sorted(allcand,key=lambda u: allcand[u])[0]
            changedposition=changedposition+[p1+i for (i,BP) in enumerate(cand1+cand2) if BP!=(c1+c2)[i]]
            new=list(s0)
            new[p1:p1+3]=cand1[0:3]
            new[p1+3:p1+6]=cand2[0:3]
            news0="".join(new)
            newsx=news0[x:l0-3+x]
            newpt=findfirststart(newsx,pt)
            assert(len(s0)==len(news0))
            if verbose:
                print "We removed start codon at position " + str(pt*3+x) + ". "
            s0=news0
            sx=newsx
            pt=newpt
        else:
            bummer=True
            if verbose:
                print "Bummer ! We are unable to remove start codon at position " + str(pt*3+x) + ". Let's skip this one and go on."
            pt=findfirststart(sx,pt+1)
            continue
    # We have removed all start codons in frame x
    assert (str(Seq(s0).translate())==str(prot0))
    assert (str(Seq(sx).translate()).count('*')<=nstop) # Actually it is expected to be = and not <
    assert (bummer or findfirststart(sx)==-1)
    return (s0,changedposition,[i*3+x for (i,A) in enumerate(codonfold(sx)) if isstart(A)])


def removeFShotspots2frame(sequence,frame,maxlrun,verbose=True):
    """Try to remove runs of 3 base pair or more doing only synonymous changes in main frame, and without creating start and stop codon in alternative frame."""
    """Input is a python string and not a BioSeq object"""
    import itertools
    from sys import path 
    path.append('/Users/antoine/code/bioinfo')
    import tunestopfs
    
    # Check consistency of inputs
    assert(type(sequence)==str)
    assert (frame>0) and (frame<3)
    assert(maxlrun>1)

    # Record initial number of starts and stops in alternative frame (+frame)
    protein=Seq(sequence).translate()
    l0=len(sequence)
    sx=sequence[frame:l0-3+frame]
    initnbstarts=countstart(sx)
    initnbstops=str(Seq(sx).translate()).count('*') 
    
    # Find all runs
    initsequence=sequence
    nbrepeats = tunestopfs.frameshiftability(sequence)
    startrepeats = [i+1-maxlrun for (i,x) in enumerate(nbrepeats) if x==maxlrun]
    
    # Try to eliminate them 
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
        lcontext=len(suffixe)+len(prefixe)+len(codons)*3
        allpossible = [SynonymousCodons[k]+[k] for k in codons]
        allcombination = itertools.imap(lambda combination: reduce(lambda x,y: x+y,combination),itertools.product(*allpossible))
        allowedcombination = itertools.ifilter(lambda combination: countstart((prefixe+combination+suffixe)[frame:lcontext-3+frame])<=initnbstarts and str(Seq((prefixe+combination+suffixe)[frame:lcontext-3+frame]).translate()).count('*')<=initnbstops,allcombination)
        bestcombination = min(allowedcombination, key = lambda combination: tunestopfs.frameshiftability_score(prefixe+combination+suffixe))
        #print [x==sequence[firstpos-firstpos%3+i] for (i,x) in enumerate(bestcombination)].count(False)
        #could be use as a criteria to refine above metric: when equals in frameshiftability, we could choose the combination that minimizes the number of (synonymous) BPS made.
        l=list(sequence)
        l[firstpos-firstpos%3:lastpos-lastpos%3+3]=bestcombination
        sequence=str('').join(l)
    assert(str(Seq(initsequence).translate())==str(Seq(sequence).translate()))
    return(sequence,tunestopfs.frameshiftability_score(sequence))
