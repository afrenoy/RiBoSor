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
    
