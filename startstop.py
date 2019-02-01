from Bio import SeqIO
from Bio.Seq import Seq

import codons as codonsfun
from functools import reduce
SynonymousCodons = codonsfun.SynonymousCodons
allnonsyn = codonsfun.NonSynonymousCodons


def optimize(s0):
    newseq = list(s0)
    for i in range(0, len(s0), 3):
        newcodon = codonsfun.mfsc[s0[i:i + 3]]
        newseq[i:i + 3] = newcodon
    s1 = ''.join(newseq)
    return s1


def removestopinframepx(s0, x, verbose=True):
    """Take a sequence and remove stop codons in alternative frame with only synonymous changes in main frame"""
    """The sequence is taken as a BioSeq object and not as a python string"""
    """Returns the modified sequence and the number of remaining stops"""
    import itertools

    # assert(type(s0)==str) #NON!!
    assert (x > 0) and (x < 3)
    prot0 = s0.translate()
    l0 = len(s0)

    sx = s0[x:l0 - 3 + x]
    tx = sx.translate()
    pt = tx.find('*')

    bummer = False

    changedposition = []

    # While there remains some stop codons
    while(pt > -1):
        # Find the involved codons in s0
        p1 = pt * 3
        c1 = str(s0[p1:p1 + 3])
        c2 = str(s0[p1 + 3:p1 + 6])
        # Try to change them with synonymous
        candfound = False
        allcand = dict()
        for ((cand1, cand2), rarity_score) in codonsfun.generate_synonymous_combinations([c1, c2]):
            new = list(s0)
            new[p1:p1 + 3] = cand1[0:3]
            new[p1 + 3:p1 + 6] = cand2[0:3]
            news0 = Seq("".join(new))
            # Check whether the stop codon has been removed in new sequence
            newsx = news0[x:l0 - 3 + x]
            newtx = newsx.translate()
            newpt = newtx.find('*', pt)
            if (newpt > pt) or (newpt == -1):
                candfound = True
                # Calculate the number of BPS
                allcand[(cand1, cand2)] = ([(c1 + c2)[i] == BP for (i, BP) in enumerate(cand1 + cand2)].count(False), rarity_score)
        if candfound:  # We found at least one candidate
            # Find the best one (less BPS and less rare codons) among all possible candidates
            (cand1, cand2) = sorted(allcand, key=lambda u: float(allcand[u][0]) + allcand[u][1] * 10.)[0]
            changedposition = changedposition + [p1 + i for (i, BP) in enumerate(cand1 + cand2) if BP != (c1 + c2)[i]]
            new = list(s0)
            new[p1:p1 + 3] = cand1[0:3]
            new[p1 + 3:p1 + 6] = cand2[0:3]
            news0 = Seq("".join(new))
            newsx = news0[x:l0 - 3 + x]
            newtx = newsx.translate()
            newpt = newtx.find('*', pt)
            assert(len(s0) == len(news0))
            if verbose:
                print("We removed stop codon at position " + str(pt * 3 + x) + ". ")
            s0 = news0
            sx = newsx
            tx = newtx
            pt = newpt
        else:
            bummer = True
            if verbose:
                print("Bummer ! We are unable to remove stop codon at position " + str(pt * 3 + x) + ". Let's skip this one and go on.")
            pt = tx.find('*', pt + 1)
            continue
    # We have removed all stop codons in frame x
    assert (str(s0.translate()) == str(prot0))
    assert (bummer or ('*' not in tx))
    return (s0, changedposition, [i * 3 + x for (i, A) in enumerate(tx) if A == '*'])


def codonfold(sx):
    """Transform a string into a list of codons"""
    assert (len(sx) % 3 == 0)
    return [sx[i:i + 3] for i in range(0, len(sx), 3)]


def isstart(codon):
    return codon == 'ATG' or codon == 'GTG' or codon == 'TTG'


def findfirststart(sx, m=0):
    sfx = codonfold(sx)
    assert(len(sfx) >= m)
    i = m
    while(i < len(sfx)):
        if isstart(sfx[i]):
            return i
        i = i + 1
    return -1

# def firststart(sfx,off=0):
#   if len(sfx)==0:
#      return -1
#   elif isstart(sfx[0]):
#      return off
#   else:
#      return firststart(sfx[1:],off+1)

# def findfirststart(sx, m=0):
#   sfx=codonfold(sx)
#   assert(len(sfx)>=m)
#   cand=firststart(sfx[m:])
#   if cand>-1:
#      return cand+m
#   else:
#      return -1


def countstart(sx):
    sfx = codonfold(sx)
    return len([c for c in sfx if isstart(c)])


def removestartinframepx(s0, x, verbose=True):
    """Take a sequence and remove start codons in alternative frame with only synonymous changes in main frame, without creating stop codons"""
    """The sequence is taken as a python string and NOT a BioSeq object"""
    """Returns the modified sequence and the number of remaining starts"""
    import itertools

    assert(type(s0) == str)
    assert (x > 0) and (x < 3)
    prot0 = Seq(s0).translate()
    l0 = len(s0)
    sx = s0[x:l0 - 3 + x]
    pt = findfirststart(sx)

    nstop = str(Seq(sx).translate()).count('*')

    bummer = False

    changedposition = []

    # While there remains some start codons
    while(pt > -1):
        # Find the involved codons in s0
        p1 = pt * 3
        c1 = s0[p1:p1 + 3]
        c2 = s0[p1 + 3:p1 + 6]
        # Try to change them with synonymous
        candfound = False
        allcand = dict()
        for ((cand1, cand2), rarity_score) in codonsfun.generate_synonymous_combinations([c1, c2]):
            new = list(s0)
            new[p1:p1 + 3] = cand1[0:3]
            new[p1 + 3:p1 + 6] = cand2[0:3]
            news0 = "".join(new)
            # Check whether the start codon has been removed in new sequence
            newsx = news0[x:l0 - 3 + x]
            newpt = findfirststart(newsx, pt)
            newnstop = str(Seq(newsx).translate()).count('*')
            if ((newpt > pt) or (newpt == -1)) and (newnstop <= nstop):  # We removed one start and did not create one stop
                candfound = True
                # Calculate the number of BPS
                allcand[(cand1, cand2)] = ([(c1 + c2)[i] == BP for (i, BP) in enumerate(cand1 + cand2)].count(False), rarity_score)
        if candfound:  # We found at least one candidate
            # Find the best one (less rare codons and less BPS) among all possible candidates
            (cand1, cand2) = sorted(allcand, key=lambda u: allcand[u][0] + int(allcand[u][1] * 10))[0]
            changedposition = changedposition + [p1 + i for (i, BP) in enumerate(cand1 + cand2) if BP != (c1 + c2)[i]]
            new = list(s0)
            new[p1:p1 + 3] = cand1[0:3]
            new[p1 + 3:p1 + 6] = cand2[0:3]
            news0 = "".join(new)
            newsx = news0[x:l0 - 3 + x]
            newpt = findfirststart(newsx, pt)
            newnstop = str(Seq(newsx).translate()).count('*')
            if newnstop < nstop:
                nstop = newnstop
            assert(len(s0) == len(news0))
            if verbose:
                print("We removed start codon at position " + str(pt * 3 + x) + ". ")
            s0 = news0
            sx = newsx
            pt = newpt
        else:
            bummer = True
            if verbose:
                print("Bummer ! We are unable to remove start codon at position " + str(pt * 3 + x) + ". Let's skip this one and go on.")
            pt = findfirststart(sx, pt + 1)
            continue
    # We have removed all start codons in frame x
    assert (str(Seq(s0).translate()) == str(prot0))
    assert (str(Seq(sx).translate()).count('*') <= nstop)  # Actually it is expected to be = and not <
    assert (bummer or findfirststart(sx) == -1)
    return (s0, changedposition, [i * 3 + x for (i, A) in enumerate(codonfold(sx)) if isstart(A)])


def removeFShotspots2frame(sequence, frame, maxlrun, begconserve, endconserve):
    """Try to remove runs of 3 base pair or more doing only synonymous changes in main frame, and without creating start and stop codon in alternative frame."""
    """Input is a python string and not a BioSeq object"""
    import itertools
    from sys import path
    path.append('/Users/antoine/code/bioinfo')
    import tunestopfs

    # Check consistency of inputs
    assert(type(sequence) == str)
    assert (frame > 0) and (frame < 3)
    assert(maxlrun > 1)
    assert(begconserve < endconserve)

    # Record initial number of starts and stops in alternative frame (+frame)
    protein = Seq(sequence).translate()
    l0 = len(sequence)
    sx = sequence[frame:l0 - 3 + frame]
    initnbstarts = countstart(sx)
    initnbstops = str(Seq(sx).translate()).count('*')

    # Find all runs
    initsequence = sequence
    nbrepeats = tunestopfs.frameshiftability(sequence)
    startrepeats = [i + 1 - maxlrun for (i, x) in enumerate(nbrepeats) if x == maxlrun]

    # Try to eliminate them
    for firstpos in startrepeats:

        # Find end position
        nucl = sequence[firstpos]
        i = firstpos + 1
        while (i <= len(sequence) - 1 and sequence[i] == nucl):
            i = i + 1
        lastpos = i - 1

        # If this run involves the part of the sequence that should be conserved (RBS-start) then we skip it
        if (firstpos - firstpos % 3 - 3 <= endconserve) and (lastpos - lastpos % 3 + 6 >= begconserve):
            continue

        # Find the involved codons in 'main' frame
        codons = [sequence[k:k + 3] for k in range(firstpos - firstpos % 3, lastpos - lastpos % 3 + 1, 3)]

        # Find the preceding and following codons because: they can also be used to avoid creating start/stop in alternative frame, and will be used to avoid creating a new frameshift hotspot.
        prefixe = ''
        suffixe = ''
        if firstpos >= 3:
            prefixe = sequence[firstpos - firstpos % 3 - 3:firstpos - firstpos % 3]
        if lastpos <= len(sequence) - 4:
            suffixe = sequence[lastpos - lastpos % 3 + 3:lastpos - lastpos % 3 + 6]

        # Record the number of start and stop in alternative frame of local sequence before making changes
        localseqbefore = prefixe + reduce(lambda x, y: x + y, codons) + suffixe
        lcontext = len(suffixe) + len(prefixe) + len(codons) * 3
        assert(len(localseqbefore) == lcontext)
        nbstartsbefore = countstart(localseqbefore[frame:lcontext - 3 + frame])
        nbstopsbefore = str(Seq(localseqbefore[frame:lcontext - 3 + frame]).translate()).count('*')

        # Generate all the possible combinations of synonymous codons and the corresponding local sequence
        allcombination = map(lambda combination: (combination[1], reduce(lambda x, y: x + y, combination[0])), codonsfun.generate_synonymous_combinations(codons))

        # Only keep the combinations that do not create start/stop in alternative frame
        allowedcombination = filter(lambda combination: countstart((prefixe + combination[1] + suffixe)[frame:lcontext - 3 + frame]) <= nbstartsbefore and str(Seq((prefixe + combination[1] + suffixe)[frame:lcontext - 3 + frame]).translate()).count('*') <= nbstopsbefore, allcombination)

        # Sort them according to potentiallity for frame shifts and smaller use or rare codons
        bestcombination = min(allowedcombination, key=lambda combination: tunestopfs.frameshiftability_score(prefixe + combination[1] + suffixe) * 100 + int(combination[0] * 10.))
        # could be use as a criteria to refine above metric: when equals in frameshiftability, we could choose the combination that minimizes the number of (synonymous) BPS made.
        # print [x==sequence[firstpos-firstpos%3+i] for (i,x) in enumerate(bestcombination)].count(False)

        # Take the best combination and replace appropriate nucleotides in sequence
        l = list(sequence)
        l[firstpos - firstpos % 3:lastpos - lastpos % 3 + 3] = bestcombination[1]
        sequence = str('').join(l)

    # Assert that we only made synonymous changes
    assert(len(sequence) == l0)
    assert(str(Seq(initsequence).translate()) == str(Seq(sequence).translate()))

    # Assert that we did not created start/stop in alternative frame
    finalnbstarts = countstart(sequence[frame:l0 - 3 + frame])
    finalnbstops = str(Seq(sequence[frame:l0 - 3 + frame]).translate()).count('*')
    assert(finalnbstarts <= initnbstarts)  # warning: not tolerant
    assert(finalnbstops <= initnbstops)  # warning: not tolerant

    # return the new sequence and the remaining hotspots
    newnbrepeats = tunestopfs.frameshiftability(sequence)
    remaininghotspots = [i + 1 - maxlrun for (i, x) in enumerate(newnbrepeats) if x == maxlrun]
    return(sequence, remaininghotspots)


def removerarecodonsinframepx(sequence, frame, maxlrun, rarethreshold=8., verbose=True):
    '''Try to remove rare codons in alternative frame with the following constraints:
        Do not add starts in alternative frame
        Do not add stops in alternative frame
        Do not create fs hotspots
        Do not add rare codons in main frame
        '''
    from sys import path
    import tunestopfs
    assert(type(sequence) == str)
    assert (frame > 0) and (frame < 3)
    protein = Seq(sequence).translate()
    l0 = len(sequence)
    sx = sequence[frame:l0 - 3 + frame]
    initnbstarts = countstart(sx)
    initnbstops = str(Seq(sx).translate()).count('*')
    initfs = tunestopfs.frameshiftability(sequence)
    initfsscore = sum([x for x in initfs if x >= maxlrun])

    changedposition = []
    remaining = []

    # pt=findfirststart(sx)
    codons = [sx[3 * i:3 * i + 3] for i in range(0, int(len(sx) / 3))]
    listerare = [i for i in range(0, int(len(sx) / 3)) if codonsfun.CodonUsage[codons[i]] < rarethreshold]
    # sumscore=sum[codonsfun.CodonUsage[codons[i]] for i in listrare]
    nbrare = len(listerare)
    while(nbrare > 0):  # for irare in listerare:
        # print ' '
        irare = listerare[0]
        c1 = sequence[3 * irare:3 * irare + 3]
        c2 = sequence[3 * irare + 3:3 * irare + 6]
        # print codons[irare] + ' is rare ' + c1 + ' ' + c2
        candfound = False
        old_rarity_score = codonsfun.compute_rarity_score([c1, c2])  # Rarity score of the six bp sequence in main frame
        old_usage_framepx = codonsfun.CodonUsage[codons[irare]]  # Codon usage of the involved rare codon in alternative frame
        allcand = dict()
        for ((cand1, cand2), rarity_score) in codonsfun.generate_synonymous_combinations([c1, c2]):  # Warning: rarity_score is calculated in the frame of GalK.
            # print cand1 + ' ' + cand2 + ' ' + str(irare)
            new = list(sequence)
            new[3 * irare:3 * irare + 3] = cand1[0:3]
            new[3 * irare + 3:3 * irare + 6] = cand2[0:3]
            newsequence = "".join(new)
            newsx = newsequence[frame:l0 - 3 + frame]
            # print newsequence + ' ' + newsx
            # Check that this candidate does not add start/stop codons or fs hotspots or rare codons in galK frame or rare codon upstream in alternative frame
            if countstart(newsx) > initnbstarts:
                # print 'no (start)'
                continue
            if str(Seq(newsx).translate()).count('*') > initnbstops:
                # print 'no (stop)'
                continue
            newfs = tunestopfs.frameshiftability(newsequence)
            newfsscore = sum([x for x in newfs if x >= maxlrun])
            if newfsscore > initfsscore:
                # print 'no (fs)'
                continue
            if rarity_score > old_rarity_score:
                # print 'no (rare in main frame)'
                continue
            newcodons = [newsx[3 * i:3 * i + 3] for i in range(0, int(len(newsx) / 3))]
            newlisterare = [i for i in range(0, int(len(newsx) / 3)) if codonsfun.CodonUsage[newcodons[i]] < rarethreshold]
            diff = set(newlisterare) - set(remaining)
            if len(diff) > 0 and min(diff) < irare:  # Are we adding rare codons in candidate before the one we are trying to remove?
                # print 'no (rare in alt frame)'
                continue
            # This candidate is valid
            candfound = True
            altframecodon = newsx[3 * irare:3 * irare + 3]
            # print 'yes ' + str(codonsfun.CodonUsage[altframecodon])
            # Record the number of BPS and the rarity score in alternative frame
            allcand[(cand1, cand2)] = ([(c1 + c2)[i] == BP for (i, BP) in enumerate(cand1 + cand2)].count(False), codonsfun.CodonUsage[altframecodon])

        if candfound:  # We found at least one candidate
            # Find the best one (less BPS and less rare codons) among all possible candidates
            sortedcandidates = sorted(allcand, key=lambda u: allcand[u][0] - int(allcand[u][1] * 10.))
            (cand1, cand2) = sortedcandidates[0]
            new_usage_framepx = allcand[sortedcandidates[0]][1]
            new = list(sequence)
            new[3 * irare:3 * irare + 3] = cand1[0:3]
            new[3 * irare + 3:3 * irare + 6] = cand2[0:3]
            newsequence = "".join(new)
            newsx = newsequence[frame:l0 - 3 + frame]
            assert(len(sequence) == len(newsequence))
            assert (str(Seq(sequence).translate()) == str(Seq(newsequence).translate()))

            newcodons = [newsx[3 * i:3 * i + 3] for i in range(0, int(len(newsx) / 3))]
            newlisterare = set([i for i in range(0, int(len(newsx) / 3)) if codonsfun.CodonUsage[newcodons[i]] < rarethreshold]) - set(remaining)

            # Is this candidate better than original sequence ?
            if new_usage_framepx > old_usage_framepx:
                sequence = newsequence
                codons = newcodons
                # listerare.remove(irare)
                listerare = sorted(list(newlisterare))
                nbrare = len(listerare)
                changedposition = changedposition + [3 * irare + i for (i, BP) in enumerate(cand1 + cand2) if BP != (c1 + c2)[i]]
                # print sequence + ' is new sequence'
                continue

        # No valid candidate found
        remaining = remaining + [irare]
        listerare.remove(irare)
        nbrare = len(listerare)
        if verbose:
            print("unable to remove rare codon " + codons[irare] + " at position " + str(3 * irare + frame))
    return (sequence, changedposition, [3 * x + frame for x in remaining])
