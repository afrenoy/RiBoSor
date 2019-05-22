""" Functions to remove start codons, stop codons, mono-nucleotide repeats and rare codons, using genetic code redundancy
These functions can be used in the main reading or in an alternative reading frame
"""

from Bio import SeqIO
from Bio.Seq import Seq
from functools import reduce

""" Data and simple functions to work with the genetic code
"""


NonSynonymousCodons = {
    'TTC': ['TAT', 'TAC', 'TGG', 'TGG', 'ATC', 'ATA', 'ATT', 'TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT'],
    'TTT': ['TAT', 'TAC', 'TGG', 'TGG', 'ATC', 'ATA', 'ATT', 'TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT'],
    'TTG': ['ATC', 'ATA', 'ATT', 'ATG', 'GTC', 'GTA', 'GTG', 'GTT', 'TTT', 'TTC'],
    'TTA': ['ATC', 'ATA', 'ATT', 'ATG', 'GTC', 'GTA', 'GTG', 'GTT', 'TTT', 'TTC'],
    'CTT': ['ATC', 'ATA', 'ATT', 'ATG', 'GTC', 'GTA', 'GTG', 'GTT', 'TTT', 'TTC'],
    'CTC': ['ATC', 'ATA', 'ATT', 'ATG', 'GTC', 'GTA', 'GTG', 'GTT', 'TTT', 'TTC'],
    'CTA': ['ATC', 'ATA', 'ATT', 'ATG', 'GTC', 'GTA', 'GTG', 'GTT', 'TTT', 'TTC'],
    'CTG': ['ATC', 'ATA', 'ATT', 'ATG', 'GTC', 'GTA', 'GTG', 'GTT', 'TTT', 'TTC'],
    'ATT': ['GTC', 'GTA', 'GTG', 'GTT', 'TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT', 'TTT', 'TTC', 'TGT', 'TGC'],
    'ATC': ['GTC', 'GTA', 'GTG', 'GTT', 'TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT', 'TTT', 'TTC', 'TGT', 'TGC'],
    'ATA': ['GTC', 'GTA', 'GTG', 'GTT', 'TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT', 'TTT', 'TTC', 'TGT', 'TGC'],
    'ATG': ['GTC', 'GTA', 'GTG', 'GTT', 'TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT', 'TTT', 'TTC', 'TGT', 'TGC'],
    'GTT': ['ATC', 'ATA', 'ATT', 'TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT', 'ATG', 'GCT', 'GCC', 'GCA', 'GCG'],
    'GTC': ['ATC', 'ATA', 'ATT', 'TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT', 'ATG', 'GCT', 'GCC', 'GCA', 'GCG'],
    'GTA': ['ATC', 'ATA', 'ATT', 'TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT', 'ATG', 'GCT', 'GCC', 'GCA', 'GCG'],
    'GTG': ['ATC', 'ATA', 'ATT', 'TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT', 'ATG', 'GCT', 'GCC', 'GCA', 'GCG'],
    'TCT': ['ACC', 'ACT', 'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'AAC', 'AAT', 'GAC', 'GAT'],
    'TCC': ['ACC', 'ACT', 'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'AAC', 'AAT', 'GAC', 'GAT'],
    'TCA': ['ACC', 'ACT', 'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'AAC', 'AAT', 'GAC', 'GAT'],
    'TCG': ['ACC', 'ACT', 'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'AAC', 'AAT', 'GAC', 'GAT'],
    'CCT': ['GCT', 'GCC', 'GCA', 'GCG', 'GAC', 'GAT', 'CAG', 'CAA', 'GAA', 'GAG'],
    'CCC': ['GCT', 'GCC', 'GCA', 'GCG', 'GAC', 'GAT', 'CAG', 'CAA', 'GAA', 'GAG'],
    'CCA': ['GCT', 'GCC', 'GCA', 'GCG', 'GAC', 'GAT', 'CAG', 'CAA', 'GAA', 'GAG'],
    'CCG': ['GCT', 'GCC', 'GCA', 'GCG', 'GAC', 'GAT', 'CAG', 'CAA', 'GAA', 'GAG'],
    'ACT': ['TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG', 'GTC', 'GTA', 'GTG', 'GTT', 'GAA', 'GAG', 'GGG', 'GGC', 'GGA', 'GGT'],
    'ACC': ['TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG', 'GTC', 'GTA', 'GTG', 'GTT', 'GAA', 'GAG', 'GGG', 'GGC', 'GGA', 'GGT'],
    'ACA': ['TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG', 'GTC', 'GTA', 'GTG', 'GTT', 'GAA', 'GAG', 'GGG', 'GGC', 'GGA', 'GGT'],
    'ACG': ['TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG', 'GTC', 'GTA', 'GTG', 'GTT', 'GAA', 'GAG', 'GGG', 'GGC', 'GGA', 'GGT'],
    'GCT': ['GTC', 'GTA', 'GTG', 'GTT', 'GGG', 'GGC', 'GGA', 'GGT', 'TGT', 'TGC', 'ACC', 'ACT', 'ACA', 'ACG'],
    'GCC': ['GTC', 'GTA', 'GTG', 'GTT', 'GGG', 'GGC', 'GGA', 'GGT', 'TGT', 'TGC', 'ACC', 'ACT', 'ACA', 'ACG'],
    'GCA': ['GTC', 'GTA', 'GTG', 'GTT', 'GGG', 'GGC', 'GGA', 'GGT', 'TGT', 'TGC', 'ACC', 'ACT', 'ACA', 'ACG'],
    'GCG': ['GTC', 'GTA', 'GTG', 'GTT', 'GGG', 'GGC', 'GGA', 'GGT', 'TGT', 'TGC', 'ACC', 'ACT', 'ACA', 'ACG'],
    'TAT': ['TGG', 'TGG', 'ATG', 'TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT', 'ATC', 'ATA', 'ATT'],
    'TAC': ['TGG', 'TGG', 'ATG', 'TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT', 'ATC', 'ATA', 'ATT'],
    'CAT': ['TAT', 'TAC', 'AGG', 'CGC', 'CGA', 'CGT', 'AGA', 'CGC', 'AAC', 'AAT', 'GAA', 'GAG'],
    'CAC': ['TAT', 'TAC', 'AGG', 'CGC', 'CGA', 'CGT', 'AGA', 'CGC', 'AAC', 'AAT', 'GAA', 'GAG'],
    'CAA': ['AGG', 'CGC', 'CGA', 'CGT', 'AGA', 'CGC', 'AAC', 'AAT', 'GAC', 'GAT', 'CAT', 'CAC'],
    'CAG': ['AGG', 'CGC', 'CGA', 'CGT', 'AGA', 'CGC', 'AAC', 'AAT', 'GAC', 'GAT', 'CAT', 'CAC'],
    'AAT': ['TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG', 'CAT', 'CAC', 'GAA', 'GAG', 'GGG', 'GGC', 'GGA', 'GGT'],
    'AAC': ['TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG', 'CAT', 'CAC', 'GAA', 'GAG', 'GGG', 'GGC', 'GGA', 'GGT'],
    'AAA': ['AGG', 'CGC', 'CGA', 'CGT', 'AGA', 'CGC', 'GAA', 'GAG', 'CAG', 'CAA', 'TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG'],
    'AAG': ['AGG', 'CGC', 'CGA', 'CGT', 'AGA', 'CGC', 'GAA', 'GAG', 'CAG', 'CAA', 'TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG'],
    'GAT': ['AAC', 'AAT', 'GAA', 'GAG', 'TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG', 'CAG', 'CAA'],
    'GAC': ['AAC', 'AAT', 'GAA', 'GAG', 'TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG', 'CAG', 'CAA'],
    'GAA': ['AAA', 'AAG', 'GCT', 'GCC', 'GCA', 'GCG', 'AAC', 'AAT', 'TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG'],
    'GAG': ['AAA', 'AAG', 'GCT', 'GCC', 'GCA', 'GCG', 'AAC', 'AAT', 'TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG'],
    'TGT': ['GCT', 'GCC', 'GCA', 'GCG', 'ATC', 'ATA', 'ATT', 'TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT', 'TTT', 'TTC'],
    'TGC': ['GCT', 'GCC', 'GCA', 'GCG', 'ATC', 'ATA', 'ATT', 'TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT', 'TTT', 'TTC'],
    'TGG': ['TAT', 'TAC', 'TTT', 'TTC', 'ACC', 'ACT', 'ACA', 'ACG', 'CAT', 'CAC'],
    'CGT': ['AAA', 'AAG', 'CAG', 'CAA', 'CAT', 'CAC', 'GAA', 'GAG'],
    'CGC': ['AAA', 'AAG', 'CAG', 'CAA', 'CAT', 'CAC', 'GAA', 'GAG'],
    'CGA': ['AAA', 'AAG', 'CAG', 'CAA', 'CAT', 'CAC', 'GAA', 'GAG'],
    'CGG': ['AAA', 'AAG', 'CAG', 'CAA', 'CAT', 'CAC', 'GAA', 'GAG'],
    'AGA': ['AAA', 'AAG', 'CAG', 'CAA', 'CAT', 'CAC', 'GAA', 'GAG'],
    'AGG': ['AAA', 'AAG', 'CAG', 'CAA', 'CAT', 'CAC', 'GAA', 'GAG'],
    'AGT': ['ACC', 'ACT', 'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'AAC', 'AAT', 'GAC', 'GAT'],
    'AGG': ['ACC', 'ACT', 'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'AAC', 'AAT', 'GAC', 'GAT'],
    'GGT': ['GCT', 'GCC', 'GCA', 'GCG', 'AAC', 'AAT', 'TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG', 'GAC', 'GAT'],
    'GGC': ['GCT', 'GCC', 'GCA', 'GCG', 'AAC', 'AAT', 'TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG', 'GAC', 'GAT'],
    'GGA': ['GCT', 'GCC', 'GCA', 'GCG', 'AAC', 'AAT', 'TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG', 'GAC', 'GAT'],
    'GGG': ['GCT', 'GCC', 'GCA', 'GCG', 'AAC', 'AAT', 'TCT', 'TCA', 'TCG', 'TCC', 'AGT', 'AGG', 'GAC', 'GAT'],
    'ATG': ['TTG', 'TTA', 'CTC', 'CTG', 'CTA', 'CTT', 'GTC', 'GTA', 'GTG', 'GTT', 'ATC', 'ATA', 'ATT', 'CAG', 'CAA']
}


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
    'GGT': ['GGC', 'GGG', 'GGA'],
    'GGG': ['GGC', 'GGT', 'GGA'],
    'GGA': ['GGC', 'GGT', 'GGG'],
    'GGC': ['GGT', 'GGG', 'GGA'],
    'ATC': ['ATT', 'ATA'],
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
    'CGA': ['CGT', 'CGC', 'CGG', 'AGA', 'AGG'],
    'CGC': ['CGT', 'CGG', 'CGA', 'AGA', 'AGG'],
    'CGG': ['CGT', 'CGC', 'CGA', 'AGA', 'AGG'],
    'CGT': ['CGC', 'CGG', 'CGA', 'AGA', 'AGG'],
    'AGG': ['CGT', 'CGC', 'CGG', 'CGA', 'AGA'],
    'AGA': ['CGT', 'CGC', 'CGG', 'CGA', 'AGG'],
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


CodonUsage = {
    'GGG': 10.99,
    'GGA': 7.92,
    'GGT': 24.85,
    'GGC': 29.44,
    'GAG': 17.79,
    'GAA': 39.58,
    'GAT': 32.22,
    'GAC': 19.05,
    'GTG': 26.25,
    'GTA': 10.88,
    'GTT': 18.38,
    'GTC': 15.22,
    'GCG': 33.66,
    'GCA': 20.28,
    'GCT': 15.34,
    'GCC': 25.51,
    'AGG': 1.22,
    'AGA': 2.05,
    'CGG': 5.38,
    'CGA': 3.53,
    'CGT': 21.02,
    'CGC': 22.02,
    'AAG': 10.21,
    'AAA': 33.62,
    'AAT': 17.62,
    'AAC': 21.67,
    'ATG': 27.77,
    'ATA': 4.28,
    'ATT': 30.40,
    'ATC': 25.00,
    'ACG': 14.37,
    'ACA': 7.02,
    'ACT': 8.92,
    'ACC': 23.38,
    'TGG': 15.31,
    'TGT': 5.18,
    'TGC': 6.44,
    'TAG': 0.23,
    'TAA': 2.02,
    'TGA': 0.90,
    'TAT': 16.32,
    'TAC': 12.27,
    'TTT': 22.40,
    'TTC': 16.59,
    'AGT': 8.71,
    'AGC': 16.03,
    'TCG': 8.92,
    'TCA': 7.13,
    'TCT': 8.50,
    'TCC': 8.59,
    'CAG': 28.84,
    'CAA': 15.45,
    'CAT': 12.90,
    'CAC': 9.72,
    'TTG': 13.72,
    'TTA': 13.89,
    'CTG': 52.82,
    'CTA': 3.85,
    'CTT': 11.04,
    'CTC': 11.04,
    'CCG': 23.27,
    'CCA': 8.52,
    'CCT': 7.04,
    'CCC': 5.52,
}


def compute_rarity_score(x):
    """Compute the rarity score of a list of codons
    The higer this score, the more rare codons are used
    """
    return sum([1 + (8. - CodonUsage[c]) / 8. for c in x if (CodonUsage[c] <= 8.)])


def generate_synonymous_combinations(codons_list):
    """Try all possible combinations of sets of codons,
    taking one codon of each set,
    and sort them according to the rarity score
    with the combinantion using the less rare codons first
    """
    import itertools
    score = dict()
    for combination in itertools.product(*[SynonymousCodons[i] + [i] for i in codons_list]):
        score[combination] = compute_rarity_score(combination)
    return sorted(score.items(), key=lambda x: x[1])


mfsc = dict()  # Most Frequent Synonymous Codon
for codon in CodonUsage:
    if len(SynonymousCodons[codon]) == 0:
        mfsc[codon] = codon
        continue
    candidate = max(SynonymousCodons[codon], key=lambda sc: CodonUsage[sc])
    if CodonUsage[candidate] < 8.0 and CodonUsage[codon] > CodonUsage[candidate]:
        mfsc[codon] = codon
    else:
        mfsc[codon] = candidate


def optimize(s0):
    """Codon optimization: replace each codon by the most common alternative"""
    newseq = list(s0)
    for i in range(0, len(s0), 3):
        newcodon = mfsc[s0[i:i + 3]]
        newseq[i:i + 3] = newcodon
    s1 = ''.join(newseq)
    return s1


def removestopinframepx(s0, x, verbose=True):
    """Take a sequence and remove stop codons in alternative frame with only synonymous changes in main frame
    The sequence is taken as a BioSeq object and not as a python string
    Returns the modified sequence and the number of remaining stops
    """

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
        for ((cand1, cand2), rarity_score) in generate_synonymous_combinations([c1, c2]):
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
    """Is this codon a start?"""
    return codon == 'ATG' or codon == 'GTG' or codon == 'TTG'


def findfirststart(sx, m=0):
    """Find the next start codon starting from position m"""
    sfx = codonfold(sx)
    assert(len(sfx) >= m)
    i = m
    while(i < len(sfx)):
        if isstart(sfx[i]):
            return i
        i = i + 1
    return -1


def countstart(sx):
    """How many start codons does the sequence contain?"""
    sfx = codonfold(sx)
    return len([c for c in sfx if isstart(c)])


def removestartinframepx(s0, x, verbose=True):
    """Take a sequence and remove start codons in alternative frame with only synonymous changes in main frame, without creating stop codons
    The sequence is taken as a python string and NOT a BioSeq object
    Returns the modified sequence and the number of remaining starts
    """

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
        for ((cand1, cand2), rarity_score) in generate_synonymous_combinations([c1, c2]):
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


def frameshiftability(sequence):
    """Measure the susceptibility for frameshifts due to runs of repeated nucleotides
    The sequence is taken as a python string
    """
    n = 2
    nbrepeats = [0 for nucl in sequence]
    for (i, nucl) in enumerate(sequence):
        if i == 0:
            nbrepeats[i] = 1
            continue
        if (sequence[i-1] == sequence[i]):
            nbrepeats[i] = nbrepeats[i-1]+1
        else:
            nbrepeats[i] = 1
    return nbrepeats


def frameshiftability_score(sequence):
    nbrepeats = frameshiftability(sequence)
    return sum(nbrepeats)-len(sequence)


def removeFShotspots2frame(sequence, frame, maxlrun, begconserve, endconserve):
    """Try to remove runs of 3 base pair or more doing only synonymous changes in main frame, and without creating start and stop codon in alternative frame.
    Input is a python string and not a BioSeq object
    """

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
    nbrepeats = frameshiftability(sequence)
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
        allcombination = map(lambda combination: (combination[1], reduce(lambda x, y: x + y, combination[0])), generate_synonymous_combinations(codons))

        # Only keep the combinations that do not create start/stop in alternative frame
        allowedcombination = filter(lambda combination: countstart((prefixe + combination[1] + suffixe)[frame:lcontext - 3 + frame]) <= nbstartsbefore and str(Seq((prefixe + combination[1] + suffixe)[frame:lcontext - 3 + frame]).translate()).count('*') <= nbstopsbefore, allcombination)

        # Sort them according to potentiallity for frame shifts and smaller use or rare codons
        bestcombination = min(allowedcombination, key=lambda combination: frameshiftability_score(prefixe + combination[1] + suffixe) * 100 + int(combination[0] * 10.))
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
    newnbrepeats = frameshiftability(sequence)
    remaininghotspots = [i + 1 - maxlrun for (i, x) in enumerate(newnbrepeats) if x == maxlrun]
    return(sequence, remaininghotspots)


def removerarecodonsinframepx(sequence, frame, maxlrun, rarethreshold=8., verbose=True):
    """Try to remove rare codons in alternative frame with the following constraints:
    Do not add starts in alternative frame
    Do not add stops in alternative frame
    Do not create fs hotspots
    Do not add rare codons in main frame
    """
    assert(type(sequence) == str)
    assert (frame > 0) and (frame < 3)
    protein = Seq(sequence).translate()
    l0 = len(sequence)
    sx = sequence[frame:l0 - 3 + frame]
    initnbstarts = countstart(sx)
    initnbstops = str(Seq(sx).translate()).count('*')
    initfs = frameshiftability(sequence)
    initfsscore = sum([x for x in initfs if x >= maxlrun])

    changedposition = []
    remaining = []

    # pt=findfirststart(sx)
    codons = [sx[3 * i:3 * i + 3] for i in range(0, int(len(sx) / 3))]
    listerare = [i for i in range(0, int(len(sx) / 3)) if CodonUsage[codons[i]] < rarethreshold]
    # sumscore=sum[CodonUsage[codons[i]] for i in listrare]
    nbrare = len(listerare)
    while(nbrare > 0):  # for irare in listerare:
        # print ' '
        irare = listerare[0]
        c1 = sequence[3 * irare:3 * irare + 3]
        c2 = sequence[3 * irare + 3:3 * irare + 6]
        # print codons[irare] + ' is rare ' + c1 + ' ' + c2
        candfound = False
        old_rarity_score = compute_rarity_score([c1, c2])  # Rarity score of the six bp sequence in main frame
        old_usage_framepx = CodonUsage[codons[irare]]  # Codon usage of the involved rare codon in alternative frame
        allcand = dict()
        for ((cand1, cand2), rarity_score) in generate_synonymous_combinations([c1, c2]):  # Warning: rarity_score is calculated in the frame of GalK.
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
            newfs = frameshiftability(newsequence)
            newfsscore = sum([x for x in newfs if x >= maxlrun])
            if newfsscore > initfsscore:
                # print 'no (fs)'
                continue
            if rarity_score > old_rarity_score:
                # print 'no (rare in main frame)'
                continue
            newcodons = [newsx[3 * i:3 * i + 3] for i in range(0, int(len(newsx) / 3))]
            newlisterare = [i for i in range(0, int(len(newsx) / 3)) if CodonUsage[newcodons[i]] < rarethreshold]
            diff = set(newlisterare) - set(remaining)
            if len(diff) > 0 and min(diff) < irare:  # Are we adding rare codons in candidate before the one we are trying to remove?
                # print 'no (rare in alt frame)'
                continue
            # This candidate is valid
            candfound = True
            altframecodon = newsx[3 * irare:3 * irare + 3]
            # Record the number of BPS and the rarity score in alternative frame
            allcand[(cand1, cand2)] = ([(c1 + c2)[i] == BP for (i, BP) in enumerate(cand1 + cand2)].count(False), CodonUsage[altframecodon])

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
            newlisterare = set([i for i in range(0, int(len(newsx) / 3)) if CodonUsage[newcodons[i]] < rarethreshold]) - set(remaining)

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


