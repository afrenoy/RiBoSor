#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import argparse
from functools import reduce

def evaluateRBS(sequence, maxdist, minspace=3, maxspace=7, atg = 'ATG', rbs = 'AGGAGG', len_rbs=6):
    """ Evaluate whether a given sequence is a good RBS"""
    disttable = dict()
    relativefirstpostochange = dict()
    newsequence = dict()
    for len_spacer in range(minspace, maxspace+1):
        matchrbs = [sequence[i] == x for (i, x) in enumerate(rbs)]
        drbs = matchrbs.count(False)
        matchatg = [sequence[i+len_rbs+len_spacer] == x for (i, x) in enumerate(atg)]
        datg = matchatg.count(False)
        disttable[len_spacer] = (drbs+datg)
        if drbs > 0:
            relativefirstpostochange[len_spacer] = matchrbs.index(False)
            newsequence[len_spacer] = list(sequence[:len_rbs+len_spacer+3])
            newsequence[len_spacer][relativefirstpostochange[len_spacer]] = rbs[relativefirstpostochange[len_spacer]]
        elif datg > 0:
            relativefirstpostochange[len_spacer] = len_rbs + len_spacer+matchatg.index(False)
            newsequence[len_spacer] = list(sequence[:len_rbs+len_spacer+3])
            newsequence[len_spacer][relativefirstpostochange[len_spacer]] = atg[matchatg.index(False)]
        else:
            relativefirstpostochange[len_spacer] = -1
            newsequence[len_spacer] = list(sequence[:len_rbs+len_spacer+3])
    return [(len_spacer, "".join(newsequence[len_spacer]), disttable[len_spacer], relativefirstpostochange[len_spacer]) for len_spacer in list(disttable.keys()) if disttable[len_spacer] <= maxdist]


def treatRBS(origgenesequence, pos, rbssequence, len_spacer, nb_nonsyn, pos_nonsyn, differ, save, detaildir, removestart=2, removefs=True, optimize_beginning=True, len_rbs=6):
    """ Once a candidate RBS has been found, try to make it suitable for expression of a downstream protein
    This function reports all changes it makes in a text file
    These changes aim at: removing rare codons in both existing and new reading frames, removing stop codons downstream the new start in the new reading frame, removing the alternative start codons in the new reading frame, removing mono-nucleotide repeats (mutagenic motifs)
    """
    # In one text file per RBS, output the detail of every changed nucleotide
    if detaildir:
        rbsfile = open("%s/%d.txt" % (detaildir, pos), "w")
        print("-- Original sequence\n%s" % str(origgenesequence), file=rbsfile)
    else:
        rbsfile = None
    # Compute the new sequence from the original sequence and the RBS-START
    posfirstcds = pos+len_rbs+len_spacer+3
    assert len(rbssequence) == len_rbs+len_spacer+3  # posfirstcds-pos
    # We do not use rbssequence entirely because we do not need the end
    genesequence = origgenesequence[:pos] + rbssequence + origgenesequence[posfirstcds:]
    shift = posfirstcds % 3
    if shift == 0:
        # Then this RBS is not usefull because in frame
        return
    fr = posfirstcds-shift+3  # The position in the upstream gene where we will start considering making synonymous change to eliminate STOP codons in our downstream new frame
    # We do not start at posfirtscds-shift otherwise we would risk to modify the RBS we just created.
    if not (len(genesequence[fr:]) % 3 == 0):
        assert(len(origgenesequence) % 3 == len(genesequence) % 3)
        if rbsfile:
            print("Warning: the original sequence length is not a multiple of 3 (pseudogene?), cutting...", file=rbsfile)
        genesequence = genesequence[:len(genesequence) - len(genesequence) % 3]
        origgenesequence = origgenesequence[:len(origgenesequence) - len(origgenesequence) % 3]

    # Output the synonymous changes we made to create this RBS-START
    if rbsfile:
        print("-- Creation of the RBS:", file=rbsfile)
        for i in differ:
            print("  At position %d replaced %s by %s (synonymous) to create a perfect RBS-START"
                  % (i, origgenesequence[i], genesequence[i]), file=rbsfile)

    # Report the possible non synonymous change made to create a perfect RBS-START
    if nb_nonsyn == 1:
        posnonsyn = pos+pos_nonsyn
        origcodon = origgenesequence[posnonsyn:posnonsyn+3].translate()
        newcodon = genesequence[posnonsyn:posnonsyn+3].translate()
        if rbsfile:
            print("  At position %d replaced %s by %s (non synonymous, %s -> %s) to create a perfect RBS-START"
                  % (posnonsyn, origgenesequence[posnonsyn], genesequence[posnonsyn], origcodon, newcodon), file=rbsfile)

    # Make synonymous changes in existing gene to optimize expression by removing rare codons
    # In some cases, codon optimization is also a good way of removing cryptit promoters
    optimized_end_galK = readingframes.optimize(str(genesequence[fr+6:]))
    if optimize_beginning and (fr > 25):
        optimized_beg_galK = readingframes.optimize(str(genesequence[:fr-24]))
        optimized_genesequence = Seq(optimized_beg_galK)+genesequence[fr-24:]
        optimized_genesequence = optimized_genesequence[:fr+6]+Seq(optimized_end_galK)
    else:
        optimized_genesequence = genesequence[:fr+6]+Seq(optimized_end_galK)
    if rbsfile:
        print("-- Sequence after optimizing the frame of galK\n%s" % optimized_genesequence, file=rbsfile)

    # Remove the stop codons that could exist in our new frame (aph, after the ATG) without modifying what is coded by the existing gene in main frame (GalK)
    (end_genesequence_wostop, changedpositionstop, remaining_stops) = readingframes.removestopinframepx(optimized_genesequence[fr:], shift, False)

    # Remove the start codons that could exist in our new frame (aph) without modifying what is coded by the existing gene in main frame (GalK)
    if (removestart == 1):  # We only remove starts after the ATG of aph
        (end_genesequence_wostopstart, changedpositionstart, remaining_starts) = readingframes.removestartinframepx(str(end_genesequence_wostop), shift, False)
        modified_genesequence = optimized_genesequence[:fr]+Seq(end_genesequence_wostopstart)
    elif (removestart == 2):  # We also remove starts before the ATG of aph
        (end_genesequence_wostopstart, rel_changedpositionstart, rel_remaining_starts) = readingframes.removestartinframepx(str(end_genesequence_wostop), shift, False)
        (start_genesequence_wostart, abs_changedpositionstart, abs_remaining_starts) = readingframes.removestartinframepx(str(optimized_genesequence[:pos-pos % 3]), shift, False)
        modified_genesequence = Seq(start_genesequence_wostart) + optimized_genesequence[pos-pos % 3:fr] + Seq(end_genesequence_wostopstart)
        changedpositionstart = [x-fr for x in abs_changedpositionstart]+rel_changedpositionstart
        remaining_starts = [x-fr for x in abs_remaining_starts]+rel_remaining_starts
    else:  # We do not remove starts at all
        modified_genesequence = optimized_genesequence[:fr] + end_genesequence_wostop
        changedpositionstart = []
        remaining_starts = []

    if rbsfile:
        # Output information about the synonymous changes we made to remove the stop codons
        print('-- Remove stop codons in new frame', file=rbsfile)
        for p in changedpositionstop:
            # It is possible that one change is reverted by the next one if changing two close positions, in which case we do not need to report
            if str(optimized_genesequence[fr+p]) != str(modified_genesequence[fr+p]):
                print("  At position %d replaced %s by %s (synonymous) to eliminate a STOP codon"
                      % (fr+p, optimized_genesequence[fr+p], modified_genesequence[fr+p]), file=rbsfile)

        # Output information about potential remaining STOP codons we were not able to remove
        for r in remaining_stops:
            print("At position %d unable to remove a STOP codon" % (fr+r), file=rbsfile)

        # Output information about the synonymous changes we made to remove the start codons
        print("-- Remove start codons in new frame", file=rbsfile)
        for p in changedpositionstart:
            if str(optimized_genesequence[fr+p]) != str(modified_genesequence[fr+p]):
                print("  At position %d replaced %s by %s (synonymous) to eliminate a START codon"
                      % (fr+p, optimized_genesequence[fr+p], modified_genesequence[fr+p]), file=rbsfile)

        # Output information about potential remaining START codons we were not able to remove
        for r in remaining_starts:
            print("At position %s unable to remove a START codon" % (fr+r), file=rbsfile)

        if removefs:
            print("-- Remove FS hotspots (single nucleotide repeats)", file=rbsfile)
        # Remove FS hotspots (everywhere, not only in the overlapping sequence)
        (nofs_sequence, remaininghotspots) = readingframes.removeFShotspots2frame(str(modified_genesequence), shift, 4, pos, pos+len_rbs+len_spacer+3)
        if rbsfile:
            for pn in [i for (i, x) in enumerate(str(modified_genesequence)) if not x == nofs_sequence[i]]:
                if str(nofs_sequence[pn]) != str(modified_genesequence[pn]):
                    print("  At position %d replaced %s by %s (synonymous) to eliminate a FS hotspot"
                          % (pn, modified_genesequence[pn], nofs_sequence[pn]), file=rbsfile)

            # Output information about potential FS hotspots we were not able to remove
            for pn in remaininghotspots:
                print("At position %d unable to remove a FS hotspot" % pn, file=rbsfile)

        modified_genesequence = Seq(nofs_sequence)

    # Remove rare codons in alternative frame
    (end_genesequence_worare, changedpositionrare, remaining_rare) = readingframes.removerarecodonsinframepx(str(modified_genesequence)[fr:], shift, 4, 8.0, False)
    modified_genesequence = modified_genesequence[:fr] + Seq(end_genesequence_worare)
    if rbsfile:
        # Output information about changes we made to remove these codons
        print("-- Remove rare codons", file=rbsfile)
        for p in changedpositionrare:
            if str(optimized_genesequence[fr+p]) != str(modified_genesequence[fr+p]):
                print("  At position %d replaced %s by %s (synonymous) to eliminate a rare codon"
                      % (fr+p, optimized_genesequence[fr+p], modified_genesequence[fr+p]), file=rbsfile)
        for r in remaining_rare:
            print("At position %d unable to remove a rare codon" % (fr+r), file=rbsfile)

    if rbsfile:
        print("-- Final sequence\n%s" % modified_genesequence, file=rbsfile)
        rbsfile.close()

    # Compute other information we want about this overlap (number of syn/nonsyn changes, absolute and relative size, ...)
    distAA = [str(modified_genesequence.translate())[i] == x for (i, x) in enumerate(str(origgenesequence.translate()))].count(False)
    distBP = [modified_genesequence[i] == x for (i, x) in enumerate(origgenesequence)].count(False)
    sizeoverlap = (len(genesequence)-(pos+len_rbs+len_spacer))
    proportion = float(sizeoverlap)/len(genesequence)
    if rbsfile:
        print('%d,%d,%.18f,%d,%d,%d,%d,%d,%s,%s' % (pos, sizeoverlap, proportion, shift, len_spacer, distBP, distAA, len(remaining_stops), 'ff266af', str(modified_genesequence)), file=save)
    else:  # Quiet mode, do not output the full modified sequence
        print('%d,%d,%.18f,%d,%d,%d,%d,%d,%s' % (pos, sizeoverlap, proportion, shift, len_spacer, distBP, distAA, len(remaining_stops), 'ff266af'), file=save)


def findRBS(genesequence, save, detaildir=None):
    """ Try to create RBS in a coding sequence making only synonymous changes
    This is achieved by local bruteforce: at each position, we consider all possible synonymous substitutions in the next 18 nucleotides
    """
    lseq = len(genesequence)
    local = 18  # RBS (6) + spacer (<=7) + ATG (3) rounded to the next codon
    for pos in range(0, lseq-local):  # For each position in the gene, we try to start an RBS there
        totreat = dict()
        shift = pos % 3
        # We look at the 18 next nucleotides
        candidate = genesequence[pos-shift:pos-shift+18]
        assert (len(candidate) == 18)
        # We transform the list of 18 nucleotides into a list of 6 codons
        codonliste = [str(candidate[k:k+3]) for k in range(0, local, 3)]
        # We look at each possible combination of synonymous changes
        for (synonymous_combination, rarity_score) in readingframes.generate_synonymous_combinations(codonliste):
            syn_candidate = Seq(reduce(lambda x, y: x+y, synonymous_combination))
            differ = set([i+pos-shift for (i, x) in enumerate(candidate) if syn_candidate[i] != x])
            rbs_candidate = syn_candidate[shift:]
            nbps = [candidate[i] == x for (i, x) in enumerate(syn_candidate)].count(False)
            assert (len(differ) == nbps)
            # We calculate the score of this candidate. We have to be carreful that it can happens it is a good candidate for several different RBS configurations (different spacer sizes), so we want to find all the possible configurations
            # We allow RBS+spacer+ATG that have one nucleotide different from the consensus (can be turned into consensus making one non synonymous change)
            allconfigs = evaluateRBS(str(rbs_candidate), 1)
            # nb_nonsyn is 0 or 1, so pos_nonsyn is only one number
            for (len_spacer, mod_rbs_candidate, nb_nonsyn, pos_nonsyn) in allconfigs:
                if (
                    mod_rbs_candidate not in totreat
                    or totreat[mod_rbs_candidate][1] > nb_nonsyn
                    or ((totreat[mod_rbs_candidate][1] == nb_nonsyn) and (totreat[mod_rbs_candidate][3] > nbps))
                ):
                    totreat[mod_rbs_candidate] = (len_spacer, nb_nonsyn, pos_nonsyn, nbps, differ, rarity_score)
                    # TODO: maybe index by rbs_candidate and not mod_rbs_candidate
        # Treat the best RBS we found
        if totreat:
            besttotreat = sorted(totreat, key=lambda x: totreat[x][1]*1000 + int(totreat[x][5]*30) + totreat[x][3])[0]
            treatRBS(genesequence, pos, besttotreat, totreat[besttotreat][0], totreat[besttotreat][1], totreat[besttotreat][2], totreat[besttotreat][4], save, detaildir, removestart=2, removefs=True)


def main():
    parser = argparse.ArgumentParser(description="The RiBoSor creates alternative reading frames in an existing coding sequence, using genetic code redundancy")
    parser.add_argument('inputfile', type=str, help='input fasta file (containing one or several coding sequences)')
    parser.add_argument('--outputdir', type=str, help='output directory (will be inferred from the fasta file if not specified)', metavar="dir")
    parser.add_argument('--summary', action='store_true', default=False, help='only produce summary csv file (no detail by candidate), useful for screening')
    args = parser.parse_args()

    if not ".fasta" in args.inputfile:
        print("Warning: input should be a fasta file")

    if args.outputdir:
        outputdir = args.outputdir
    else:
        outputdir = os.path.splitext(args.inputfile)[0]
    os.makedirs(outputdir, exist_ok=True)

    records = list(SeqIO.parse(open(args.inputfile, "r"), "fasta"))
    for nseq, fasta_record in enumerate(records):
        bioseq = fasta_record.seq
        if len(records) > 1:
            seqname = fasta_record.name
        else:
            seqname = os.path.splitext(args.inputfile)[0]
        seqdir = os.path.join(outputdir, seqname)+'.details'
        if args.summary:
            seqdir = None
        else:
            os.mkdir(seqdir)
        outputfilename = os.path.join(outputdir,seqname)+'.csv'
        outputfile = open(outputfilename, "x")
        print("Start position, Length of the Overlap, Fraction protected, Frame, Spacer, Number of base pair changed, Number of amino acid changed, Number of remaining STOPS, Algorithm, New sequence", file=outputfile)
        try:
            findRBS(bioseq, outputfile, seqdir)
        except:
            print("Error while processing record %d of file %s" % (nseq, args.inputfile))

if __name__ == '__main__':
    import readingframes 
    #allsyn = readingframes.SynonymousCodons
    #allnonsyn = readingframes.NonSynonymousCodons
    main()
else:
    from ribosor import readingframes 
    #allsyn = readingframes.SynonymousCodons
    #allnonsyn = readingframes.NonSynonymousCodons
# This is rather ugly, but this permits to directly download and use this program as a script instead of installing it as a package
