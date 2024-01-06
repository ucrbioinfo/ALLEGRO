#!/usr/bin/env python2.7

import os
import argparse
import random
import subprocess
from Bio.Seq import Seq
from Bio.Restriction import Analysis, RestrictionBatch
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1:
            return
        yield start
        start += len(sub)  # use start += 1 to find overlapping matches


def calcGCContent(seq):
    Gcount = seq.count("G")
    Ccount = seq.count("C")
    GCcontent = 100 * (float(Gcount + Ccount) / len(seq))

    return GCcontent


def calcSelfComplementarity(fwd, backbone_regions, STEM_LEN=4):
    rvs = str(Seq(fwd).reverse_complement())
    L = len(fwd) - STEM_LEN - 1

    tmp = backbone_regions.strip().split(",")
    backbone_regions = [str(Seq(el).reverse_complement()) for el in tmp]

    folding = 0

    for i in range(0, len(fwd) - STEM_LEN):
        if calcGCContent(fwd[i:i + STEM_LEN]) >= 50:
            if fwd[i:i + STEM_LEN] in rvs[0:(L - i)] or any([fwd[i:i + STEM_LEN] in item for item in backbone_regions]):
                folding += 1

    return folding


def findRestrictionSites(sequence, restr_batch):
    mySeq = Seq(sequence, IUPACAmbiguousDNA())
    rb = RestrictionBatch(restr_batch)
    analyze = Analysis(rb, mySeq)

    return analyze.full()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genome", type=str,
                        help="Path to the bowtie index of the genome to search for off-targets.")
    parser.add_argument("--type", type=str, default="Cas9",
                        choices=["Cas9", "Cpf1"], help="Type of CRISPR.")
    parser.add_argument("--PAM", type=str, default="NGG",
                        help="The Protospacer Adjacent Motif where N is allowed eg. NGG")
    parser.add_argument("--g_len", type=int, default=20,
                        help="Length of the guide WITHOUT PAM eg. 20")
    parser.add_argument("--how_many", type=int, default=10,
                        help="Number of guides to find.")
    parser.add_argument("--GC_MAX", type=int, default=70,
                        help="Maximum acceptable value for GC percent eg. 70")
    parser.add_argument("--GC_MIN", type=int, default=40,
                        help="Minimum acceptable value for GC percent eg. 40")
    parser.add_argument("--backbone", type=str, default="",
                        help="Used when calculating complementarity of the guide to the backbone region. When empty self/backbone complementarity is not calculated. Can be comma separated string eg. AGGCTAGTCCGT,ATGCTGGAA")
    parser.add_argument("--restrict", type=str, default="",
                        help="Guides will be filtered for supplied comma separated list of restriction enzymes eg. BbsI,EcoRI or by default not filtered at all.")
    args = parser.parse_args()

    curr_dir = os.path.dirname(__file__)
    bowtie = os.path.join(curr_dir, "/bowtie/bowtie")

    allsequences = set({})
    print("Target sequence\tGC content (%)")

    while len(allsequences) < args.how_many:

        Ns = find_all(args.PAM, "N")
        PAM = list(args.PAM)
        for i in Ns:
            PAM[i] = "".join(random.choice("ACTG"))
        PAM = "".join(PAM)

        if args.type == "Cas9":
            seq = "".join(random.choice("ACTG") for i in range(args.g_len)) + PAM
        else:
            seq = PAM + "".join(random.choice("ACTG") for i in range(args.g_len))

        if seq in allsequences:
            continue

        gcperc = calcGCContent(seq)
        # validate on GC content
        if gcperc > args.GC_MAX or gcperc < args.GC_MIN:
            continue

        # bowtie against genome
        command = "%s -v 3 -m 50 -k 50 %s -c %s" % (bowtie, args.genome, seq)
        prog = subprocess.Popen(command, shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        result = prog.communicate()[0]
        if result != "":
            continue

        # Complementarity
        if args.backbone != "":
            if args.type == "Cas9":
                seq_noPAM = seq[len(PAM):]
            else:
                seq_noPAM = seq[:-len(PAM)]
            folding = calcSelfComplementarity(seq_noPAM, args.backbone)
            if folding > 0:
                continue

        # Restriction Sites
        if args.restrict != "":
            resSites = findRestrictionSites(seq, args.restrict.split(","))
            res = len(resSites.values()[0])
            if res > 0:
                continue

        allsequences.add(seq)
        print("%s\t%.0f" % (seq, gcperc))


if __name__ == '__main__':
    main()
