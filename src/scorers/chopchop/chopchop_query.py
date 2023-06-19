#!/usr/bin/env python2.7

import sys
import os
import ast
import argparse
import subprocess
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(prog="CHOPCHOP query")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--genePred_file",
                       help="Full path to genePred file from which gene names will be extracted.")
    group.add_argument("--gene_names",
                       help="You can supply comma separated list of gene names eg. rho,mtm2.")
    parser.add_argument("-o", "--outputDir", required=True,
                        help=("Output directory. It is important to not run this script without creating" +
                              " beforehand a folder for results. Cleanup after CHOPCHOP scripts will remove " +
                              "everything with the exception of .tsv files."))
    parser.add_argument("--amplican", default=None, type=int,
                        help=("This option will generate partial amplican configuration file " +
                              "containing guides, primers, amplicon sequences, identifiers and barcodes. " +
                              "amplican is an R package for anlysis of CRISPR experiments. Supply number " +
                              "of guideRNA to use per gene eg. --amplican 1"))
    parser.add_argument("--guide_dist", default=120, type=int,
                        help="Mimnimum distance between guides when using many guides per gene.")
    parser.add_argument("--bed", default=False, action='store_true',
                        help=("Creates bed file to be viewed in UCSC genome browser as custom tracks: " +
                              "https://genome.ucsc.edu/cgi-bin/hgCustom " +
                              "Each track is shown as Forward Primer --- Guide --- Reverse Primer. " +
                              "Location of the browser window is set for the first guide in the first gene."))
    args = parser.parse_known_args()

    if args[0].amplican is not None and args[0].amplican < 1:
        parser.error("Number of guides per gene in --amplican has to be greater than 0.")

    if args[0].bed and args[0].amplican is None:
        parser.error("Option --bed can only be used with option --amplican.")

    if args[0].bed and "--isoforms" in args[1]:
        parser.error("Option --bed is not yet ready with option --isoform.")

    if args[0].outputDir == "./":
        exit("Supply a folder for results.")

    if args[0].gene_names is not None:
        genes = args[0].gene_names.split(",")
    else:
        genePred = pd.read_csv(args[0].genePred_file, sep="\t")
        genes = genePred["name"].tolist()

    relative_cc = os.path.join(os.path.dirname(__file__), "chopchop.py")
    command = args[1]
    command.insert(0, relative_cc)
    command.append("-o")
    command.append(args[0].outputDir)

    if args[0].amplican is not None:
        if "-P" not in args[1] or "--makePrimers" not in args[1]:
            command.append("-P")
        if "-3" not in args[1] or "--primer3options" not in args[1]:
            command.append("-3")
            command.append(
                "PRODUCT_SIZE_MIN=150,PRODUCT_SIZE_MAX=290,PRIMER_MIN_SIZE=18,PRIMER_MAX_SIZE=25,PRIMER_OPT_SIZE=22,PRIMER_MIN_TM=57,PRIMER_MAX_TM=63,PRIMER_OPT_TM=60")

        config = pd.DataFrame(columns=(
        "ID", "Barcode", "Forward_Reads", "Reverse_Reads", "Group", "Control", "guideRNA", "Forward_Primer",
        "Reverse_Primer", "Direction", "Amplicon", "Donor"))
        previous_guide_positions = {}  # store for reference

    if args[0].bed:
        bed_file = os.path.join(args[0].outputDir, "UCSC_custom_track.bed")
        bed_handle = open(bed_file, "w")

    for gene in genes:

        gene_file = os.path.join(args[0].outputDir, gene + ".tsv")
        gene_handle = open(gene_file, "w")

        temp_comd = list(command)
        temp_comd.append("-Target")
        temp_comd.append(gene)

        c = subprocess.Popen(temp_comd, stdout=gene_handle, stderr=subprocess.STDOUT)
        c.wait()
        gene_handle.close()

        if args[0].amplican is not None:

            try:
                gene_table = pd.read_csv(gene_file, sep="\t")
                guide_loci = gene_table["Genomic location"].tolist()
                this_gene_chrom = guide_loci[0].split(":")[0]  # chrom of the guides
                guide_loci = [int(g_loc.split(":")[1]) for g_loc in guide_loci]  # genomic start position of the guides
            except Exception:
                print "No guides for " + gene
                continue

            if not previous_guide_positions.has_key(this_gene_chrom):
                previous_guide_positions[this_gene_chrom] = []

            guide_row_in_table = 0
            guide_num = 0

            while guide_num != args[0].amplican:

                if guide_row_in_table > len(guide_loci):
                    print "Only " + str(guide_num + 1) + " guides with specified distances for " + gene
                    break

                # this guide candidate has to be spaced with minimum distance to ALL previous guides
                distances = [abs(pre - guide_loci[guide_row_in_table]) > args[0].guide_dist for pre in
                             previous_guide_positions[this_gene_chrom]]
                if not all(distances):
                    guide_row_in_table += 1
                    continue

                # try getting primers
                try:
                    primer_file = os.path.join(args[0].outputDir, "primer_" + str(guide_row_in_table + 1) + ".json")
                    primer_handle = open(primer_file, "r")
                    primer_content = primer_handle.readlines()
                    primer_handle.close()
                    primer_content = ast.literal_eval(primer_content[0])
                except Exception:  # try next guide
                    guide_row_in_table += 1
                    continue

                # filter primers with offtargets
                primer_conent = [a for a in primer_content if a[9:12] == [0, 0, 0]]
                if not primer_content:
                    guide_row_in_table += 1
                    continue

                # get amplicon
                try:
                    amplicon = SeqIO.parse(
                        os.path.join(args[0].outputDir, gene + "_" + str(guide_row_in_table + 1) + ".gb"), "genbank")
                    for index, record in enumerate(amplicon):
                        amplicon = str(record.seq)
                except Exception:
                    guide_row_in_table += 1
                    continue

                # figure whether Left and RIGHT are swapped
                LR_swap = amplicon.lower().find(primer_content[0][7].lower()) == -1

                # select primers and barcodes toward minimizing number of unique barcodes
                # each barcode + F + R have to be unique and each barcode has to have unique set of fastq files
                barcodes_so_far = config["Barcode"].tolist()
                primer_left = ""
                primer_right = ""

                for br in barcodes_so_far:

                    if primer_left != "" and primer_right != "":
                        break

                    for pr in primer_content:
                        L = pr[8] if LR_swap else pr[7]
                        R = pr[7] if LR_swap else pr[8]
                        LR = L + R
                        this_br_rows = config.loc[config["Barcode"] == br]
                        this_br_pr = [a + b for a, b in zip(this_br_rows["Forward_Primer"].tolist(),
                                                            this_br_rows["Reverse_Primer"].tolist())]
                        if LR not in this_br_pr:
                            barcode = br
                            primer_left = L
                            primer_left_start = pr[1]
                            primer_right = R
                            primer_right_start = pr[3]
                            break

                # when none of the primers can't fit within existing barcodes make a new one
                # with the first primer from the list
                if primer_left == "" or primer_right == "":
                    primer_left = primer_conent[0][8] if LR_swap else primer_conent[0][7]
                    primer_left_start = primer_conent[0][3] if LR_swap else primer_conent[0][1]
                    primer_right = primer_conent[0][7] if LR_swap else primer_conent[0][8]
                    primer_right_start = primer_conent[0][1] if LR_swap else primer_conent[0][3]
                    barcode = 0 if not barcodes_so_far else max(barcodes_so_far) + 1

                # Direction of the guide
                direction = int(gene_table["Strand"][guide_row_in_table] == "-") if "--isoforms" not in args[1] else 0
                # trim amplicon sequence to the primers
                left_end = amplicon.lower().find(primer_left.lower())
                right_end = amplicon.lower().find(str(Seq(primer_right).reverse_complement()).lower()) + len(
                    primer_right)
                amplicon = amplicon[left_end:right_end].upper()

                # put to upper guide sequence in the amplicon
                guide = gene_table["Target sequence"][guide_row_in_table].upper()
                guide_left_pos = amplicon.find(guide) if direction == 0 else amplicon.find(
                    str(Seq(guide).reverse_complement()))
                guide_right_pos = guide_left_pos + len(guide)
                amplicon = amplicon[:guide_left_pos].lower() + amplicon[guide_left_pos:guide_right_pos] + amplicon[
                                                                                                          guide_right_pos:].lower()

                # add line to the config
                guide_name = gene + "_guide_" + str(guide_num + 1)
                config.loc[len(config.index)] = [guide_name,
                                                 barcode, "", "", gene, 0,
                                                 guide.upper(), primer_left.upper(), primer_right.upper(),
                                                 direction,
                                                 amplicon]
                previous_guide_positions[this_gene_chrom].append(guide_loci[guide_row_in_table])

                if args[0].bed:
                    primer_right_end = primer_right_start + len(primer_right)

                    if genes.index(
                            gene) == 0 and guide_num == 0:  # for the first gene only put the browser it's chromosome and add header
                        bed_handle.write("browser position {0}:{1}-{2}\n".format(this_gene_chrom, primer_left_start - 1,
                                                                                 primer_right_end - 1))
                        header = ("""track name=chopchop_query description=\"""" + (" ").join(sys.argv) +
                                  """\" visibility="pack" itemRgb="Off"\n""")
                        bed_handle.write(header)

                    bed_line = "{0}\t{1}\t{2}\t{3}\t0\t{4}\t{1}\t{2}\t0\t3\t{5}\t{6}\n".format(
                        this_gene_chrom, primer_left_start - 1, primer_right_end - 1, guide_name,
                        gene_table["Strand"][guide_row_in_table] if "--isoforms" not in args[1] else ".",
                                         str(len(primer_left)) + "," + str(len(guide)) + "," + str(len(primer_right)),
                                         "0," + str(guide_loci[guide_row_in_table] - primer_left_start) + "," + str(
                                             primer_right_start - primer_left_start))
                    bed_handle.write(bed_line)

                # before next while condition check
                guide_row_in_table += 1
                guide_num += 1

        offt_file = os.path.join(args[0].outputDir, "offtargetsTable.csv")
        if "--offtargetsTable" in args[1] and os.path.isfile(offt_file):
            os.rename(offt_file, os.path.join(args[0].outputDir, "%s_offtargets.tsv" % gene))

        print "Finished for " + gene

    if args[0].amplican is not None:
        config["Barcode"] = ["Barcode_" + str(int(i)) for i in config["Barcode"].tolist()]
        config["Control"] = config["Control"].astype(int)
        config["Direction"] = config["Direction"].astype(int)
        config.to_csv(os.path.join(args[0].outputDir, "amplican_config.tsv"), sep='\t')

    if args[0].bed:
        bed_handle.close()

    # cleanup after CHOPCHOP temp files
    for the_file in os.listdir(args[0].outputDir):
        file_path = os.path.join(args[0].outputDir, the_file)
        try:
            if os.path.isfile(file_path) and os.path.splitext(the_file)[1] not in [".tsv", ".bed", ".fastq"]:
                os.unlink(file_path)
        except Exception as e:
            print(e)


if __name__ == "__main__":
    main()
