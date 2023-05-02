#Calculates the Rule set 2 score for the given 30-mer
#Input: 1. 30mer sgRNA+context sequence, NNNN[sgRNA sequence]NGGNNN
#       2. Amino acid cut position, for full model prediction only
#       3. Percent peptide, for full model prediction only
#Output: Rule set 2 score

import sys
import pickle
import pandas
import argparse
import model_comparison


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seq',
        type=str,
        help='30-mer')
    return parser


# def calculate_doench_score(seq_list):
args = get_parser().parse_args()
seq = args.seq.upper()

model_file = 'saved_models/V3_model_nopos.pickle'

try:
    with open(model_file, 'rb') as f:
        model = pickle.load(f)
except:
    raise Exception("could not find model stored to file %s" % model_file)

# if seq[25:27] == 'GG':
scores = model_comparison.predict([seq], model=model)
print(scores)
# return scores
        # print 'Rule set 2 score: %.4f'% (score)
    # else:
    #     print >> sys.stderr, 'Calculates on-target scores for sgRNAs with NGG PAM only.'