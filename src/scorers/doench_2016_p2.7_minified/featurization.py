import pandas
import numpy as np
import Bio.SeqUtils.MeltingTemp as Tm
import itertools


def featurize_data(data, learn_options):
    '''
    assumes that data contains the 30mer
    returns set of features from which one can make a kernel for each one
    '''

    feature_sets = {}

    if learn_options["nuc_features"]:
        # spectrum kernels (position-independent) and weighted degree kernels (position-dependent)
        get_all_order_nuc_features(data['30mer'], feature_sets, learn_options, learn_options["order"], max_index_to_use=30)

    if learn_options["gc_features"]:
        gc_above_10, gc_below_10, gc_count = gc_features(data)
        feature_sets['gc_above_10'] = pandas.DataFrame(gc_above_10)
        feature_sets['gc_below_10'] = pandas.DataFrame(gc_below_10)
        feature_sets['gc_count'] = pandas.DataFrame(gc_count)

    if learn_options["include_NGGX_interaction"]:
        feature_sets["NGGX"] = NGGX_interaction_feature(data)

    if learn_options["include_Tm"]:
        feature_sets["Tm"] = Tm_feature(data)

    check_feature_set_dimensions(feature_sets)

    return feature_sets


def check_feature_set_dimensions(feature_sets):
    '''
    Ensure the # of people is the same in each feature set
    '''
    N = None
    for ft in feature_sets.keys():
        N2 = feature_sets[ft].shape[0]
        if N is None:
            N = N2
        else:
            assert N == N2, "# of individuals do not match up across feature sets"


def NGGX_interaction_feature(data):
    '''
    assuming 30-mer, grab the NGGX _ _ positions, and make a one-hot
    encoding of the NX nucleotides yielding 4x4=16 features
    '''
    sequence = data['30mer'].values
    feat_NX = pandas.DataFrame()
    # check that GG is where we think
    for seq in sequence:
        if seq[25:27] != "GG":
            raise Exception("expected GG but found %s" % seq[25:27])
        NX = seq[24]+seq[27]
        NX_onehot = nucleotide_features(NX,order=2, include_pos_independent=False, max_index_to_use=2, prefix="NGGX")
        # NX_onehot[:] = np.random.rand(NX_onehot.shape[0]) ##TESTING RANDOM FEATURE
        feat_NX = pandas.concat([feat_NX, NX_onehot], axis=1, sort=True)
    return feat_NX.T


def get_all_order_nuc_features(data, feature_sets, learn_options, maxorder, max_index_to_use, prefix=""):
    for order in range(1, maxorder+1):
        #print "\t\tconstructing order %s features" % order
        nuc_features_pd, nuc_features_pi = apply_nucleotide_features(data, order, learn_options["num_proc"],
                                                                     include_pos_independent=True, max_index_to_use=max_index_to_use, prefix=prefix)
        feature_sets['%s_nuc_pd_Order%i' % (prefix, order)] = nuc_features_pd
        if learn_options['include_pi_nuc_feat']:
            feature_sets['%s_nuc_pi_Order%i' % (prefix, order)] = nuc_features_pi
        #print "\t\t\t\t\t\t\tdone"


def countGC(s):
    '''
    GC content for only the 20mer, as per the Doench paper/code
    '''
    assert len(s) == 30, "seems to assume 30mer"
    return len(s[5:25].replace('A', '').replace('T', ''))


def Tm_feature(data):
    '''
    assuming '30-mer'is a key
    get melting temperature features from:
        0-the 30-mer ("global Tm")
        1-the Tm (melting temperature) of the DNA:RNA hybrid from positions 16 - 20 of the sgRNA, i.e. the 5nts immediately proximal of the NGG PAM
        2-the Tm of the DNA:RNA hybrid from position 8 - 15 (i.e. 8 nt)
        3-the Tm of the DNA:RNA hybrid from position 3 - 7  (i.e. 5 nt)
    '''
    sequence = data['30mer'].values
    featarray = np.ones((sequence.shape[0],4))
    for i, seq in enumerate(sequence):
        if seq[25:27] != "GG":
            raise Exception("expected GG but found %s" % seq[25:27])
        rna = False
        featarray[i,0] = Tm.Tm_staluc(seq, rna=rna)        #30mer Tm
        featarray[i,1] = Tm.Tm_staluc(seq[20:25], rna=rna) #5nts immediately proximal of the NGG PAM
        featarray[i,2] = Tm.Tm_staluc(seq[12:20], rna=rna)   #8-mer
        featarray[i,3] = Tm.Tm_staluc(seq[7:12], rna=rna)      #5-mer

    feat = pandas.DataFrame(featarray, index=data.index, columns=["Tm global_%s" % rna, "5mer_end_%s" %rna, "8mer_middle_%s" %rna, "5mer_start_%s" %rna])

    return feat


def gc_features(data):
    gc_count = data['30mer'].apply(countGC)
    gc_count.name = 'GC count'
    gc_above_10 = (gc_count > 10)*1
    gc_above_10.name = 'GC > 10'
    gc_below_10 = (gc_count < 10)*1
    gc_below_10.name = 'GC < 10'
    return gc_above_10, gc_below_10, gc_count


def apply_nucleotide_features(seq_data_frame, order, num_proc, include_pos_independent, max_index_to_use, prefix=""):

    fast = True

    if fast:
        feat_pd = seq_data_frame.apply(nucleotide_features, args=(order, include_pos_independent, max_index_to_use, prefix, 'pos_dependent'))
        feat_pi = seq_data_frame.apply(nucleotide_features, args=(order, include_pos_independent, max_index_to_use, prefix, 'pos_independent'))
    else:
        # old, much slower code
        feat_pd = pandas.DataFrame()
        feat_pi = pandas.DataFrame()
        for s in seq_data_frame.values:
            assert not (s==''), "string is empty"
            feat_pd_tmp, feat_pi_tmp = nucleotide_features(s, order, include_pos_independent, max_index_to_use, prefix=prefix)
            feat_pd = pandas.concat([feat_pd,feat_pd_tmp], axis=1)
            feat_pi = pandas.concat([feat_pi,feat_pi_tmp], axis=1)
        
        feat_pd = feat_pd.T
        feat_pi = feat_pi.T  
    
    return feat_pd, feat_pi


def nucleotide_features(s, order, include_pos_independent, max_index_to_use, prefix="", feature_type='all'):
    '''
    compute position-specific order-mer features for the 4-letter alphabet
    (e.g. for a sequence of length 30, there are 30*4 single nucleotide features
          and (30-1)*4^2=464 double nucleotide features
    '''
    if max_index_to_use is not None:
        s = s[:max_index_to_use]
    #assert(len(s)==30, "length not 30")
    #s = s[:30] #cut-off at thirty to clean up extra data that they accidentally left in, and were instructed to ignore in this way
    raw_alphabet = ['A', 'T', 'C', 'G']
    alphabet = ["".join(i) for i in itertools.product(raw_alphabet, repeat=order)]
    features_pos_dependent = np.zeros(len(alphabet)*(len(s)-(order-1)))
    features_pos_independent = np.zeros(np.power(len(raw_alphabet),order))


    #for position in range(0, len(s)-order, 1): JENN 9/4/2014 failing when len(s)=2
    for position in range(0, len(s)-order+1, 1):
        nucl = s[position:position+order]
        features_pos_dependent[alphabet.index(nucl)+(position*len(alphabet))] = 1.0
        features_pos_independent[alphabet.index(nucl)] += 1.0
    index_dependent = ['%s_pd.Order%d_P%d' % (prefix, order,i) for i in range(len(features_pos_dependent))]

    if feature_type == 'all' or feature_type == 'pos_independent':
        if include_pos_independent:
            index_independent = ['%s_pi.Order%d_P%d' % (prefix, order,i) for i in range(len(features_pos_independent))]
            if feature_type == 'all':
                return pandas.Series(features_pos_dependent,index=index_dependent), pandas.Series(features_pos_independent,index=index_independent)
            else:
                return pandas.Series(features_pos_independent, index=index_independent)
    
    if np.any(np.isnan(features_pos_dependent)): raise Exception("found nan features in features_pos_dependent")
    if np.any(np.isnan(features_pos_independent)): raise Exception("found nan features in features_pos_independent")
    
    return pandas.Series(features_pos_dependent, index=index_dependent)
