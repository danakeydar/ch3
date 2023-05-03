

class Conf:

    dir_hg19 = '../res/hg19/'
    checkpoint_dir = ''
    numSurrounding = 400  # per side of CpG i.e. total is x2
    chrArr = [str(i) for i in range(1,23)]
    chrArr.extend(['X', 'Y'])
    suffix = ''

    ### YOUR SETTINGS - START ###

    filename_sequence = '../full/probe_to_surrounding_seq_one_hot_formatted.csv'
    filename_expression ='../full/expression_sampled.csv'
    filename_dist = '../full/distances_full_after_fix_26_2.csv'
    filename_labels = '../full/labels.csv'

    validation_portion_subjects = 0.1
    validation_portion_probes = 0.1
    train_portion_probes = 0.7

    ### YOUR SETTINGS - END ###

    # Below conf files are intended for use ONLY in dataProcessor, not in model code
    probeToSurroundingSeqFilePrefixAll = '../res/probe_to_surroundingSeq_'
    probeToSurroundingSeqFilePrefixChr = '../res/interims/probe_to_surroundingSeq_'
    probeToOneHotMtrxFilePrefixChr = '../res/probeToOneHotMtrx_'
    probeToOneHotMtrxFilePrefixAll = '../res/probeToOneHotMtrxAll'+str(suffix)
    probeToOneHotPrefixAll = '../res/probeToOneHotAll'+str(suffix)
    probeToOneHotPrefixChr = '../res/probeToOneHotChr_'+str(suffix)
    numBases = 5
    dfDistances = '../res/distances.csv'
    dfMethylName = 'combined_CH3'
    dfMethyl = '../res/BRCA_CA_normal_methyl.csv'
    # dfExpression = '../res/BRCA_CA_normal_expressi.csv'

    # changes
    dfExpression = ""

    numSampleInputCpgs = 4
    numInputCpgs = 5000

    epochs = 2
    batch_size = 32
    num_steps = 50000


class ConfSample(Conf):

    numSurrounding = 400 #per side
    suffix = ''

    filename_sequence = '../sampled/probe_to_surrounding_seq_one_hot_formatted_sampled_10.csv'
    filename_expression = '../sampled/expression_sampled_10.csv'
    filename_dist = '../sampled/distances_cpg_sampled_10.csv'
    filename_labels = '../sampled/labels_sampled_10.csv'

    #filename_sequence = '../old/probeToOneHotAll_sample_mini.csv'
    #filename_expression = '../old/e_sample_mini.csv'
    #filename_dist = '../old/d_sample_mini.csv'
    #filename_labels = '../old/ch3_sample_mini.csv'

    validation_portion_subjects = 0.5
    validation_portion_probes = 0.5
    train_portion_probes = 0.7

    probeToOneHotPrefixAll = './res/probeToOneHotAll_sample' + str(suffix)
    numBases = 5 #4
    dfDistances = './res/distances_sample_withDist_10k_closest_gene.csv'

    numSampleInputCpgs = 4

    epochs = 80
    batch_size = 7


sample = True
exp_n_genes = 17996
dist_n_genes = 17996

if sample:
    Conf = ConfSample
    exp_n_genes = 20157#20174
    dist_n_genes = 1823#1822