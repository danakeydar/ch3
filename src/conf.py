class Conf:

    dir_hg19 = '../res/hg19/'
    checkpoint_dir = ''
    numSurrounding = 400  # per side of CpG i.e. total is x2
    chrArr = [str(i) for i in range(1,23)]
    chrArr.extend(['X', 'Y'])
    suffix = ''

    ### YOUR SETTINGS - START ###
    # filename_sequence = '../full/probeToOneHotAll_sample_mini.csv'
    # filename_expression = '../full/e_sample_mini.csv'
    # filename_dist = '../full/d_sample_mini.csv'
    # filename_labels = '../full/ch3_sample_mini.csv'

    # filename_sequence = '../full/probe_to_surrounding_seq_one_hot_formatted.csv'
    # filename_expression = '../full/expression.csv'
    # filename_dist = '../full/distances_final_1_4.csv'
    # filename_labels = '../full/labels.csv'  # lessDecimals.csv'

    filename_sequence = '../full/probe_to_surrounding_seq_one_hot_formatted.csv'
    filename_expression = '../full/expression.csv'
    filename_dist = '../full/distances_final_1_4_formatted.csv'
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
    dfExpression = '../res/BRCA_CA_normal_expressi.csv'

    numSampleInputCpgs = 4
    numInputCpgs = 5000

    epochs = 2
    batch_size = 32
    num_steps = 50000
    exp_n_genes = 0#17996
    dist_n_genes = 0#17996

class ConfSample(Conf):

    numSurrounding = 400 #per side
    suffix = ''
    PATH =  '../sampled/sampled_12_04/'
    filename_sequence = PATH + 'sampled_seq.csv'
    filename_expression = PATH + 'sampled_exp.csv'
    filename_dist = PATH + 'sampled_dist.csv'
    filename_labels = PATH + 'sampled_labels.csv'

    # filename_sequence = '../full/probeToOneHotAll_sample_mini.csv'
    # filename_expression = '../full/e_sample_mini.csv'
    # filename_dist = '../full/d_sample_mini.csv'
    # filename_labels = '../full/ch3_sample_mini.csv'

    #filename_sequence = '../sampled/sampled_seq.csv'
    #filename_expression = '../sampled/sampled_exp.csv'
    #filename_dist = '../sampled/sampled_dist.csv'
    #filename_labels = '../sampled/sampled_labels.csv'

    validation_portion_subjects = 0.5
    validation_portion_probes = 0.5
    train_portion_probes = 0.7

    probeToOneHotPrefixAll = '../res/probeToOneHotAll_sample' + str(suffix)
    numBases = 5 #4
    dfDistances = '../res/distances_sample_withDist_10k_closest_gene.csv'

    numSampleInputCpgs = 4

    epochs = 20
    batch_size = 7
    exp_n_genes = 0 #998
    dist_n_genes = 0 #1822

sample = True
if sample:
    Conf = ConfSample