from Bio import SeqIO, SearchIO
import pandas as pd
import os
from os import path
from code.HmmPy import *
import time
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import random
sns.set_theme()
sns.set_style("whitegrid")
from code.HmmPy import *

cur = os.getcwd()
path_swiss_prot = path.join('data', 'part_1', 'swiss_prot',
                            'uniprot_sprot.fasta')


def swiss_prot_parser():
    if 'uniprot_sprot.csv' not in os.listdir(
            path.join('data', 'part_1', 'swiss_prot')):
        swissprot = list(SeqIO.parse(path_swiss_prot, "fasta"))
        swissprot = [[str(x.seq), x.id, x.description] for x in swissprot]

        columns_swissprot = ['seq', 'id', 'description']
        swissprot_df = pd.DataFrame(swissprot, columns=columns_swissprot)
        swissprot_df['ids'] = swissprot_df['id'].apply(
            lambda x: x.replace('|', '-').split('-')[1])
        swissprot_df.to_csv(cur + '\\data\\swiss_prot\\' + 'uniprot_sprot.csv')

    else:
        swissprot_df = pd.read_csv(
            path.join('data', 'part_1', 'swiss_prot', 'uniprot_sprot.csv'))

    return swissprot_df


def parse_psiblast():
    dirs_to_parse = [
        path.join('data', 'part_1', 'PSSMs', 'PSSM_{}'.format(a), 'to_parse')
        for a in ['C', 'M', 'O']
    ]
    #dirs_to_parse = [cur + '\\data\\PSSMs\\PSSM_' + a + '\\to_parse' for
    #a in ['C', 'M', 'O']]
    dirs_parsed = [
        path.join('data', 'part_1', 'PSSMs', 'PSSM_{}'.format(a), 'parsed')
        for a in ['C', 'M', 'O']
    ]
    #dirs_parsed = [cur + '\\data\\PSSMs\\PSSM_' + a + '\\parsed' for a in ['C', 'M', 'O']]

    parsed_dfs = {}

    for i, dir in enumerate(dirs_to_parse):
        files = os.listdir(dir)
        for filename in files:
            filename = filename.split('.')

            parsed_dfs['.'.join(filename)] = psiblast_parser(
                dir + '\\', filename[0], filename[1], dirs_parsed[i])

    return parsed_dfs


def psiblast_parser(dir_to_parse, filename, extension, dir_parsed):
    if filename + '.csv' not in os.listdir(dir_parsed):
        # print("Parsing {} ...".format(filename))
        hits_dict = {}
        blast_records = SearchIO.parse(
            path.join(dir_to_parse, filename + '.' + extension), 'blast-xml')
        for blast_record in blast_records:
            for rec in blast_record.hits:
                hits_dict.setdefault("ids", []).append(rec.id.split('|')[1])
                hits_dict.setdefault("e_value", []).append(rec.hsps[0].evalue)
                hits_dict.setdefault("hit_start",
                                     []).append(rec.hsps[0].hit_range[0] + 1)
                hits_dict.setdefault("hit_end",
                                     []).append(rec.hsps[0].hit_range[1])
                hits_dict.setdefault("bitscore",
                                     []).append(rec.hsps[0].bitscore)

        lengths = []
        with open(path.join(dir_to_parse, filename + '.' + extension),
                  'r') as f:
            for line in f:
                if (line[2:11] == '<Hit_len>'):
                    lengths.append(line.split('>')[1].split('<')[0])
        hits_dict["protein_length"] = lengths

        hits_df = pd.DataFrame.from_dict(hits_dict)
        hits_df.to_csv(path.join(dir_parsed, filename + '.csv'))

    else:
        hits_df = pd.read_csv(path.join(dir_parsed, filename + '.csv'))

    return hits_df


def metrics_sequences(df, gt):
    gt_acc = gt.accession.to_list()
    df = df.drop_duplicates(subset=['ids'])
    df_ids = df.ids.to_list()
    swissprot_df = swiss_prot_parser()
    len_swissprot = len(swissprot_df)
    len_df = len(df)

    # intersection between gt.ids and df.ids
    true_positives = (df['ids'].isin(gt_acc)).sum()
    # rows in swissprot_df but not in gt.ids and df.ids
    true_negatives = ((~swissprot_df['ids'].isin(gt_acc)) &
                      (~swissprot_df['ids'].isin(df_ids))).sum()
    false_positives = len_df - true_positives
    #swissprot_df['ids'].apply(lambda x: 1 if x not in df_ids else 0).sum() - true_negatives
    false_negatives = len_swissprot - len_df - true_negatives

    accuracy = (true_positives + true_negatives) / (
        true_positives + true_negatives + false_positives + false_negatives)

    precision = true_positives / (true_positives + false_positives)

    recall = true_positives / (true_positives + false_negatives)

    specificity = true_negatives / (true_negatives + false_positives)

    balanced_accuracy = (recall + specificity) / 2

    mcc = ((true_positives * true_negatives) -
           (false_positives * false_negatives)) / (np.sqrt(
               (true_positives + false_positives)) * np.sqrt(
                   (true_positives + false_negatives)) * np.sqrt(
                       (true_negatives + false_positives)) * np.sqrt(
                           (true_negatives + false_negatives)))

    f1_score = (2 * true_positives) / (2 * true_positives + false_positives +
                                       false_negatives)

    return [
        len_df, accuracy, precision, recall, specificity, balanced_accuracy,
        mcc, f1_score
    ]


def metrics_8(gt,
              h_threshold=False,
              hi_threshold=False,
              p_threshold=False,
              smart_update=True):
    """ 
    Compute all the various metrics for our models (point 8 of the project)
    - gt: dataframe containing ground truth proteins;
    - smart_update: if True, then the function doesn't recompute from scratch all the statistics,
    if those are already found in the metrics8.csv file.
    - h_threshold = hmm Sequence E-value;
    - hi_threshold = hmm i-Value;
    - p_threshold = PsiBlast E-Value;
    """
    parsed_tblouts, parsed_domtblouts = parse_hmms()
    parsed_psiblast = parse_psiblast()

    metrics = []
    index_metrics = list(parsed_domtblouts.keys()) + list(
        parsed_psiblast.keys())
    columns_metrics = [
        'n_hits', 'accuracy', 'precision', 'recall', 'specificity',
        'balanced_accuracy', 'mcc', 'f1_score'
    ]

    if h_threshold < 1 and h_threshold != False:
        h_threshold_name = str(h_threshold).split('.')[1]
    else:
        h_threshold_name = h_threshold

    if hi_threshold < 1 and hi_threshold != False:
        hi_threshold_name = str(hi_threshold).split('.')[1]
    else:
        hi_threshold_name = hi_threshold

    if p_threshold < 1 and p_threshold != False:
        p_threshold_name = str(p_threshold).split('.')[1]
    else:
        p_threshold_name = p_threshold

    # if 'metrics_8_{0}_{1}_{2}.csv'.format(h_threshold_name, hi_threshold_name,
    # p_threshold_name) in os.listdir(cur + '\\data\\metrics'):
    if 'metrics_8_{0}_{1}_{2}.csv'.format(h_threshold_name, hi_threshold_name,
                                          p_threshold_name) in os.listdir(
                                              path.join(
                                                  'data', 'part_1',
                                                  'metrics')):
        # old_metrics_df = pd.read_csv(cur +
        # '\\data\\metrics\\metrics_8_{0}_{1}_{2}.csv'.format(h_threshold_name,
        # hi_threshold_name, p_threshold_name), index_col=0)
        old_metrics_df = pd.read_csv(path.join(
            'data', 'part_1', 'metrics',
            'metrics_8_{0}_{1}_{2}.csv'.format(h_threshold_name,
                                               hi_threshold_name,
                                               p_threshold_name)),
                                     index_col=0)
    else:
        old_metrics_df = pd.DataFrame()

    for df in parsed_domtblouts.keys():
        if h_threshold:
            parsed_domtblouts[df] = parsed_domtblouts[df][
                parsed_domtblouts[df]['E-value'] < h_threshold]

        if hi_threshold:
            parsed_domtblouts[df] = parsed_domtblouts[df][
                parsed_domtblouts[df]['i-Evalue'] < hi_threshold]

        if smart_update and not old_metrics_df.empty:
            #If smartupdate = True, then compute only new entries
            if (df not in old_metrics_df.index.to_list()):
                # print("Computing metrics for: {}".format(df))
                metrics.append(metrics_sequences(parsed_domtblouts[df], gt))
            else:
                # print("Recycling metrics for {}".format(df))
                metrics.append(list(old_metrics_df.loc[df].values))
        else:
            # if smartupdate = False, we recompute from scratch all the metrics
            # print("Computing metrics for: {}".format(df))
            metrics.append(metrics_sequences(parsed_domtblouts[df], gt))

    for df in parsed_psiblast.keys():
        if p_threshold:
            parsed_psiblast[df] = parsed_psiblast[df][
                parsed_psiblast[df]['e_value'].astype(
                    dtype='float64') < p_threshold]

        if smart_update and not old_metrics_df.empty:
            #If smartupdate = True, then compute only new entries
            if (df not in old_metrics_df.index.to_list()):
                # print("Computing metrics for: {}".format(df))
                metrics.append(metrics_sequences(parsed_psiblast[df], gt))
            else:
                #  print("Recycling metrics for {}".format(df))
                metrics.append(list(old_metrics_df.loc[df].values))
        else:
            # if smartupdate = False, we recompute from scratch all the metrics
            # print("Computing metrics for: {}".format(df))
            metrics.append(metrics_sequences(parsed_psiblast[df], gt))

    metrics_df = pd.DataFrame(metrics, index_metrics, columns_metrics)

    metrics_df.to_csv(
        path.join(
            'data', 'part_1', 'metrics',
            'metrics_8_{0}_{1}_{2}.csv'.format(h_threshold_name,
                                               hi_threshold_name,
                                               p_threshold_name)))

    return metrics_df, parsed_tblouts, parsed_domtblouts, parsed_psiblast


def dfs_to_dicts(parsed_domtblouts, parsed_psiblast):
    # from dfs to dicts with no repeated columns -> fast af
    list_dfs = list(parsed_domtblouts.keys()) + list(parsed_psiblast.keys())
    list_dfs_not_repeated = []

    for df in parsed_domtblouts.keys():
        df_not_repeated = {}
        for i, row in parsed_domtblouts[df].iterrows():
            # if 'domtblout' in df:
            df_not_repeated.setdefault(row['ids'], []).append([
                row['from_ali_coord'], row['to_ali_coord'], row['target_len']
            ])
            # else:
            #     df_not_repeated.setdefault(row['ids'],[]).append([row['hit_start'], row['hit_end'], row['protein_length']])

        list_dfs_not_repeated.append(df_not_repeated)

    for df in parsed_psiblast.keys():
        df_not_repeated = {}
        for i, row in parsed_psiblast[df].iterrows():
            # if 'domtblout' in df:
            #     df_not_repeated.setdefault(row['ids'],[]).append([row['from_ali_coord'], row['to_ali_coord'], row['target_len']])
            # else:
            df_not_repeated.setdefault(row['ids'], []).append(
                [row['hit_start'], row['hit_end'], row['protein_length']])

        list_dfs_not_repeated.append(df_not_repeated)

    dict_dfs = dict(zip(list_dfs, list_dfs_not_repeated))

    return dict_dfs


def new_logic_case_1(a, b):  #, tp, tn, fp, fn):
    global tp
    global tn
    global fp
    global fn

    if a == 0 and b == 1:
        fp += 1

    elif a == 1 and b == 1:
        tp += 1

    elif a == 1 and b == 0:
        fn += 1

    else:
        tn += 1


def new_logic_case_2(a, b):  #, tp_, tn_, fp_, fn_):
    global tp_
    global tn_
    global fp_
    global fn_

    if a == 0 and b == 1:
        fp_ += 1

    elif a == 1 and b == 1:
        tp_ += 1

    elif a == 1 and b == 0:
        fn_ += 1

    else:
        tn_ += 1


def compute_con_matrix_9(parsed_domtblouts, parsed_psiblast, gt):
    new_logic_case_1_vectorized = np.vectorize(new_logic_case_1)
    new_logic_case_2_vectorized = np.vectorize(new_logic_case_2)

    gt_acc = gt.accession.to_list()

    dict_dfs = dfs_to_dicts(parsed_domtblouts, parsed_psiblast)

    conf_matrix = []

    #conf_df = pd.DataFrame(columns=['true_positives', 'true_negatives', 'false_positives', 'false_negatives'])

    for k, (df_name, df_dict) in enumerate(dict_dfs.items()):
        gt_int_df = [x for x in gt_acc
                     if x in list(df_dict.keys())]  # ids both in gt and df
        df_not_gt = [x for x in list(df_dict.keys())
                     if x not in gt_acc]  # ids in df and not in gt
        gt_not_df = [x for x in gt_acc if x not in list(df_dict.keys())
                     ]  # ids in gt and not in df

        # this is the only way to declare the variables and let new_logic function modify them
        global tp
        global tn
        global fp
        global fn
        global tp_
        global tn_
        global fp_
        global fn_
        tp = 0
        tn = 0
        fp = 0
        fn = 0
        tp_ = 0
        tn_ = 0
        fp_ = 0
        fn_ = 0

        for j, (ids, lists) in enumerate(df_dict.items()):
            # case 1 of Riccà's schema
            if ids in gt_int_df:
                row_df = np.zeros((int(lists[0][2]), ), dtype=int)

                for i, el in enumerate(lists):
                    row_df[int(el[0]):int(el[1])] = 1

                row_gt = np.zeros((int(lists[0][2]), ), dtype=int)
                row_gt[gt[gt.accession == ids].iloc[0, 1]:gt[
                    gt.accession == ids].iloc[0, 2]] = 1

                new_logic_case_1_vectorized(row_gt, row_df)

            # case 2 of Riccà's schema
            elif ids in df_not_gt:
                row_df = np.zeros((int(lists[0][2]), ), dtype=int)

                for i, el in enumerate(lists):
                    row_df[int(el[0]):int(el[1])] = 1

                row_gt = np.zeros((int(lists[0][2]), ), dtype=int)

                new_logic_case_2_vectorized(row_gt, row_df)

        # case 3 of Riccà's schema
        if gt_not_df != []:
            for ids in gt_not_df:
                row_gt = np.zeros((gt[gt.accession == ids].iloc[0, 3], ),
                                  dtype=int)
                row_gt[gt[gt.accession == ids].iloc[0, 1]:gt[
                    gt.accession == ids].iloc[0, 2]] = 1

                row_df = np.zeros((gt[gt.accession == ids].iloc[0, 3], ),
                                  dtype=int)

                new_logic_case_1_vectorized(row_gt, row_df)  #, tp, tn, fp, fn)

        conf_matrix.append([tp + tp_, tn + tn_, fp + fp_, fn + fn_])

    conf_df = pd.DataFrame(conf_matrix,
                           index=list(dict_dfs.keys()),
                           columns=[
                               'true_positives', 'true_negatives',
                               'false_positives', 'false_negatives'
                           ])

    return conf_df


def metrics_computation(df):
    metrics_list = []

    for i, row in df.iterrows():
        true_positives = row['true_positives']
        true_negatives = row['true_negatives']
        false_positives = row['false_positives']
        false_negatives = row['false_negatives']

        accuracy = (true_positives +
                    true_negatives) / (true_positives + true_negatives +
                                       false_positives + false_negatives)

        precision = true_positives / (true_positives + false_positives)

        recall = true_positives / (true_positives + false_negatives)

        specificity = true_negatives / (true_negatives + false_positives)

        balanced_accuracy = (recall + specificity) / 2

        mcc = ((true_positives * true_negatives) -
               (false_positives * false_negatives)) / (np.sqrt(
                   (true_positives + false_positives)) * np.sqrt(
                       (true_positives + false_negatives)) * np.sqrt(
                           (true_negatives + false_positives)) * np.sqrt(
                               (true_negatives + false_negatives)))

        f1_score = (2 * true_positives) / (2 * true_positives +
                                           false_positives + false_negatives)

        metrics_list.append([
            accuracy, precision, recall, specificity, balanced_accuracy, mcc,
            f1_score
        ])

    return metrics_list


def metrics_9(parsed_domtblouts,
              parsed_psiblast,
              gt,
              h_threshold=False,
              hi_threshold=False,
              p_threshold=False):
    index_metrics = list(parsed_domtblouts.keys()) + list(
        parsed_psiblast.keys())
    columns_metrics = [
        'accuracy', 'precision', 'recall', 'specificity', 'balanced_accuracy',
        'mcc', 'f1_score'
    ]

    metrics_list = []

    if h_threshold < 1 and h_threshold != False:
        h_threshold_name = str(h_threshold).split('.')[1]
    else:
        h_threshold_name = h_threshold

    if hi_threshold < 1 and hi_threshold != False:
        hi_threshold_name = str(hi_threshold).split('.')[1]
    else:
        hi_threshold_name = hi_threshold

    if p_threshold < 1 and p_threshold != False:
        p_threshold_name = str(p_threshold).split('.')[1]
    else:
        p_threshold_name = p_threshold

    for df in parsed_domtblouts.keys():
        if h_threshold:
            parsed_domtblouts[df] = parsed_domtblouts[df][
                parsed_domtblouts[df]['E-value'] < h_threshold]

        if hi_threshold:
            parsed_domtblouts[df] = parsed_domtblouts[df][
                parsed_domtblouts[df]['E-value'] < hi_threshold]

    for df in parsed_psiblast.keys():
        if p_threshold:
            parsed_psiblast[df] = parsed_psiblast[df][
                parsed_psiblast[df]['e_value'].astype(
                    dtype='float64') < p_threshold]

    conf_df = compute_con_matrix_9(parsed_domtblouts, parsed_psiblast, gt)

    # at the beginning: if file already there and it is not to be modified -> read it instead of computing it
    # if 'metrics_9_{0}_{1}_{2}.csv'.format(h_threshold_name, hi_threshold_name,
    # p_threshold_name) in os.listdir(cur + '\\data\\metrics'):
    if 'metrics_9_{0}_{1}_{2}.csv'.format(h_threshold_name, hi_threshold_name,
                                          p_threshold_name) in os.listdir(
                                              path.join(
                                                  'data', 'part_1',
                                                  'metrics')):
        # old_metrics_df = pd.read_csv(cur +
        # '\\data\\metrics\\metrics_9_{0}_{1}_{2}.csv'.format(h_threshold_name,
        # hi_threshold_name, p_threshold_name), index_col=0)
        old_metrics_df = pd.read_csv(path.join(
            'data', 'part_1', 'metrics',
            'metrics_9_{0}_{1}_{2}.csv'.format(h_threshold_name,
                                               hi_threshold_name,
                                               p_threshold_name)),
                                     index_col=0)

        if old_metrics_df.index.to_list() == index_metrics:
            return old_metrics_df, conf_df

        # computation metrics for each df
        else:
            metrics_list = metrics_computation(conf_df)
            metrics_df = pd.DataFrame(metrics_list,
                                      index=conf_df.index,
                                      columns=columns_metrics)

            metrics_df.to_csv(
                path.join(
                    'data', 'part_1', 'metrics',
                    'metrics_9_{0}_{1}_{2}.csv'.format(h_threshold_name,
                                                       hi_threshold_name,
                                                       p_threshold_name)))

            return metrics_df, conf_df

    else:
        metrics_list = metrics_computation(conf_df)
        metrics_df = pd.DataFrame(metrics_list,
                                  index=conf_df.index,
                                  columns=columns_metrics)

        metrics_df.to_csv(
            path.join(
                'data', 'part_1', 'metrics',
                'metrics_9_{0}_{1}_{2}.csv'.format(h_threshold_name,
                                                   hi_threshold_name,
                                                   p_threshold_name)))

        return metrics_df, conf_df


def plot_metrics(metrics_df):
    col = metrics_df.columns.to_list()
    x = np.arange(metrics_df.shape[0])
    c = np.random.rand(metrics_df.shape[0], 3)
    for i in col:
        fig, ax = plt.subplots(figsize=(20, 10))
        ax.bar(metrics_df.index, list(metrics_df[i]), color=c)
        ax.set_xticks(x)
        ax.set_xticklabels(metrics_df.index, rotation=65)
        if i == 'n_hits':
            ax.set_ylim((0, 100))
        else:
            ax.set_ylim((0, 1))
        plt.title(i)
        #fig.savefig(i+".png")


def plot_metrics_models(metrics_df, num=8):
    if num == 8:
        metrics_df = metrics_df.iloc[:, 1:]
    col = metrics_df.index.to_list()
    x = np.arange(metrics_df.shape[1])
    c = np.random.rand(metrics_df.shape[1], 3)
    for i in col:
        fig, ax = plt.subplots(figsize=(20, 10))
        ax.bar(metrics_df.columns, list(metrics_df.loc[i, :]), color=c)
        ax.set_xticks(x)
        ax.set_xticklabels(metrics_df.columns)
        ax.set_ylim((0, 1))
        plt.title(i)
        #fig.savefig(i+".png")


def plot_metrics_summary(metrics_df,
                         num=8,
                         threshold_hmms_e_value=None,
                         threshold_hmms_i_e_value=None,
                         threshold_pssm_e_value=None):
    if num == 8:
        metrics_df = metrics_df.iloc[:, 1:]
        
    x = np.arange(metrics_df.shape[1])

    list_bars = [
        metrics_df.iloc[i, :].to_list() for i in range(len(metrics_df))
    ]
    offsets = np.arange(len(metrics_df)) - 2
    color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'][:len(metrics_df)]

    width = 0.15

    plt.figure(figsize=(15, 10))

    ax = plt.subplot(111)
    for i, el in enumerate(list_bars):
        ax.bar(x + offsets[i] * width,
               el,
               width=width,
               color=color[i],
               align='center')

    ax.set_xticks(x)
    ax.set_xticklabels(metrics_df.columns.to_list())

    ax.legend([convert(el.split('.')[0]) for el in metrics_df.index.to_list()])

    plt.savefig(
        path.join(
            'data', 'part_1', 'metrics', 'best_models_{0}_{1}_{2}_{3}.png'.format(
                remove_dot(threshold_hmms_e_value),
                remove_dot(threshold_hmms_i_e_value),
                remove_dot(threshold_pssm_e_value), num)))

    plt.show()


def remove_dot(x):
    x = str(x)
    x = x.split('.')[1]

    return x


def convert(l):
    l = l.split('_')
    sigla = []
    if l[0] == "out":
        sigla.append("P")
    elif l[0] == "hmmsearch":
        sigla.append("H")
    
    sigla.append(l[2])
    if l[3] == "1":
        sigla.append("1")
    elif l[3] == "3":
        sigla.append("2")
    elif l[3] == "4":
        sigla.append("3")
    
    if "denoised1" in l:
        sigla.append("d")
    
    if sigla[0]=="P" and l[-1]=="1iterations":
        sigla.append("1_it")

    return "".join(sigla)