# Functions dedicated to the structural caraterization part of the project.
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import pandas as pd
import os
from os import path
import itertools
from scipy.cluster.hierarchy import linkage, dendrogram
sns.set_theme()
sns.set_style("whitegrid")


def generate_pdb_df(sifts_path, model_output_path):
    model_prots_df = pd.read_csv(model_output_path)
    model_prots = list(model_prots_df.ids.values)

    # Import pdb - uniprot relation file
    sifts_tsv = pd.read_csv(sifts_path, sep='\t', header=1)
    sifts_tsv.columns = list(map(lambda x: x.lower(),
                                 sifts_tsv.columns.values))
    pdb_df = sifts_tsv.loc[sifts_tsv.sp_primary.isin(model_prots), [
        'pdb', 'sp_primary', 'chain', 'pdb_beg', 'pdb_end', 'sp_beg', 'sp_end'
    ]].copy()
    pdb_df = pdb_df.reset_index()

    return pdb_df


def filter_pdb_db(pdb_db):
    denom = (pdb_db.sp_end - pdb_db.sp_beg + 1).astype(int)
    num = (pdb_db.pdb_end.astype(int) - pdb_db.pdb_beg.astype(int) + 1)
    coverage = num / denom
    pdb_db['coverage'] = coverage
    pdb_db = pdb_db[pdb_db.coverage > 0.8]

    return pdb_db


def parse_tmalign_out(filename):
    """
    Takes as input the path of the .out file provided by TM-Align.
    Returns both the RMSD and the TMSCORE (the first one provided in the .out file) """
    lines = []
    # temp_path =
    # ".\\data\\part_2\\original_datasets\\family_structures\\temp"
    temp_path = path.join('data', 'part_2', 'original_datasets',
                          'family_structures', 'temp')
    with open(temp_path + '\\' + filename, 'r') as f:
        for line in f:
            for line in itertools.islice(f, 15, 21):
                lines.append(line)
    RMSD = lines[0].split(',')[1].split('=')[1].strip()
    TMSCORE = lines[2].split('=')[1].split('(')[0].strip()

    return float(RMSD), float(TMSCORE)


def create_rmsd_matrix(best_model):
    """
    Takes all the files inside the temp folder, reads the rmsd values from them, 
    stores them into a matrix and saves it into a .csv.
    If the csv already exists, the function just reads it."""

    # model_path =
    # ".\\data\\part_2\\original_datasets\\family_structures\\pdbs_{}".format(best_model)
    model_path = path.join('data', 'part_2', 'original_datasets',
                           'family_structures', 'pdbs_{}'.format(best_model))

    # Check if the rmsds_....csv file is already present. If so, just read it.
    if 'rmsds_{}.csv'.format(best_model) in os.listdir(model_path):
        # rmsds_df = pd.read_csv(model_path + '\\' + 'rmsds_' + best_model +
        #                        '.csv',
        #                        index_col=0)
        rmsds_df = pd.read_csv(path.join(model_path,
                                         'rmsds_{}.csv'.format(best_model)),
                               index_col=0)
        return rmsds_df

    else:
        rmsd_dict = {}

        # dir_to_parse = ".\\data\\part_2\\original_datasets\\family_structures\\temp"
        dir_to_parse = path.join('data', 'part_2', 'original_datasets',
                                 'family_structures', 'temp')
        files = os.listdir(dir_to_parse)
        for filename in files:
            ls = filename.split("_")
            o1 = ls[0].split(".")[0]
            o2 = ls[1].split(".")[0]
            rmsd_dict.setdefault(o1, {})
            rmsd_dict[o1][o2], _ = parse_tmalign_out(filename)

        rmsd_df = pd.DataFrame.from_dict(rmsd_dict)
        rmsd_df.index = rmsd_df.index.map(lambda x: x[3:])
        rmsd_df.columns = rmsd_df.columns.map(lambda x: x[3:])
        # rmsd_df.to_csv(model_path + '\\' + 'rmsds_' + best_model + '.csv')
        rmsd_df.to_csv(path.join(model_path,
                                 'rmsds_{}.csv'.format(best_model)))

        return rmsd_df


def create_tmscores_matrix(best_model):
    """
    Takes all the files inside the temp folder, reads the tmscore values from them, 
   stores them into a matrix and saves it into a .csv.
    If the csv already exists, the function just reads it."""

    # model_path = ".\\data\\part_2\\original_datasets\\family_structures\\pdbs_{}".format(
    #     best_model)
    model_path = path.join('data', 'part_2', 'original_datasets',
                           'family_structures', 'pdbs_{}'.format(best_model))

    if 'tmscores_{}.csv'.format(best_model) in os.listdir(model_path):
        # tmscore_df = pd.read_csv(model_path + '\\' + 'tmscores_' + best_model +
        #                          '.csv',
        #                          index_col=0)
        tmscore_df = pd.read_csv(path.join(
            model_path, 'tmscores_{}.csv'.format(best_model)),
                                 index_col=0)
        return tmscore_df
    else:

        tmscore_dict = {}

        # dir_to_parse = ".\\data\\part_2\\original_datasets\\family_structures\\temp"
        dir_to_parse = path.join('data', 'part_2', 'original_datasets',
                                 'family_structures', 'temp')
        files = os.listdir(dir_to_parse)
        for filename in files:
            ls = filename.split("_")
            o1 = ls[0].split(".")[0]
            o2 = ls[1].split(".")[0]
            tmscore_dict.setdefault(o1, {})
            _, tmscore_dict[o1][o2] = parse_tmalign_out(filename)

        tmscore_df = pd.DataFrame.from_dict(tmscore_dict)
        tmscore_df.index = tmscore_df.index.map(lambda x: x[3:])
        tmscore_df.columns = tmscore_df.columns.map(lambda x: x[3:])
        # tmscore_df.to_csv(model_path + '\\' + 'tmscores_' + best_model +
        #                   '.csv')
        tmscore_df.to_csv(
            path.join(model_path, 'tmscores_{}.csv'.format(best_model)))
        return tmscore_df


def clear_temp_folder():
    # temp_dir = ".\\data\\part_2\\original_datasets\\family_structures\\temp"
    temp_dir = path.join('data', 'part_2', 'original_datasets',
                         'family_structures', 'temp')
    files = os.listdir(temp_dir)
    for filename in files:
        if filename.split('.')[-1] == 'out':
            # os.remove(temp_dir + '\\' + filename)
            os.remove(path.join(temp_dir, filename))
    return


def plot_heatmap(path_mat, filename, header=''):
    rmsdmatrix = pd.read_csv(path.join(path_mat, filename + '.csv'), index_col=0)
    col = list(rmsdmatrix.columns)
    col2 = {}
    for i in range(len(col)):
        col2.setdefault(col[i], col[i])
    rmsdmatrix = rmsdmatrix.rename(columns=col2)
    rmsdmatrix = rmsdmatrix.rename(index=col2)
    plt.figure(figsize=(20, 15))
    ax = sns.heatmap(rmsdmatrix,
                     cmap="Blues",
                     annot=False,
                     linewidths=.3,
                     linecolor="black",
                     vmin=0)
    plt.title(header)
    plt.savefig(path.join(path_mat, header + '.png'))
    plt.show()


def plot_dendogram(matrix, header, save_path):
    z = linkage(matrix, 'ward')
    fig = plt.figure(figsize=(20, 10))
    dn = dendrogram(z, labels=matrix.index, leaf_rotation=90)
    plt.title(header)
    plt.ylabel("Distance")
    plt.xlabel("PDBs")
    plt.savefig(save_path)