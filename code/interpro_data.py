import pandas as pd
import urllib
import json
import time
import os
from os import path

cur = os.getcwd()

metadata_columns = [
    'key', 'accession', 'name', 'source_database', 'length',
    'source_organism_taxID', 'source_organism_scientificName',
    'source_organism_fullName'
]
entries_columns = [
    'key', 'accession', 'entry_protein_locations_fragments_start',
    'entry_protein_locations_fragments_end',
    'entry_protein_locations_fragments_dc-status',
    'entry_protein_locations_model', 'entry_protein_locations_score',
    'protein_length', 'source_database', 'entry_type', 'entry_integrated'
]


def json_to_df(df_metadata, df_entries, json, key):
    i = key
    for el in json['results']:
        meta = []

        for key in el['metadata'].keys():
            if key != 'source_organism':
                meta.append(el['metadata'][key])

            else:
                for key_2 in el['metadata'][key].keys():
                    meta.append(el['metadata'][key][key_2])

        entry = []

        for key in el['entries'][0]:
            if key != 'entry_protein_locations':
                entry.append(el['entries'][0][key])

            else:
                for key_2 in el['entries'][0][key][0].keys():
                    if key_2 != 'fragments':
                        entry.append(el['entries'][0][key][0][key_2])

                    else:
                        for key_3 in el['entries'][0][key][0][key_2][0].keys():
                            entry.append(
                                el['entries'][0][key][0][key_2][0][key_3])

        df_metadata = df_metadata.append(dict(zip(metadata_columns,
                                                  [i] + meta)),
                                         ignore_index=True)
        df_entries = df_entries.append(dict(zip(entries_columns, [i] + entry)),
                                       ignore_index=True)

        i += 1

    return df_metadata, df_entries


def get_next_json(url):
    response = urllib.request.urlopen(url)
    data = json.loads(response.read())
    return data


def get_df(url, team):
    data = get_next_json(url)

    metadata = pd.DataFrame(columns=metadata_columns)
    entries = pd.DataFrame(columns=entries_columns)

    while data['next'] is not None:
        try:
            if len(metadata) > 0:
                i = metadata.key.to_list()[-1] + 1
            else:
                i = 0

            metadata, entries = json_to_df(metadata, entries, data, i)

            data = get_next_json(data['next'])

            if data['next'] is None:
                if len(metadata) > 0:
                    i = metadata.key.to_list()[-1] + 1
                else:
                    i = 0

                metadata, entries = json_to_df(metadata, entries, data, i)

        except:
            checkpoint(metadata, entries, flag, team)
            time.sleep(900)

            if len(metadata) > 0:
                i = metadata.key.to_list()[-1] + 1
            else:
                i = 0

            remaining = qty - i
            print('Remaining:', remaining, 'proteins')

            metadata, entries = json_to_df(metadata, entries, data, i)

            data = get_next_json(data['next'])

            if data['next'] is None:
                if len(metadata) > 0:
                    i = metadata.key.to_list()[-1] + 1
                else:
                    i = 0

                metadata, entries = json_to_df(metadata, entries, data, i)

    checkpoint(metadata, entries, team)

    return metadata, entries


def get_data(url, team):
    metadata, entries = get_df(url, team)
    #print('Got data')

    gt = metadata.copy()
    gt = gt[['accession']]
    gt['start'] = entries['entry_protein_locations_fragments_start']
    gt['end'] = entries['entry_protein_locations_fragments_end']
    gt['length'] = entries['protein_length']
    # gt.to_csv(cur+'\\data_team_'+str(team)+'\\ground_truth\\ground_truth.csv')
    gt.to_csv(path.join('data', 'part_1', 'ground_truth', 'ground_truth.csv'))

    return metadata, entries, gt


def checkpoint(df_metadata, df_entries, team):
    # df_metadata.to_csv(cur + '\\data_team_' + str(team) +
    #                    '\\ground_truth\\metadata\\metadata.csv')
    df_metadata.to_csv(
        path.join('data', 'part_1', 'ground_truth', 'metadata',
                  'metadata.csv'))
    # df_entries.to_csv(cur + '\\data_team_' + str(team) +
    #                   '\\ground_truth\\entries\\entries.csv')
    df_entries.to_csv(
        path.join('data', 'part_1', 'ground_truth', 'entries', 'entries.csv'))
