"""
Parser scripts to parse hmp.lq.phylotype.filter.tab containing
a sample per colum that can come from any of the sample centers.
Rows are bacteria down to various depths of their family tree,
with their abundance counts for each sample
"""

import os
dir = os.path.dirname(__file__)
import csv
from preprocessors import *

COLUMN_TID = 0
COLUMN_BODY_SITE = 1


def parse_dataset(csv_file):
    """
    Returns a transposed numpy array of the dataset with abundance levels transformed to int
    These are the samples that can be found in hmp.lq.phylotype.filter.tab

    One sample is arranged per column, and a row represents a bacteria.
    Cells contain abundance counts.
    :param csv_file:
    :return:
    """

    fd = open(csv_file, 'rb')
    csv_reader = csv.reader(fd, delimiter='\t')

    # Pop the TID and STSite rows
    matrix = np.array(csv_reader.next())
    matrix = np.vstack([matrix, csv_reader.next()])

    # Read in entire matrix and transform abundance values
    # these are in the sub block from index row:2 col:1
    for row in csv_reader:
        bacteria = row[0]
        abundances = []
        for x in row[1:]:
            val = 0
            if x != '':
                val = int(x)
            abundances.append(val)
        new_row = np.array([bacteria] + abundances)
        matrix = np.vstack([matrix, new_row])

    fd.close()

    # Return a transposed matrix as
    # we work with samples as rows
    return matrix.T

def dataset_at_bodyset(csv_file, bodysite):
    dataset = parse_dataset(csv_file)
    # Return row 1 with headers + sample at desired bodysite
    return np.vstack([dataset[:1], dataset[dataset[:,COLUMN_BODY_SITE] == bodysite]])


def save_sample(bodysite='Stool'):
    """
    Helper function to save a stool sub-sample file to the data folder
    :return:
    """
    faust_file = os.path.join(dir, '../../data/hmp.lq.phylotype.filter.tab')
    dataset = dataset_at_bodyset(faust_file, bodysite)

    file_name = os.path.join(dir, '../../data/%s.tab' % bodysite)

    fd = open(file_name, 'wb')

    csv_writer = csv.writer(fd, delimiter='\t')
    csv_writer.writerows(dataset)

    fd.close()


def get_dataset(bodysite='Stool'):
    """
    Read the daset from file
    :param file_name:
    :return:
    """
    file_name = os.path.join(dir, '../../data/%s.tab' % bodysite)
    fd = open(file_name, 'rb')
    csv_reader = csv.reader(fd, delimiter='\t')
    # header rows
    matrix = np.array(csv_reader.next())
    matrix = np.vstack([matrix, csv_reader.next()])
    for index, row in enumerate(csv_reader):
        id_cols = row[:2]
        abundances = [int(x) for x in row[2:]]
        new_row = np.array(id_cols + abundances)
        matrix = np.vstack([matrix, new_row])
    return data_cleaning(matrix)
