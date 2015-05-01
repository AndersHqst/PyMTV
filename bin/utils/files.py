
def parse_dat_file(path):
    """
    Parses a .dat file to a dataset of binary vectors
    represented by integers.
    The expected input file list binary position separated
    by spaces
    ex: 0 4 42 189\n
    :param path:
    :return:
    """
    D = []
    for line in open(path):
        clean_line = line.replace('\n', '')
        cmps = clean_line.split(' ')
        if '' in cmps:
            cmps.remove('')
        if 0 < len(cmps):
            bin_positions = [int(cmp) for cmp in cmps]
            val = 0
            for bp in bin_positions:
                val = val | 2 ** bp
            D.append(val)
    return D

def parse_header_file(path):
    """
    Parses a header file where headers can correspond to
    attribute indeces. Header files are expected to be as the
    same format as .dat files with header names seperated eith a space
    :param path:
    :return:
    """
    with open(path) as fd:

        lines = fd.readlines()

        if len(lines) != 1:
            print 'Header files should only contains 1 line. Ignoring header file: ', path
            return None

        clean_line = lines[0].replace('\n', '')

        return clean_line.split(' ')


def write_dat_file(file_name, data):
    """
    Write itemsets to .dat file format. Each item is represented by
    the positions of ones in binary
    :param file_name:
    :param data:
    :return:
    """
    import itemsets
    with open(file_name, 'wb') as fd:
        for X in data:
            line = ' '.join([str(x) for x in itemsets.to_index_list(X)]) + '\n'
            fd.write(line)


def write_tab_file(file_name, dataset):
    """
    Writes a dataset to tab seperated .tab file
    :param file_name:
    :param dataset:
    :return:
    """
    import csv
    with open(file_name, 'wb') as fd:
        csv_writer = csv.writer(fd, delimiter='\t')
        csv_writer.writerows(dataset)
