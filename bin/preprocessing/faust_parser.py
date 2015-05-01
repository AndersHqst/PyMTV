#
# Specialized scripts to parse the Faust results file: faust_results_dots.csv
#

import os
dir = os.path.dirname(__file__)

COLUMN_NODE_IDENTIFIER_1 = "Node identifier 1"
COLUMN_NODE_IDENTIFIER_2 = "Node  identifier 2" # Note double space is in the original file
COLUMN_CLADE_1 = "Clade 1"
COLUMN_CLADE_2 = "Clade 2"
COLUMN_BODY_SITE_1 = "Bodysite 1"
COLUMN_BODY_SITE_2 = "Bodysite 2"
COLUMN_NUMBER_OF_SUPPORTING_METHODS = "Number of dataset- or method-specific networks supporting relationship"
COLUMN_DIRECTION = "Direction"
COLUMN_RELATIONSHIP_BETWEEN_DIFFERENT_BODY_SITES = "Relationship between the same clade in different bodysites"
COLUMN_RELATIONSHIP_SAME_BODYSITE = "Relationship within the same bodysite"


class FaustResult(object):
    def __init__(self):
        super(FaustResult, self).__init__()

        # row for resutl in the faust result file
        self.id = 0

        # Body site and clade name separeted by |
        self.node_identifier_1 = ""
        self.node_identifier_2 = ""

        # Body site, ex Stool, tongue
        self.body_site_1 = ""
        self.body_site_2 = ""

        # Bacterial clade, ex Staphylococcaceae-Gemella
        self.clade_1 = ""
        self.clade_2 = ""

        # The number of methods from Fausts pipeline that supports
        # the relationship
        self.number_of_supporting_methods = 0

        # Direction of the relations ship,
        # Numbers > 0 signal positive relationship
        self.direction = 0

        ### Note the two following values seem fishy, they are rarely set in the original data???
        # True if the relationship is between bodysites
        self.relationship_between_different_body_site = False

        # True if the relationship is in the same bodysite
        self.relationship_same_body_site = False


    # def to_csv(self):
        # return self.clade_1 + ',' +

    def __str__(self):

        s = 'Bodysites: %s-%s, ' % (self.body_site_1, self.body_site_2)
        s += 'clades: %s %s, ' % (self.clade_1, self.clade_2)
        s += 'direction: %d, ' % self.direction
        s += 'supporting methods: %d' % self.number_of_supporting_methods

        return s



def clean_faust_result(path='faust_results.csv', out='faust_results_dots.csv'):
    """
    Helper script to convert commas in strings
     to dots
    :return:
    """
    import csv
    fr = open(path, 'rb')
    csv_reader = csv.reader(fr)
    f_out = open(out, 'wb')
    csv_writer = csv.writer(f_out)
    for row in csv_reader:
        clean_row = [val.replace(',','.') for val in row]
        csv_writer.writerow(clean_row)
    f_out.close()
    fr.close()


def faust_results(csv_file):
    """
    Returns a list of parsed FaustResult objects
    :param csv_file:
    :return:
    """
    import csv

    # List of FaustResultObjects
    results = []

    # Parse the file
    fd = open(csv_file, 'rb')
    csv_reader = csv.reader(fd)

    # Get headers
    headers = csv_reader.next()

    for index, row in enumerate(csv_reader):
        try:
            faust_result = FaustResult()
            faust_result.id = index
            faust_result.node_identifier_1 = row[headers.index(COLUMN_NODE_IDENTIFIER_1)]
            faust_result.node_identifier_2 = row[headers.index(COLUMN_NODE_IDENTIFIER_2)]

            faust_result.body_site_1 = row[headers.index(COLUMN_BODY_SITE_1)]
            faust_result.body_site_2 = row[headers.index(COLUMN_BODY_SITE_2)]

            # Replace '-' with a '|' in the clade names, this is the format used elsewhre
            faust_result.clade_1 = row[headers.index(COLUMN_CLADE_1)].replace('-', '|')
            faust_result.clade_2 = row[headers.index(COLUMN_CLADE_2)].replace('-', '|')

            try:
                faust_result.number_of_supporting_methods = int(row[headers.index(COLUMN_NUMBER_OF_SUPPORTING_METHODS)])
            except ValueError, e:
                pass


            try:
                faust_result.direction = int(row[headers.index(COLUMN_DIRECTION)])
            except ValueError, e:
                pass

            try:
                faust_result.relationship_between_different_body_site = int(row[headers.index(COLUMN_RELATIONSHIP_BETWEEN_DIFFERENT_BODY_SITES)]) != 0
            except ValueError, e:
                pass

            try:
                faust_result.relationship_same_body_site = int(row[headers.index(COLUMN_RELATIONSHIP_BETWEEN_DIFFERENT_BODY_SITES)]) != 0
            except ValueError, e:
                pass

            results.append(faust_result)
        except Exception, e:
            print e
            print 'Error parsing result line: ', row
            print 'Partially parsed FaustResult: ', faust_result

    fd.close()

    return results


def filtered_results(bodysite=None, same_bodysite_only=False):

    file_name = os.path.join(dir, '../../data/faust_results_dots.csv')

    filtered_results = faust_results(file_name)

    if not (bodysite is None):
        filtered_results = filter(lambda faust_result: faust_result.body_site_1 == bodysite or faust_result.body_site_2 == bodysite, filtered_results)

    if same_bodysite_only:
        filtered_results = filter(lambda faust_result: faust_result.body_site_1 == faust_result.body_site_2, filtered_results)

    return filtered_results


def results(bodysite='Stool'):
    """
    Helper method to get inter-site Stool results
    :return:
    """
    return filtered_results(bodysite=bodysite, same_bodysite_only=True)




