"""
Class and specialized functions to build and work with a bacteria family tree.
"""

import numpy as np
from utils.dataset_helpers import abundance_matrix

NODE_GREEN_COLOR = '#43AC43'
NODE_RED_COLOR = '#C57B7B'
class Node(object):
    def __init__(self, tree):
        super(Node, self).__init__()
        self.parent = None
        self.children = []
        self.name = ""
        self.depth = 0
        self.clades = ""
        # Numpy array of the abundance colum
        # w.r.t. the clade at a given leaf
        self.abundances = None
        self.tree = tree
        self.tag = False
        self.meta_data = ''


    def to_xml(self):
        xml = ''
        size = 1
        subtree_abundance = sum(self.tree.abundance_column_in_subtree(self))

        if self.is_root():

            tag_name = 'Process'


            # Attributes on root node?
            xml = '<%s name="nullRoot" size="%s" total_abundance="%f">' % (tag_name, size, float(subtree_abundance))
            for child in self.children:
                xml += child.to_xml()
            xml += '</%s>' % tag_name

            return xml

        elif self.is_leaf():
            short_name = self.name.split('|')[-1]
            tag_name = 'depth_%d' % self.depth
            xml = '<%s name="%s" path="%s" size="%s" total_abundance="%f"/>' % (tag_name, short_name, self.clades, size, float(subtree_abundance))
            return xml

        else:
            short_name = self.name.split('|')[-1]
            tag_name = 'depth_%d' % self.depth
            xml = '<%s name="%s" path="%s" size="%s" total_abundance="%f">' % (tag_name, short_name, self.clades, size, float(subtree_abundance))
            for child in self.children:
                xml += child.to_xml()
            xml += '</%s>' % tag_name

            return xml


    def digraph_to_ancestor(self, ancestor, co_excluded=[]):

        color = NODE_GREEN_COLOR
        if self.name in co_excluded:
            color = NODE_RED_COLOR

        # Node names are not unique so use the label to set the node text
        unique_node_name = self.name.replace("|","").replace('/','')
        node_label = self.name.split('|')[-1].replace('/','')
        s = '\t%s [label="%s", color="%s", style=filled];\n' % (unique_node_name, node_label, color)

        node = self
        while node != ancestor:
            if not node.tag:
                parent_unique = node.parent.name.replace("|","").replace('/','')
                child_unique = node.name.replace("|","").replace('/','')
                s += '\t%s [label="%s"];\n' % (parent_unique, node.parent.name.split('|')[-1].replace('/',''))
                s += '\t%s [label="%s"];\n' % (child_unique, node.name.split('|')[-1].replace('/',''))
                s += '\t%s -> %s;\n' % (parent_unique, child_unique)
            node.tag = True
            node = node.parent

        return s


    def is_ancestor(self, node):
        """
        Returns Ture if the passed in node is an ancestor of
        this node
        :param node:
        :return:
        """

        parent = self.parent

        while not (parent is None):
            if parent == node:
                return True

        return False


    def is_leaf(self):
        """
        Note: only leaf nodes hold a abundance column
        :return: True if this is a leaf node
        """
        return (self.abundances is None) == False


    def is_root(self):
        return self.parent is None


    def is_in_lineage(self, other_node):
        """
        Returns true if the node and the passed in node
        share lineage
        :param other_node:
        :return:
        """
        return self.name in other_node.clades or other_node.name in self.clades


    def __str__(self):
        return self.name


class Tree(object):

    def __init__(self, ds):
        """
        Initialize a phylogenetic Tree for a given dataset.
        :param ds: Dataset
        :return: An initialized Tree
        """
        super(Tree, self).__init__()
        self.root = Node(self)
        self.root.name = 'Root'
        "Dataset related to the datasets"
        self.ds = ds
        # Submatrix of ds with bacteria abundances
        self.bacteria_abundances = None
        self.nodes = {}
        self.__build(self.ds)


    def __add_clades(self, node, clades, abundances, depth=0):
        """
        Recursively adds a list of clades to current datasets.
        :param node:
        :param clades:
        :return:
        """
        if depth < len(clades):

            # Unique name for the node includes the parent name
            # if not roow in the phylogenetic tree
            if depth == 0:
                name = clades[depth]
            else:
                name = clades[depth-1] + '|' + clades[depth]

            if not name in [n.name for n in node.children]:
                assert not (name in self.nodes), 'Attempting to add node twice'

                child = Node(self)
                child.name = name
                child.parent = node
                child.depth = depth + 1
                child.clades = '|'.join(clades[:depth+1])
                self.nodes[name] = child
                node.children.append(child)
                self.__add_clades(child, clades, abundances, depth + 1)
            else:
                child = [n for n in node.children if n.name == name][0]
                self.__add_clades(child, clades, abundances, depth + 1)

        # This is a leaf
        else:
            node.abundances = abundances


    def __build(self, ds):
        """
        Build the bacteria family datasets from a dataset.
        The data set is expected to be in the format of the datasets
        that can be created with bacteria parser, i.e. fir row from
        index [0][2:] hold bacterial clades
        row [1][2:] hold bacteria abundances
        :param ds: Daset with all know clades as columns.
        :return: Bacteria family datasets
        """
        self.ds = ds

        self.bacteria_abundances = abundance_matrix(np.array(ds))

        # Headers with clades
        bacteria_clades = ds[0][2:]

        # Insert all clades in the datasets
        for index, clade in enumerate(bacteria_clades):
            names = clade.split('|')
            if len(names) > 0:
                abundances = self.bacteria_abundances.T[index]
                # pass abundance as a column
                self.__add_clades(self.root, names, abundances.T)
            else:
                print 'bad column: ', names
        return self.root


    def abundance_column_in_subtree(self, node, column=None):
        """
        Returns the abundance column for a subtree.
        This is the sum of all abundance columns in the subtree

        :param node: Subtree root node
        :param column: Current subtree column
        :return: Column with abundances for the subtree
        """
        if column is None:
            # init 0 column with as many dimensions as in the abundance matrix
            column = np.array([0 for val in range(len(self.bacteria_abundances))])
            column = column.T

        # Leaf node, return the abundance column
        if not (node.abundances is None):
            return column + node.abundances

        # This is not a leaf, so add all children
        for child in node.children:
            column = self.abundance_column_in_subtree(child, column)

        return column

    def nodes_at_depth(self, node, depth, nodes=None):
        """
        Return nodes at at a given depth. Leafs at lower depth are ignored
        :param node: Current node
        :param depth: Depth of desired nodes
        :return: Nodes at a given depth
        """
        if nodes is None:
            nodes = []

        if depth == 0:
            return nodes + [node]
        else:
            for child in node.children:
                nodes = self.nodes_at_depth(child, depth - 1, nodes)
        return nodes


    def nodes_at_max_depth(self, max_depth, node, depth=0, nodes=None):
        """
        Return nodes at a max depth, or leafs if these are before max depth.
        :param max_depth: Max depth
        :param node: Starting node
        :param depth: Optional. Current depth
        :param nodes: Optional. Current nodes found
        :return:
        """
        if nodes is None:
            nodes = []

        if max_depth <= depth or node.is_leaf():
            return nodes + [node]
        else:
            for child in node.children:
                nodes = self.nodes_at_max_depth(max_depth, child, depth + 1, nodes)
        return nodes


    def nodes_for_clades_or_leaf(self, clade_names, node=None, nodes=None):
        """
        Returns all nodes in the tree that are either in the clade names
        or are leafs
        :param clade_names: A list of clade names
        :param node: Current node
        :param nodes: Current list of nodes
        :return:
        """
        if nodes is None:
            nodes = []
        if node is None:
            node = self.root

        if node.name in clade_names or node.is_leaf():
            return nodes + [node]
        else:
            for child in node.children:
                nodes = self.nodes_for_clades_or_leaf(clade_names, child, nodes)
        return nodes


    def dataset_at_depth(self, depth):
        """
        Return a dataset at at a given depth. Leafs at lower depths are ignored
        :param depth: Depth of desired nodes
        :return: Dataset with nodes at desired depth
        """
        nodes = self.nodes_at_depth(self.root, depth)
        return self.dataset_for_nodes(nodes)


    def dataset_for_clades_or_leaf(self, clade_names):
        """
        Returns a dataset that cuts the tree at the provided clades,
        and leaf nodes for subtrees where a clade name is not found
        :param clade_names: Clades that must be included
        :return:
        """
        nodes = self.nodes_for_clades_or_leaf(clade_names)
        return self.dataset_for_nodes(nodes)


    def dataset_at_leaves(self):
        """
        Build a dataset using only the leaf nodes
        :return:
        """
        nodes = [node for node in self.nodes if node.is_leaf()]
        return self.dataset_for_nodes(nodes)


    def dataset_for_all_nodes(self):
        """
        Builds a dataset for all nodes
        :return:
        """
        nodes = self.nodes.values()
        return self.dataset_for_nodes(nodes)


    def dataset_at_max_depth(self, max_depth):
        """
        Returns a data w.r.t a max depth. Subtrees at max deth will have the
        abundances merge, and the clade name is create for the reached node.
        :param max_depth: Max depth in bacteria datasets.
        :return: Dataset at max_depth in the datasets
        :return:
        """
        nodes = self.nodes_at_max_depth(max_depth, self.root)
        return self.dataset_for_nodes(nodes)


    def dataset_for_clades(self, clade_names):
        """
        Returns a dataset of merged abundances for the clade names provided
        :param clade_names: A list of clade names.
        :return:
        """
        nodes = []
        for name in clade_names:
            if name in self.nodes:
                nodes.append(self.nodes[name])
            else:
                print 'Warning, no known clade with name: ', name

        return self.dataset_for_nodes(nodes)


    def dataset_for_nodes(self, nodes):
        """
        Returns a subset of the dataset used to build the tree, where leaf abundances are merged.
        :return: Dataset
        """

        abundance_columns = None
        headers = []

        for node in nodes:

            if not node.is_root():

                # The new header name is the clade at the node
                # 'Bacteria|Firmicutes|' from 'Bacteria|Firmicutes|Bacilli' if we
                # cut at some node, not leaf
                headers.append(node.name)

                # Get the merged abundance columns for the subtree from the node
                column = self.abundance_column_in_subtree(node)

                # A bit of juggling here, the column is a row here
                column = np.array([column]).T

                if abundance_columns is None:
                    abundance_columns = column
                else:
                    abundance_columns = np.hstack((abundance_columns, column))

        # attach header row on top of abundances
        headers = np.array(headers)
        dataset = np.vstack((headers, abundance_columns))

        # Attach sample columns, on left side
        left_columns = np.array(self.ds)[:,0:2]
        dataset = np.hstack((left_columns, dataset))


        return dataset


    def abundance_for_clade(self, clade_name):
        """
        Get the abundance row at a given clade.
        Important that calde names are given with the parent name
        if any ex: Firmicutes|Erysipelotrichi
        :param clade_name:
        :return:
        """
        assert clade_name and not (clade_name == ''), 'No clade name provided'
        assert clade_name in self.nodes, 'No node for clade_name'

        node = self.node_for_clade_name(clade_name)

        return self.abundance_column_in_subtree(node)


    def count_nodes(self, node=None, depth=0, count=0, count_depth=None):
        """
        Returns the number of nodes in the datasets
        :param node: Current node
        :param depth: Current depth
        :param count: Current count
        :param count_depth: Optional max depth for the count
        :return: Number of nodes in the datasets
        """
        if node is None:
            node = self.root

        for child in node.children:
            count = self.count_nodes(child, depth+1, count, count_depth)
        if count_depth is None or depth <= count_depth:
            return 1 + count
        return count


    def count_leafs(self, node=None, count=0):
        """
        Returns the number of leafs
        :param node: Current node
        :param count: Current count
        :return: Number of leafs in the datasets
        """
        if node is None:
            node = self.root

        if node.is_leaf():
            return 1 + count
        else :
            for child in node.children:
                count = self.count_leafs(child, count)

        return count


    def node_for_clade_name(self, clade_name):
        assert clade_name and not (clade_name == ''), 'No clade name provided'

        return self.nodes[clade_name]


    def nodes_have_same_lineage(self, node_a, node_b):
        """
        Returns true if node a and node b are in the same lineage
        :param node_a:
        :param node_b:
        :return:
        """
        return node_a.is_in_lineage(node_b)


    def has_clade(self, clade_name):
        return clade_name in self.nodes


    def clear_node_tags(self):
        for node in self.nodes.values():
            node.tag = False
            node.meta_data = ''


    def dot_graph_for_clades(self, clades, co_excluded=[]):
        """
        Specialized method to produce a digraph
        file for the dot program in graphviz.

        By simply passing in a list of clades tuple
        a digraph tree is constructed with a least common
        ancestor.

        The list of clades passed in should be (pattern, co-exclusion)
        where co-exclusion should be None if it is not there.
        A co-exclusion node will be colored with red.

        :param clades:
        :return:
        """
        # given a list of clades, find common ancestors,
        # and construct subtree from involved nodes

        self.clear_node_tags()

        # Targeted nodes
        nodes = []
        for clade in clades:
            nodes.append(self.nodes[clade])


        # find minimum shared depth
        min_depth = nodes[0].depth
        for node in nodes:
            if node.depth < min_depth:
                min_depth = node.depth

        nodes_at_equal_depth = []
        # Move to share min depth
        for node in nodes:
            equal_depth_ancestor = node
            while equal_depth_ancestor.depth != min_depth:
                equal_depth_ancestor = equal_depth_ancestor.parent
            nodes_at_equal_depth.append(equal_depth_ancestor)

        # Initialize set of least common ancestors
        # with set of parent nodes
        least_common_ancestor = set(nodes_at_equal_depth)

        # Bubble up the tree to find least common ancestor
        while 1 < len(least_common_ancestor):
            in_common_ancestors = least_common_ancestor.copy()
            for node in in_common_ancestors:
                least_common_ancestor.remove(node)
                least_common_ancestor.add(node.parent)

        ancestor = least_common_ancestor.pop()
        digraph = 'digraph G { \n'
        for node in nodes:
            digraph += node.digraph_to_ancestor(ancestor, co_excluded)
        digraph += '\n }'

        self.clear_node_tags()

        return digraph


    def simple_node_name(self, clade):
        chunks = clade.split('|')
        name = ''
        if chunks[-1] == 'unclassified':
            name = '_'.join(chunks)
        else:
            name = chunks[-1]
        name = name.replace('/', '_or_')
        return name


    def dot_graph_for_summary(self, summary):
        """
        Specialized method to produce a graph
        file for the dot program in graphviz.

        Pass in a summary as jagged list of tuples.
        Tuples contain a list of the items in a pattern,
        and a list of co-excluded items
        Each pattern will be assigned a unique color,
        and all nodes the same color.

        :param Summary: Summary [[a,b,c], [c,d], ...]
        :return: A Graph string in the dot syntax
        """
        from itertools import combinations
        import random

        base_colors = ['#800000', '#ff0000', '#ffA500', '#ffff00', '#808000', '#800080', '#ff00ff', '#00ff00', '#008000', '#000080', '#0000ff', '#00ffff', '#008080', '#000000', '#c0c0c0']

        # Create set of all unique names
        # we have been provided.

        graph = 'graph G { \n'

        for index, (pattern, co_excluded) in enumerate(summary):

            graph += '\n'
            graph += '\n\t//Pattern: %i\n' % index

            # Check if this is a co-exclusion pattern
            is_cooccurrence_pattern = False
            for item in pattern:
                if not item in co_excluded:
                    is_cooccurrence_pattern = True
                    break

            for item in pattern:
                name = self.simple_node_name(item)

                node_color = 'white'
                if item in co_excluded and is_cooccurrence_pattern:
                    node_color = NODE_RED_COLOR
                    graph += '\t%s [color="%s", style=filled]\n' % (name, node_color)
                else:
                    graph += '\t%s\n' % name


            edge_color = ''
            if len(base_colors) > 0:
                edge_color = base_colors.pop()
            else:
                r = lambda: random.randint(0,255)
                edge_color = ('#%02X%02X%02X' % (r(),r(),r()))

            label_color = NODE_GREEN_COLOR
            if not is_cooccurrence_pattern:
                label_color = 'red'

            graph += '\tedge [color=\"%s\", penwidth=3, label="%d", fontcolor=\"%s\", fontname=\"Helvetica bold\", fontsize=24]\n' % (edge_color, index, label_color)
            graph += '\t' + self.simple_node_name(pattern[0])
            for item in pattern[1:]:
                graph += ' -- %s' % self.simple_node_name(item)


        #
        # Use this code to construct complete subgraphs
        #
        # for index, pattern in enumerate(clean_summary):
        #     graph += '\n'
        #     graph += '//Pattern: %i\n' % index
        #
        #     # Color the edge, use a random color if out of
        #     # base colors
        #     hex_color = ''
        #     if len(base_colors) > 0:
        #         hex_color = base_colors.pop()
        #     else:
        #         r = lambda: random.randint(0,255)
        #         hex_color = ('#%02X%02X%02X' % (r(),r(),r()))
        #     graph += 'edge [color=\"%s\", penwidth=3]\n' % hex_color
        #     graph += pattern[0]
        #     for item in pattern[1:]:
        #         graph += ' -- %s' % item
        #     # for comb in combinations(pattern, 2):
        #     #     graph += '%s -- %s\n' % (comb[0], comb[1])

        graph += '}'

        return graph


def test_tree():
    # Create a toy dataset, with the same .tab format as we use
    ds = [
        ['TID', 'STSite', 'Bacteria|Firmicutes|Bacilli', 'Bacteria|Firmicutes|Erysipelotrichi', 'Bacteria|Fusobacteria', 'unclassified'],
        [700107520, 'Tongue_dorsum',   0               ,               1                       ,            2           ,     7         ],
        [700107521, 'Tongue_dorsum',   10              ,              11                       ,           12           ,    19         ],
        [700107521, 'Tongue_dorsum',   20              ,              21                       ,           22           ,    99         ],

    ]

    tree = Tree()
    tree.build_bacteria_family_tree(ds)

    sub_ds = tree.dataset_at_max_depth(2)
    # Test be inspecting the result..
    print sub_ds

# test_tree()





