"""
Simple classes to build a graph that splits its nodes into disjoint components
"""
import itemsets

class Component(object):

    def __init__(self):
        super(Component, self).__init__()
        self.model = None

    def __str__(self):
        return 'model: %s' % str(self.model)


class Graph(object):

    def __init__(self):
        super(Graph, self).__init__()

        # Nodes in a component
        self.components = []


    def add_nodes(self, itemset, new_model):
        """
        Add a new node to the graph. Components will be re-computed, and
        the itemset added to one component. The model passed in
        will be setup with the itemsets in the new component.

        It is the responsibility of the caller to
        create and init the new_model.

        :param itemset:
        :param itemset: A model
        :return: The new model updated according to the graph component
        """

        # Find all intersecting and disjoint components
        intersecting_components = []
        disjoint_components = []
        while 0 < len(self.components):
            component = self.components.pop()
            if component.model.union_of_C & itemset != 0:
                intersecting_components.append(component)
            else:
                disjoint_components.append(component)

        # Create a new component as a merge of C of all intersecting
        # components with the new model
        new_component = Component()
        new_component.model = new_model
        new_component.model.C = [itemset]
        new_component.model.I = itemsets.singletons_of_itemset(itemset)
        for intersecting_component in intersecting_components:
                new_component.model.C += intersecting_component.model.C
                new_component.model.I =  new_component.model.I.union(intersecting_component.model.I)

        # Set the union of C on the new model
        new_component.model.union_of_C = itemsets.union_of_itemsets(new_component.model.C)

        self.components = disjoint_components + [new_component]

        return new_model, self.components


    def disjoint_itemsets(self):
        """
        Iterator for disjoint summaries C
        """
        for component in self.components:
            yield component.model.C


    def model_iterator(self):
        """
        Iterator for disjoint summaries C
        """
        for component in self.components:
            yield component.model

