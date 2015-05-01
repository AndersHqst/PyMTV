from unittest import TestCase

from graph import Graph
import test_helper
from utils import charitems


class TestGraph(TestCase):

    def test_add_nodes(self):

        ab = charitems.to_binary('ab')
        cd = charitems.to_binary('cd')
        bc = charitems.to_binary('bc')

        G = Graph()
        assert len(G.components) == 0

        # Add first model with nodes to the graph
        test_model1 = test_helper.init_simple_model()
        G.add_nodes(ab, test_model1)
        assert len(G.components) == 1
        assert G.components[0].model == test_model1
        assert len(G.components[0].model.C) == 1
        assert ab in G.components[0].model.C

        # Added disjoint model
        test_model2 = test_helper.init_simple_model()
        G.add_nodes(cd, test_model2)
        assert len(G.components) == 2
        assert test_model2 in [model for model in G.model_iterator()]
        assert len(G.components[0].model.C) == 1
        assert len(G.components[1].model.C) == 1
        for model in G.model_iterator():
            assert len(model.C) == 1
            assert ab in model.C or cd in model.C
            if ab in model.C:
                assert not (cd in model.C)
            elif cd in model.C:
                assert not (ab in model.C)
            else:
                assert False


        # Add model that will join all existing components
        test_model3 = test_helper.init_simple_model()
        G.add_nodes(bc, test_model3)
        assert len(G.components) == 1
        models = [model for model in G.model_iterator()]
        assert len(models) == 1
        assert test_model3 in models
        assert len(G.components[0].model.C) == 3
        assert ab in test_model3.C
        assert cd in test_model3.C
        assert bc in test_model3.C



