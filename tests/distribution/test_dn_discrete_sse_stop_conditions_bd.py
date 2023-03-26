import unittest

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class TestSSEStopConditionsBD(unittest.TestCase):

    @classmethod
    def setUp(cls):
        # not state-dependent (just state 0, and no transition)
        rates_t0_s0 = [
            sseobj.MacroevolStateDependentRateParameter(name="lambda", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
            sseobj.MacroevolStateDependentRateParameter(name="mu", val=0.9, event=sseobj.MacroevolEvent.EXTINCTION, states=[0])
        ]
    
        # original implementation
        matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice
        
        fig_rates_manager = sseobj.FIGRatesManager(matrix_atomic_rate_params, 1)

        cls.event_handler = sseobj.MacroevolEventHandler(fig_rates_manager)


    def test_tree_size_stop_condition_origin(self):
        """
        Test if birth-death trees have correct number of tips, starting from origin
        """
        # setting up stopping conditions
        stop_condition = "size"
        stop_condition_value = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        start_at_origin = True

        # simulation initialization
        n_sim = 10
        start_states_list = [0 for i in range(n_sim)]
        # seeds_list = [i+1 for i in range(n_sim)]
    
        # simulations
        dnsse = distsse.DnSSE(self.event_handler, stop_condition_value, n=n_sim, stop=stop_condition, origin=start_at_origin, start_states_list=start_states_list, condition_on_survival=True, epsilon=1e-12)
        trs = dnsse.generate()

        tr_sizes = [ann_tr.n_extant_terminal_nodes for ann_tr in trs]

        # for ann_tr in trs:
        #     print(ann_tr.tree.as_string(schema="newick", suppress_internal_taxon_labels=True))

        self.assertEqual(tr_sizes, stop_condition_value)


    def test_tree_size_stop_condition_root(self):
        """
        Test if birth-death trees have correct number of tips, starting from root
        """
        # setting up stopping conditions
        stop_condition = "size"
        stop_condition_value = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        start_at_origin = False

        # simulation initialization
        n_sim = 10
        start_states_list = [0 for i in range(n_sim)]
        # seeds_list = [i+1 for i in range(n_sim)]
    
        # simulations
        dnsse = distsse.DnSSE(self.event_handler, stop_condition_value, n=n_sim, stop=stop_condition, origin=start_at_origin, start_states_list=start_states_list, condition_on_survival=True, epsilon=1e-12)
        trs = dnsse.generate()

        tr_sizes = [ann_tr.n_extant_terminal_nodes for ann_tr in trs]

        # for ann_tr in trs:
        #     print(ann_tr.tree.as_string(schema="newick", suppress_internal_taxon_labels=True))

        self.assertEqual(tr_sizes, stop_condition_value)


    def test_tree_height_stop_condition_origin(self):
        """
        Test if birth-death trees have correct height, starting from origin
        """
        # setting up stopping conditions
        stop_condition = "age"
        stop_condition_value = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        start_at_origin = True

        # simulation initialization
        n_sim = 10
        start_states_list = [0 for i in range(n_sim)]
        # seeds_list = [i+1 for i in range(n_sim)]
    
        # simulations
        dnsse = distsse.DnSSE(self.event_handler, stop_condition_value, n=n_sim, stop=stop_condition, origin=start_at_origin, start_states_list=start_states_list, condition_on_survival=True, epsilon=1e-12)
        trs = dnsse.generate()

        tr_sizes = [ann_tr.origin_age for ann_tr in trs]

        # for ann_tr in trs:
        #     print(ann_tr.tree.as_string(schema="newick", suppress_internal_taxon_labels=True))

        for idx, tr_size in enumerate(tr_sizes):
            self.assertAlmostEqual(tr_size, stop_condition_value[idx], delta=1e-12)


    def test_tree_height_stop_condition_root(self):
        """
        Test if birth-death trees have correct height, starting from root
        """
        # setting up stopping conditions
        stop_condition = "age"
        stop_condition_value = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        start_at_origin = False

        # simulation initialization
        n_sim = 10
        start_states_list = [0 for i in range(n_sim)]
        # seeds_list = [i+1 for i in range(n_sim)]
    
        # simulations
        dnsse = distsse.DnSSE(self.event_handler, stop_condition_value, n=n_sim, stop=stop_condition, origin=start_at_origin, start_states_list=start_states_list, condition_on_survival=True, epsilon=1e-12)
        trs = dnsse.generate()

        tr_sizes = [ann_tr.root_age for ann_tr in trs]

        # for ann_tr in trs:
        #     print(ann_tr.tree.as_string(schema="newick", suppress_internal_taxon_labels=True))

        for idx, tr_size in enumerate(tr_sizes):
            self.assertAlmostEqual(tr_size, stop_condition_value[idx], delta=1e-12)


    def test_invalid_stop_condition_value(self):
        """
        Test if invalid stop condition values are correctly caught and
        the appropriate exceptions are raise:
            (1) Negative number of terminal nodes,
            (2) Number of terminal nodes being a float that cannot be
                converted to integer,
            (3) Negative tree height
        """
        # simulation initialization
        n_sim = 2
        start_states_list = [0 for i in range(n_sim)]

        stop_condition = "size"
        stop_condition_value = [1, -1]
        start_at_origin = True

        # (1) Negative number of terminal nodes
        with self.assertRaises(ec.DnInitMisspec) as exc:
            distsse.DnSSE(
                self.event_handler,
                stop_condition_value,
                n=n_sim,
                stop=stop_condition,
                origin=start_at_origin,
                start_states_list=start_states_list,
                condition_on_survival=True,
                epsilon=1e-12)       
        
        self.assertEqual(str(exc.exception),
            "ERROR: Distribution DnSSE was not properly initialized. " + \
            "Stop condition value (number of terminal nodes) cannot be " + \
            "negative. Exiting...")

        # (2) Number of terminal nodes being a float that cannot be converted to integer
        stop_condition_value = [0.5, 2]
        with self.assertRaises(ec.DnInitMisspec) as exc:
            distsse.DnSSE(
                self.event_handler,
                stop_condition_value,
                n=n_sim,
                stop=stop_condition,
                origin=start_at_origin,
                start_states_list=start_states_list,
                condition_on_survival=True,
                epsilon=1e-12)       
        
        self.assertEqual(str(exc.exception),
            "ERROR: Distribution DnSSE was not properly initialized. " + \
            "Stop condition value (number of terminal nodes) could not be " + \
            "converted to integer. Exiting...")

        # (3) Negative tree height
        stop_condition = "age"
        stop_condition_value = [0.0, -0.1]

        with self.assertRaises(ec.DnInitMisspec) as exc:
            distsse.DnSSE(
                self.event_handler,
                stop_condition_value,
                n=n_sim,
                stop=stop_condition,
                origin=start_at_origin, 
                start_states_list=start_states_list,
                condition_on_survival=True,
                epsilon=1e-12)       
        
        self.assertEqual(str(exc.exception),
            "ERROR: Distribution DnSSE was not properly initialized. " + \
            "Stop condition value (tree height) cannot be negative. " + \
            "Exiting...")


if __name__ == "__main__":
    # Assuming you opened the PhyloJunction/ (repo root) folder
    # on VSCode and that you want to run this as a standalone script,
    # i.e., "Run Without Debugging", you will need to configure your
    # launch.json to have:
    # 
    # "env": {"PYTHONPATH": "${workspaceRoot}/src/phylojunction/"}
    # 
    # and your settings.json to have:
    #   
    # "python.analysis.extraPaths": [ "${workspaceFolder}/src/phylojunction/" ]
    # 
    # If you want to run this as a standalone from PhyloJunction/
    # on the terminal, remember to add "src/phylojunction" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3 tests/distribution/test_dn_discrete_sse_stop_conditions_bd.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_stop_conditions_bd
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_stop_conditions_bd.TestSSEStopConditionsBD.test_tree_size_stop_condition_origin

    unittest.main()