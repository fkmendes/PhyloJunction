import unittest
import math
import statistics

# pj imports
import phylojunction.utility.helper_functions as pjh
import phylojunction.calculation.math_utils as pjmath
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class TestBDTrees(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # not state-dependent (just state 0, and no transition)
        total_n_states = 1

        # not state-dependent (just state 0, and no transition)
        rates_t0_s0 = [ sseobj.AtomicSSERateParameter(name="lambda", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
                        sseobj.AtomicSSERateParameter(name="mu", val=0.8, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]) ]

        # original implementation
        # matrix_atomic_rate_params = [ [rates_t0_s0] ] # 1D: time slices, 2D: states, 3D: parameters of state, several parameters -> matrix
        matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

        fig_rates_manager = sseobj.FIGRatesManager(matrix_atomic_rate_params, total_n_states)

        cls.event_handler = sseobj.MacroEvolEventHandler(fig_rates_manager)

    def test_expected_size_bd(self):
        """
        Test if birth-death trees have expected number of extant observable nodes
        """
        # setting up stopping conditions
        stop_condition = "age"
        stop_condition_value = [5.0] ## 5 time units
        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]
        # seeds_list = [i+1 for i in range(n_sim)]

        # simulations
        sim_batches = list()
        for i in range(n_batches):
            # print("Doing batch " + str(n_batches - i))
            sse_sim = distsse.DnSSE(self.event_handler, stop_condition_value, n=n_sim, stop=stop_condition, origin=start_at_origin,
                    start_states_list=start_states_list, epsilon=1e-12, runtime_limit=3000,
                    debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)

        # parsing results
        exp_size = pjmath.exp_extant_count_bd(1.0, 0.8, 5.0)

        # origin_heights_when_gt_0 = set()
        batch_cis = list()
        in_ci_count = 0
        for i, batch in enumerate(sim_batches):
            n_extant_obs_nodes = [ann_tr.n_extant_obs_nodes for ann_tr in batch]
            stdevs = statistics.stdev(n_extant_obs_nodes)
            sterr = stdevs / math.sqrt(n_sim)
            diff = 1.96 * sterr
            avg_size = statistics.mean(n_extant_obs_nodes)
            batch_cis = (avg_size - diff, avg_size + diff)

            if exp_size >= batch_cis[0] and exp_size <= batch_cis[1]:
                in_ci_count += 1

        print("\n95% CI includes expectation " + str(in_ci_count) + " times.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(in_ci_count, exp_count,
                                msg="Truth should be contained within 95%-CI " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)

if __name__ == '__main__':
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
    # on the terminal, remember to add "src/" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3 tests/distribution/test_dn_discrete_sse_bd.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_bd
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_bd.TestBDTrees.test_expected_size_bd

    unittest.main()