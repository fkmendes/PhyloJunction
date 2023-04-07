import unittest
import statistics
import math

# pj imports
import phylojunction.utility.helper_functions as pjh
import phylojunction.calculation.math_utils as pjmath
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class TestBiSSEIncompleteSamplingTrees(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        total_n_states = 2

        rates_s0 = [
            sseobj.DiscreteStateDependentRate(name="lambda0", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
            sseobj.DiscreteStateDependentRate(name="mu0", val=0.5, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]),
            sseobj.DiscreteStateDependentRate(name="q01", val=0.25, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[0,1])
        ]

        rates_s1 = [
            sseobj.DiscreteStateDependentRate(name="lambda1", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[1,1,1]),
            sseobj.DiscreteStateDependentRate(name="mu1", val=0.5, event=sseobj.MacroevolEvent.EXTINCTION, states=[1]),
            sseobj.DiscreteStateDependentRate(name="q10", val=0.25, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[1,0])
        ]

        rates = rates_s0 + rates_s1
        
        sampling_probs = [
            sseobj.DiscreteStateDependentProbability(name="rho", val=0.25, state=0),
            sseobj.DiscreteStateDependentProbability(name="rho", val=0.75, state=1)
        ]

        matrix_atomic_rates = [ rates ]
        
        matrix_state_dep_probs = [ sampling_probs ]

        state_dep_rate_manager = sseobj.DiscreteStateDependentParameterManager(
            matrix_atomic_rates,
            total_n_states
        )

        state_dep_prob_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_state_dep_probs, total_n_states
            )

        event_handler = sseobj.MacroevolEventHandler(state_dep_rate_manager)
        
        prob_handler = \
            sseobj.DiscreteStateDependentProbabilityHandler(
                state_dep_prob_manager
            )

        cls.sse_stash_complete = sseobj.SSEStash(event_handler)
        cls.sse_stash_incomplete = sseobj.SSEStash(event_handler, prob_handler)


    def test_bisse_sampled_terminal_counts_incomplete_sampling(self):
        """
        """

        # setting up stopping conditions
        stop_condition = "age"
        stop_condition_value = [2.75] ## 2.5 time units
        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        # start_states_list = [0 for i in range(n_sim)]
        start_states_list = [0, 1] * 50
    
        # simulations
        sim_batches_complete = list()
        for i in range(n_batches):
            sse_sim = distsse.DnSSE(
                self.sse_stash_complete,
                stop_condition_value,
                n=n_sim,
                stop=stop_condition,
                origin=start_at_origin,
                start_states_list=start_states_list,
                epsilon=1e-8,
                runtime_limit=3000,
                debug=False)

            trs = sse_sim.generate()

            sim_batches_complete.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)


        batch_avgs_extant_sampled_terminals_complete_n0 = list()
        batch_avgs_extant_sampled_terminals_complete_n1 = list()
        for i, batch in enumerate(sim_batches_complete):
            approx_n0s_sampled_complete = \
                [ann_tr.alive_state_count_dict[0] * 0.25 for ann_tr in batch]
            approx_n1s_sampled_complete = \
                [ann_tr.alive_state_count_dict[1] * 0.75 for ann_tr in batch]
            
            mean_n0_sampled_complete = \
                statistics.mean(approx_n0s_sampled_complete)
            mean_n1_sampled_complete = \
                statistics.mean(approx_n1s_sampled_complete)

            batch_avgs_extant_sampled_terminals_complete_n0.\
                append(mean_n0_sampled_complete)
            batch_avgs_extant_sampled_terminals_complete_n1.\
                append(mean_n1_sampled_complete)

        expected_avg_sampled_incomplete_n0 = \
            statistics.mean(batch_avgs_extant_sampled_terminals_complete_n0)
        expected_avg_sampled_incomplete_n1 = \
            statistics.mean(batch_avgs_extant_sampled_terminals_complete_n1)


        sim_batches_incomplete = list()
        for i in range(n_batches):
            sse_sim = distsse.DnSSE(
                self.sse_stash_incomplete,
                stop_condition_value,
                n=n_sim,
                stop=stop_condition,
                origin=start_at_origin,
                start_states_list=start_states_list,
                epsilon=1e-8,
                runtime_limit=3000,
                debug=False)

            trs = sse_sim.generate()

            sim_batches_incomplete.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)

        batch_avgs_extant_sampled_terminals_incomplete_n0 = list()
        batch_avgs_extant_sampled_terminals_incomplete_n1 = list()
        n0_ci_overlap_count = 0
        n1_ci_overlap_count = 0
        for i, batch in enumerate(sim_batches_incomplete):
            approx_n0s_sampled_incomplete = \
                [ann_tr.alive_sampled_state_count_dict[0] for ann_tr in batch]
            approx_n1s_sampled_incomplete = \
                [ann_tr.alive_sampled_state_count_dict[1] for ann_tr in batch]
            
            mean_n0_sampled_incomplete = \
                statistics.mean(approx_n0s_sampled_incomplete)
            mean_n1_sampled_incomplete = \
                statistics.mean(approx_n1s_sampled_incomplete)

            stdevs_n0 = statistics.stdev(approx_n0s_sampled_incomplete)
            stdevs_n1 = statistics.stdev(approx_n1s_sampled_incomplete)
            
            sterr_n0 = stdevs_n0 / math.sqrt(n_sim)
            sterr_n1 = stdevs_n1 / math.sqrt(n_sim)

            diff_n0 = 1.96 * sterr_n0
            diff_n1 = 1.96 * sterr_n1

            ci_n0_low = mean_n0_sampled_incomplete - diff_n0
            ci_n0_hi = mean_n0_sampled_incomplete + diff_n0
            ci_n0 = "(" + str(ci_n0_low) + ", " + str(ci_n0_hi) + ")"
            ci_n1_low = mean_n1_sampled_incomplete - diff_n1
            ci_n1_hi = mean_n1_sampled_incomplete + diff_n1
            ci_n1 = "(" + str(ci_n1_low) + ", " + str(ci_n1_hi) + ")"
            # print("ci n0 = " + ci_n0 + " ci n1 = " + ci_n1)
 
            batch_avgs_extant_sampled_terminals_incomplete_n0.\
                append(mean_n0_sampled_incomplete)
            batch_avgs_extant_sampled_terminals_incomplete_n1.\
                append(mean_n1_sampled_incomplete)

            if (expected_avg_sampled_incomplete_n0 < \
                (mean_n0_sampled_incomplete + diff_n0)) \
                and (expected_avg_sampled_incomplete_n0 > \
                (mean_n0_sampled_incomplete - diff_n0)):
                n0_ci_overlap_count += 1

            if (expected_avg_sampled_incomplete_n1 < \
                (mean_n1_sampled_incomplete + diff_n1)) \
                and (expected_avg_sampled_incomplete_n1 > \
                (mean_n1_sampled_incomplete - diff_n1)):
                n1_ci_overlap_count += 1

        # print(statistics.mean(batch_avgs_extant_sampled_terminals_complete_n0))
        # print(statistics.mean(batch_avgs_extant_sampled_terminals_complete_n1))
        # print(statistics.mean(batch_avgs_extant_sampled_terminals_incomplete_n0))
        # print(statistics.mean(batch_avgs_extant_sampled_terminals_incomplete_n1))

        print("\n95% CI for sampled taxa in state 0 includes expectation " \
            + str(n0_ci_overlap_count) + " times.")
        print("\n95% CI for sampled taxa in state 1 includes expectation " \
            + str(n1_ci_overlap_count) + " times.")
        
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.05 * exp_count)

        self.assertAlmostEqual(
            n0_ci_overlap_count, exp_count,
            msg="Truth should be contained within 95%-CI " \
            + str(exp_count) + " (+/- " + str(a_delta) \
            + ") out of 100 times.", delta=a_delta)
        
        self.assertAlmostEqual(
            n1_ci_overlap_count, exp_count,
            msg="Truth should be contained within 95%-CI " \
            + str(exp_count) + " (+/- " + str(a_delta) \
            + ") out of 100 times.", delta=a_delta)


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
    # on the terminal, remember to add "src/phylojunction" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3.9 tests/distribution/test_dn_discrete_sse_bisse_incomplete_sampling.py
    # 
    # or
    #
    # $ python3.9 -m tests.distribution.test_dn_discrete_sse_bisse_incomplete_sampling
    #
    # or 
    #
    # $ python3.9 -m unittest tests.distribution.test_dn_discrete_sse_bisse_incomplete_sampling.TestBiSSEIncompleteSamplingTrees.test_bisse_sampled_terminal_counts_incomplete_sampling

    unittest.main()