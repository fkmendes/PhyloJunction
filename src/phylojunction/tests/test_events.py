import sys
sys.path.extend(["../", "../phylojunction"]) # necessary to run it as standalone on command line (from phylojunction/ or phylojunction/calculation/)
import unittest

# pj imports
import calculation.discrete_sse as sseobj

class TestEvent(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # birth-death
        total_n_states = 1
        cls.bd_rates_t0_s0 = [ sseobj.AtomicSSERateParameter(name="lambda", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
                            sseobj.AtomicSSERateParameter(name="mu", val=0.5, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]) ]
        
        bd_matrix_atomic_rate_params = [ cls.bd_rates_t0_s0 ] # 1D: time slices, 2D: states, 3D: parameters of state, several parameters -> matrix
        bd_fig_rates_manager = sseobj.FIGRatesManager(bd_matrix_atomic_rate_params, total_n_states)
        cls.bd_event_handler = sseobj.MacroEvolEventHandler(bd_fig_rates_manager)
        cls.bd_state_representation_dict = { 0: ["nd3", "nd4", "nd5", "nd6"] }


        # BiSSE
        total_n_states = 2
        cls.bisse_rates_t0_s0 = [ sseobj.AtomicSSERateParameter(name="lambda0", val=0.5, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
                                  sseobj.AtomicSSERateParameter(name="mu0", val=0.25, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]),
                                  sseobj.AtomicSSERateParameter(name="q01", val=0.5, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[0,1]) ]
        cls.bisse_rates_t0_s1 = [ sseobj.AtomicSSERateParameter(name="lambda1", val=1.5, event=sseobj.MacroevolEvent.W_SPECIATION, states=[1,1,1]),
                                  sseobj.AtomicSSERateParameter(name="mu1", val=0.25, event=sseobj.MacroevolEvent.EXTINCTION, states=[1]),
                                  sseobj.AtomicSSERateParameter(name="q10", val=0.5, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[1,0]) ]

        bisse_rates_t0 = cls.bisse_rates_t0_s0 + cls.bisse_rates_t0_s1
        
        bisse_matrix_atomic_rate_params = [ bisse_rates_t0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice
        bisse_fig_rates_manager = sseobj.FIGRatesManager(bisse_matrix_atomic_rate_params, total_n_states)
        cls.bisse_event_handler = sseobj.MacroEvolEventHandler(bisse_fig_rates_manager)
        cls.bisse_state_representation_dict = { 0: ["nd3"], 1: ["nd4", "nd5", "nd6"] }

    def test_total_rate_single_epoch_bd(self):
        """
        Test if birth-death total rates are correctly computed for single epoch
        """
        a_time = 0.0

        total_rate_t0, state_total_rates = self.bd_event_handler.total_rate(a_time, self.bd_state_representation_dict)
        total_rate_t100, state_total_rates = self.bd_event_handler.total_rate(a_time, self.bd_state_representation_dict) # should not matter what time

        # 4 * 1.5 = 6.0
        self.assertEqual(total_rate_t0, 6.0, "Total rate at time 0.0 should be 6.0.")
        self.assertEqual(total_rate_t100, 6.0, "Total rate at time 100.0 should be 6.0.")


    def test_total_rate_single_epoch_bisse(self):
        """
        Test if global and state-conditioned total rates are correctly calculated
        """
        a_time = 0.0

        # overall rate for drawing time to next event
        total_rate_t0, state_total_rates = self.bisse_event_handler.total_rate(a_time, self.bisse_state_representation_dict)
        total_rate_t100, state_total_rates = self.bisse_event_handler.total_rate(a_time, self.bisse_state_representation_dict) # should not matter what time

        # 1 * (0.5 + 0.25 + 0.5) + 3 * (1.5 + 0.25 + 0.5) = 8.0
        self.assertEqual(total_rate_t0, 8.0, "Total rate at time 0.0 should be 8.0.")
        self.assertEqual(total_rate_t100, 8.0, "Total rate at time 100.0 should be 8.0.")

        # 0.5 + 0.25 + 0.5 = 1.25
        self.assertEqual(state_total_rates[0], 1.25, "Total rate of state 0 at time 0.0 should be 1.25.")
        self.assertEqual(state_total_rates[0], 1.25, "Total rate of state 0 at time 100.0 should be 1.25.")

        # 1.5 + 0.25 + 0.5 = 2.25
        self.assertEqual(state_total_rates[1], 2.25, "Total rate of state 1 at time 0.0 should be 2.25.")
        self.assertEqual(state_total_rates[1], 2.25, "Total rate of state 1 at time 100.0 should be 2.25.")


    def test_total_rate_two_epochs(self):
        """
        Test if total rate is correctly computed for two epochs
        """
        pass


    def test_event_sampling_bd(self):
        """
        Test if birth/death events are being randomly drawn proportionally to their weights
        """
        a_time = 0.0
        n_events = n_draws = 10000
        event_outcomes = list()
        while n_draws > 0:
            i = n_events - n_draws
            n_draws -= 1
            total_rate, state_total_rates = self.bd_event_handler.total_rate(a_time, self.bd_state_representation_dict)

            event_outcomes.append(
                self.bd_event_handler.sample_event_atomic_parameter(state_total_rates[0], a_time, [0])[0].event.value
            )

        # first [0] indexing gets lambda parameter, second [0] gets the first item inside the vectorized value
        expected_proportion_event_1 = self.bd_rates_t0_s0[0].value[0] / (self.bd_rates_t0_s0[0].value[0] + self.bd_rates_t0_s0[1].value[0])
        obs_proportion_event_1 = sum(1 for v in event_outcomes if v == 0) / n_events

        # let's see events if really drawn proportional to weights
        # print(expected_proportion_event_1)
        # print(obs_proportion_event_1) # pretty good!
        
        self.assertAlmostEqual(obs_proportion_event_1, expected_proportion_event_1, msg="Proportion of event is off by more than 0.01", delta=0.01)


    def test_event_sampling_bisse(self):
        """
        Test if state-dependent birth/death events and transition events are being
        randomly drawn proportionally to their weights and tree representation
        """
        a_time = 0.0

        n_events = n_draws = 10000
        seeds = [i for i in range(n_events)] # let's set the seed
        event_outcomes_s0, event_outcomes_s1 = list(), list()
        total_rate, state_total_rates = self.bisse_event_handler.total_rate(a_time, self.bisse_state_representation_dict)
        while n_draws > 0:
            i = n_events - n_draws

            n_draws -= 1

            event_outcomes_s0.append(
                self.bisse_event_handler.sample_event_atomic_parameter(state_total_rates[0], a_time, [0], a_seed=seeds[i])[0].event.value
            )

            event_outcomes_s1.append(
                self.bisse_event_handler.sample_event_atomic_parameter(state_total_rates[1], a_time, [1], a_seed=seeds[i])[0].event.value
            )

        # first [0] indexing gets lambda parameter, second [0] gets the first item inside the vectorized value
        expected_proportion_birth_s0 = self.bisse_rates_t0_s0[0].value[0] / state_total_rates[0]
        obs_proportion_birth_s0 = sum(1 for v in event_outcomes_s0 if v == 0) / n_events # v = 0 is the value assigned to speciation (see tpsimulator_classes)

        expected_proportion_birth_s1 = self.bisse_rates_t0_s1[0].value[0] / state_total_rates[1]
        obs_proportion_birth_s1 = sum(1 for v in event_outcomes_s1 if v == 0) / n_events # v = 0 is the value assigned to speciation

        # let's see events if really drawn proportional to weights
        # print(expected_proportion_birth_s0)
        # print(obs_proportion_birth_s0) # pretty good!
        # print(expected_proportion_birth_s1)
        # print(obs_proportion_birth_s1) # pretty good!
        
        self.assertAlmostEqual(obs_proportion_birth_s0, expected_proportion_birth_s0, msg="Proportion of event is off by more than 0.01", delta=0.01)
        self.assertAlmostEqual(obs_proportion_birth_s1, expected_proportion_birth_s1, msg="Proportion of event is off by more than 0.01", delta=0.01)

if __name__ == '__main__':
    # can be called from tests/
    # $ python3 test_events.py
    # 
    # can also be called from phylojunction/
    # $ python3 tests/test_events.py
    # or
    # $ python3 -m tests.test_events
    # or, for a specific test
    # $ python3 -m unittest tests.test_events.TestEvent.test_total_rate_single_epoch_bd
    #
    # can also be called from VS Code, if open folder is phylojuction/
    
    unittest.main()