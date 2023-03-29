import unittest

# pj imports
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.utility.exception_classes as ec


__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class TestDiscreteSSE(unittest.TestCase):
    
    def test_state_dep_prob(self):
        """
        Test that state-dependent probability parameter is correctly
        initialized
        """

        state_dep_sampling_prob_A = \
            sseobj.DiscreteStateDependentProbability(
                name="rho0_A", val=0.5, state=0
            )

        state_dep_sampling_prob_B = \
            sseobj.DiscreteStateDependentProbability(
                name="rho0_B", val=[0.5, 0.5], state=1
            )

        state_dep_sampling_prob_C = \
            sseobj.DiscreteStateDependentProbability(
                name="rho0_C", val=["0.5", "0.5"], state=0
            )

        expected_str_A = "Discrete state-dependent probability\n" \
            + "   Name:   rho0_A\n" \
            + "   Value:  0.5\n" \
            + "   State:  0\n\n"

        expected_str_B = "Discrete state-dependent probability\n" \
            + "   Name:   rho0_B\n" \
            + "   Value:  0.5, 0.5\n" \
            + "   State:  1\n\n"

        expected_str_C = "Discrete state-dependent probability\n" \
            + "   Name:   rho0_C\n" \
            + "   Value:  0.5, 0.5\n" \
            + "   State:  0\n\n"

        self.assertEqual(str(state_dep_sampling_prob_A),
            expected_str_A)

        self.assertEqual(str(state_dep_sampling_prob_B),
            expected_str_B)

        self.assertEqual(str(state_dep_sampling_prob_C),
            expected_str_C)
            

    def test_state_dep_prob_init_error(self):
        """
        Test that state-dependent probability parameter throws
        exception if value(s) < 0.0 or value(s) > 1.0
        """
        
        with self.assertRaises(ec.NotBetweenZeroAndOneError) as exc1:
            sseobj.DiscreteStateDependentProbability(
                name="rho0_A", val=-0.1, state=0
            )

        with self.assertRaises(ec.NotBetweenZeroAndOneError) as exc2:
            sseobj.DiscreteStateDependentProbability(
                name="rho0_B", val=1.1, state=0
            )

        with self.assertRaises(ec.NotBetweenZeroAndOneError) as exc3:
            sseobj.DiscreteStateDependentProbability(
                name="rho0_C", val=[0.1, -0.1], state=0
            )

        with self.assertRaises(ec.NotBetweenZeroAndOneError) as exc4:
            sseobj.DiscreteStateDependentProbability(
                name="rho0_D", val=[1.0, 1.1], state=0
            )

        self.assertEqual(str(exc1.exception),
            "ERROR: 'rho0_A' must be contained in [0,1], but one or " \
                + "more of its values was negative."
        )

        self.assertEqual(str(exc2.exception),
            "ERROR: 'rho0_B' must be contained in [0,1], but one or " \
                + "more of its values was positive."
        )

        self.assertEqual(str(exc3.exception),
            "ERROR: 'rho0_C' must be contained in [0,1], but one or " \
                + "more of its values was negative."
        )

        self.assertEqual(str(exc4.exception),
            "ERROR: 'rho0_D' must be contained in [0,1], but one or " \
                + "more of its values was positive."
        )


    def test_state_dep_param_manager(self):
        """
        """
        
        total_n_states = 2

        probs_t0 = [
            sseobj.DiscreteStateDependentProbability(name="rho0_0", val=0.1, state=0),
            sseobj.DiscreteStateDependentProbability(name="rho0_1", val=0.2, state=1),
        ]

        probs_t1 = [
            sseobj.DiscreteStateDependentProbability(name="rho1_0", val=0.2, state=0),
            sseobj.DiscreteStateDependentProbability(name="rho1_1", val=0.4, state=1),
        ]

        # 1D: time slices (i)
        # 2D: all rates from all states in i-th time slice
        matrix_state_dep_probs = [ probs_t0, probs_t1 ]

        state_dep_par_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_state_dep_probs, total_n_states,
                seed_age_for_time_slicing=2.0,
                list_time_slice_age_ends=[1.0]
            )

        params_matrix_state0 = \
            state_dep_par_manager.state_dep_params_dict[0]
        params_matrix_state1 = \
            state_dep_par_manager.state_dep_params_dict[1]
        
        params_state0_at_time11 = \
            state_dep_par_manager.state_dep_params_at_time(
                params_matrix_state0, 1.1)
        params_state1_at_time11 = \
            state_dep_par_manager.state_dep_params_at_time(
                params_matrix_state1, 1.1)

        params_state0_at_time09 = \
            state_dep_par_manager.state_dep_params_at_time(
                params_matrix_state0, 0.9)
        params_state1_at_time09 = \
            state_dep_par_manager.state_dep_params_at_time(
                params_matrix_state1, 0.9)
        
        self.assertEqual(params_state0_at_time11[0].name, "rho1_0")
        self.assertEqual(params_state1_at_time11[0].name, "rho1_1")
        self.assertEqual(params_state0_at_time09[0].name, "rho0_0")
        self.assertEqual(params_state1_at_time09[0].name, "rho0_1")


    def test_state_dep_param_manager_init_error(self):
        """
        """
        
        total_n_states = 2

        probs_t0 = [
            sseobj.DiscreteStateDependentProbability(name="rho0_0", val=0.1, state=0),
            sseobj.DiscreteStateDependentProbability(name="rho0_1", val=0.2, state=1),
        ]

        probs_t1 = [
            sseobj.DiscreteStateDependentProbability(name="rho1_0", val=0.2, state=0),
            # sseobj.DiscreteStateDependentProbability(name="rho1_1", val=0.4, state=1),
        ]

        prob_rate_t1 = [
            sseobj.DiscreteStateDependentProbability(name="rho1_0", val=0.2, state=0),
            sseobj.DiscreteStateDependentRate(name="lambda1_0", val=1.0,
            event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0])   
        ]

        # 1D: time slices (i)
        # 2D: all rates from all states in i-th time slice
        matrix_state_dep_probs_error1 = [ probs_t0, probs_t1 ]
        matrix_state_dep_probs_error2 = [ probs_t0, prob_rate_t1 ]

        # number of time slices deduced from matrix of
        # parameters does not match default number of
        # slices (as a result of not explicitly passing
        # a number of time slices
        with self.assertRaises(ec.WrongDimensionError) as exc1:
            state_dep_par_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_state_dep_probs_error1, total_n_states
            )

        # rate and prob both passed together
        # but only one type allowed
        with self.assertRaises(ec.RequireSameParameterType) as exc1:
            state_dep_par_manager = \
            sseobj.DiscreteStateDependentParameterManager(
                matrix_state_dep_probs_error2, total_n_states,
                seed_age_for_time_slicing=2.0,
                list_time_slice_age_ends=[1.0]
            )

        # TODO:
        # add test for one of the slices missing a parameter
        # associated with one of the states


    def test_state2pattern_converter(self):
        """
        Test if conversion between integers representing compound
        states and bit patterns is correct
        """
        
        n_characters = 3 # regions A, B, C
        n_states_per_char = 2 # presence/absence (bit)
        # compound-states integer coding -> bit pattern
        # Null: 0 -> 000
        # A:    1 -> 100
        # B:    2 -> 010
        # C:    3 -> 001
        # AB:   4 -> 110
        # AC:   5 -> 101
        # BC:   6 -> 011
        # ABC:  7 -> 111
        
        svc = sseobj.StateIntoPatternConverter(
            n_characters, n_states_per_char)
        
        self.assertEqual(list(svc.int2set_dict.values()),
            ["000", "100", "010", "001", "110", "101", "011", "111"])
        self.assertEqual(list(svc.set2int_dict.values()),
        ["0", "1", "2", "3", "4", "5", "6", "7"])


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
    # on the terminal, remember to add "src/" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3 tests/calculation/test_discrete_sse.py
    # 
    # or
    #
    # $ python3.9 -m tests.calculation.test_discrete_sse
    #
    # or 
    #
    # $ python3.9 -m unittest tests.calculation.test_discrete_sse.TestDiscreteSSE.test_state_dep_prob

    unittest.main()