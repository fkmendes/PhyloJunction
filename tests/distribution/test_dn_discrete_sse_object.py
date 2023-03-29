import unittest

# pj imports
import phylojunction.utility.helper_functions as pjh
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class TestDnSSEObject(unittest.TestCase):

    def test_dnsse_vectorization(self):
        """Test if DnSSE is vectorized
        
        First 5 trees must be a Yule trees (l = 1.0, mu = 0.0)
        Second 5 trees must not speciate (l = 0.0, mu = 10.0)
        """

        #########################
        # Birth-death ingredients #
        #########################
        total_n_states = 2

        l = [ 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
        lrate = sseobj.DiscreteStateDependentRate(name="lambda0", val=l, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0])
        
        mu = [ 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0, 10.0, 10.0 ]
        murate = sseobj.DiscreteStateDependentRate(name="mu0", val=mu, event=sseobj.MacroevolEvent.EXTINCTION, states=[0])
        
        rates_t0 = [ lrate, murate ]

        matrix_atomic_rate_params = [ rates_t0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice
        
        state_dep_par_manager = sseobj.DiscreteStateDependentParameterManager(matrix_atomic_rate_params, total_n_states)

        meh = sseobj.MacroevolEventHandler(state_dep_par_manager)

        ########
        # Tree #
        ########
        n_sim = 10
        stop_condition = "age"
        stop_condition_value = [4.0] # 4.0 time units
        start_at_origin = True
        start_states_list = [0 for i in range(n_sim)]

        sse_sim = distsse.DnSSE(n=n_sim, stop=stop_condition, stop_value=stop_condition_value, origin=start_at_origin, event_handler=meh,
                start_states_list=start_states_list, epsilon=1e-12, runtime_limit=3600,
                condition_on_speciation=False, condition_on_survival=False,
                debug=False)

        trs = sse_sim.generate()

        trs1to5 = [ tr for tr in trs[:5] ]
        trs6to10 = [ tr for tr in trs[5:] ]

        for tr in trs1to5:
            # first 5 trees are Yule trees and cannot have died
            self.assertFalse(tr.tree_died)
        for tr in trs6to10:
            # second 5 trees must have died
            self.assertTrue(tr.tree_died)

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
    # $ python3 tests/distribution/test_dn_discrete_sse_object.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_object
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_object.TestDnSSEObject.test_dnsse_vectorization

    unittest.main()