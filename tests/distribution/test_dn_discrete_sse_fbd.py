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

class TestFBDTrees(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # not state-dependent (just state 0, and no transition)
        total_n_states = 1

        # not state-dependent (just state 0, and no transition)
        rates_t0_s0 = [ sseobj.MacroevolStateDependentRateParameter(name="lambda", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
                        sseobj.MacroevolStateDependentRateParameter(name="mu", val=0.8, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]),
                        sseobj.MacroevolStateDependentRateParameter(name="psi", val=0.5, event=sseobj.MacroevolEvent.ANCESTOR_SAMPLING, states=[0]) ]

        matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

        fig_rates_manager = sseobj.FIGRatesManager(matrix_atomic_rate_params, total_n_states)

        cls.event_handler = sseobj.MacroevolEventHandler(fig_rates_manager)

    
    def test_tree_size_n_taxa_fbd(self):
        """
        Test if fossilized birth-death (FBD) trees simulated here have similar root ages,
        and number of tips and sampled ancestors as FBD trees simulated with FossilSim
        """
        stop_condition = "size"
        stop_condition_value = [ 50 ] # 50 living taxa

        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]

        # "expectations" from FossilSim
        n_sa_taxa_ci_width_maxtaxa_fossilsim = [
             1.04144512960121, 0.971092243225294, 0.994036066332127,  1.00002422782642,  1.10880886016191, 0.961423922289011,  1.13625670985387,  0.91230175174119,  1.11420040946364, 0.798022445310581,
             1.09304194687533, 0.892768394535425, 0.935465626311146,  1.17288070176979, 0.808028169701806, 0.932369680107502,  1.07334052301731,  1.12536054977677, 0.958558047590096,  1.07061793282086,
            0.909666925821425,  1.00244721303609, 0.953655982288645,  1.18038895194484,  1.04368073156611,  1.06320379948784, 0.972457924363392, 0.797408133347388,  1.05067995059605,  1.11708435998285,
             1.02961717009679,  1.33570250719431,  1.06585647621065, 0.903760970710463,  1.03493427219642,  1.25806973064216, 0.930249941381296, 0.888552728485135,  1.09176762056951,  1.08291704514492,
             1.08551009189828, 0.934365438552827, 0.992373210181591,  1.11768971091573,  1.00959157601775, 0.985748587729568, 0.903490950840652,   1.0383612370508, 0.968405645214411,  1.04091883431209,
            0.933088637322334,  1.10625340269267, 0.975537632136329,    1.003631580758,  1.08404827826051,  0.98947395265709, 0.773418984666904, 0.972006290692419,  1.26793129215961,  1.01821110498232,
             1.09089955746752, 0.934258652913357, 0.936365792357021, 0.926440659311644, 0.954511991226026,  1.07696570699701, 0.909499275184225, 0.966014279361485,  1.10503395970259, 0.958119799649725,
             1.05464661823498, 0.870720056371255,   1.1335962753906,  1.32135996630758,  1.04364109256917,  1.08540163594579, 0.760901727871346,   0.9042811858493, 0.905134286817805,  0.99903527314115,
            0.839226320707002, 0.925694287589084, 0.923126676486251, 0.987710196153474, 0.980410623354323, 0.953108079721275,  1.12740241667977,  1.04420156358422, 0.894646497548185, 0.997290863275285,
            0.811449528537197,  1.02125192918142,  1.05284243033448, 0.960979109697982, 0.919428067079735,  1.13161935196259,  1.04965147602406,  0.95939713750412, 0.994860868909459,  1.33529316332578
            ]

        n_sa_taxa_mean_maxtaxa_fossilsim = [
            24.33, 24.09, 25.76, 25.67, 25.87, 21.26, 30.92, 21.86, 25.27, 22.33,
            26.26, 25.93, 22.91,  26.5, 21.97, 22.36, 25.21, 25.72, 25.57, 27.15,
            21.71, 23.11, 21.63, 27.61, 27.88, 25.93, 22.85, 20.99, 26.19,  25.6,
            25.17, 25.86, 27.17, 22.53, 26.02, 27.41,  24.7, 23.67, 23.47, 27.21,
            23.06, 24.78, 20.43,  22.5,  23.5, 26.56, 26.97, 27.55, 22.14, 27.09,
             23.6, 26.68, 26.27,  28.2, 23.92, 21.64, 20.47, 26.71, 30.71, 26.49,
            27.27, 23.19, 23.85,  22.7, 23.54,  26.8, 25.84, 28.55, 26.58, 20.45,
            23.51, 22.91, 23.93, 29.03, 25.33,  28.6, 23.47, 24.13, 21.71, 22.79,
            28.23,  23.6, 23.66, 25.29, 26.58, 25.69, 23.72, 29.33, 23.11, 26.45,
            20.84, 21.94, 26.71, 19.43,  24.8, 26.98, 25.38, 26.74, 22.65, 27.33
            ]

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
    # $ python3 tests/distribution/test_dn_discrete_sse_fbd.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_fbd
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_fbd.TestFBDTrees.test_tree_size_n_taxa_fbd

    total_n_states = 1

    # not state-dependent (just state 0, and no transition)
    rates_t0_s0 = [ sseobj.MacroevolStateDependentRateParameter(name="lambda", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
                    sseobj.MacroevolStateDependentRateParameter(name="mu", val=0.5, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]),
                    sseobj.MacroevolStateDependentRateParameter(name="psi", val=0.75, event=sseobj.MacroevolEvent.ANCESTOR_SAMPLING, states=[0]) ]

    matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

    fig_rates_manager = sseobj.FIGRatesManager(matrix_atomic_rate_params, total_n_states)

    event_handler = sseobj.MacroevolEventHandler(fig_rates_manager)

    # stop_condition = "size"
    stop_condition = "age"
        
    # stop_condition_value = [ 10 ] # 10 living taxa
    stop_condition_value = [ 3.0 ] # origin age 

    start_at_origin = True

    # simulation initialization
    n_sim = 100
    start_states_list = [0 for i in range(n_sim)]

    sse_sim = distsse.DnSSE(
        event_handler,
        stop_condition_value,
        n=n_sim,
        stop=stop_condition,
        origin=start_at_origin,
        start_states_list=start_states_list,
        epsilon=1e-12,
        runtime_limit=3600,
        condition_on_speciation=True,
        condition_on_survival=True,
        debug=False)

    trs = sse_sim.generate()

    # print(trs[0].tree.as_string(schema="newick"))

    n_sa = float()
    n_leaves = float()
    root_age = float()
    origin_age = float()
    for tr in trs:
        # print(tr.tree.as_string(schema="newick"))
        # print(tr.n_sa)
        # print(str(tr.root_age) + "\n")
        n_sa += tr.n_sa
        n_leaves += tr.n_extant_terminal_nodes
        root_age += tr.root_age
        origin_age += tr.origin_age
    
    print("mean sa = " + str(n_sa / n_sim))
    print("mean leaves = " + str(n_leaves / n_sim))
    print("mean root age = " + str(root_age / n_sim))
    print("mean origin age = " + str(origin_age / n_sim))

    print("mean leaves / mean sa = " + str(n_leaves / n_sa))
    print("mean sa / mean root age = " + str(n_sa / root_age))

    # unittest.main()