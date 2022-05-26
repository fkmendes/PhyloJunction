import unittest
import pandas as pd
import io

# pj imports
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.pgm.pgm as pgm
import phylojunction.readwrite.pj_write as pjwrite
import phylojunction.distribution.dn_discrete_sse as dnsse

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class TestDataDump(unittest.TestCase):

    def test_bisse_data_dump(self):
        """Test if data is correctly outputted from pgm"""

        #########################
        # BiSSE pgm ingredients #
        #########################
        total_n_states = 2

        l0 = [ 1.0, 1.1, 0.9, 0.95, 1.05 ]
        l0rate = sseobj.AtomicSSERateParameter(name="lambda0", val=l0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0])

        mu0 = [ 0.23, 0.24, 0.25, 0.26, 0.27 ]
        mu0rate = sseobj.AtomicSSERateParameter(name="mu0", val=mu0, event=sseobj.MacroevolEvent.EXTINCTION, states=[0])
        
        q01rate = sseobj.AtomicSSERateParameter(name="q01", val=0.75, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[0,1])
        
        rates_t0_s0 = [ l0rate, mu0rate, q01rate ]

        l1rate = sseobj.AtomicSSERateParameter(name="lambda1", val=1.5, event=sseobj.MacroevolEvent.W_SPECIATION, states=[1,1,1])
        mu1rate = sseobj.AtomicSSERateParameter(name="mu1", val=0.25, event=sseobj.MacroevolEvent.EXTINCTION, states=[1])
        q10rate = sseobj.AtomicSSERateParameter(name="q10", val=0.75, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[1,0])
        
        rates_t0_s1 = [ l1rate, mu1rate, q10rate ]
        
        rates_t0 = rates_t0_s0 + rates_t0_s1

        matrix_atomic_rate_params = [ rates_t0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice
        
        fig_rates_manager = sseobj.FIGRatesManager(matrix_atomic_rate_params, total_n_states)

        meh = sseobj.MacroEvolEventHandler(fig_rates_manager)
        
        ########
        # Tree #
        ########
        n_sim = 5
        n_repl = 2
        stop_condition = "size"
        stop_condition_value = 50 # 50 living taxa
        start_at_origin = True
        start_states_list = [0 for i in range(n_sim)]

        sse_sim = dnsse.DnSSE(meh, stop_condition_value, n=n_sim, n_replicates=n_repl, stop=stop_condition, origin=start_at_origin,
                start_states_list=start_states_list, epsilon=1e-12, runtime_limit=3600,
                condition_on_speciation=True, condition_on_survival=True,
                debug=False)

        trs = sse_sim.generate()

        ####################
        # Initializing PGM #
        ####################
        bisse_pgm = pgm.ProbabilisticGraphicalModel()

        # rv
        bisse_pgm.add_node(pgm.StochasticNodePGM("l0", n_sim, value=l0, sampled_from="Log-normal"))
        bisse_pgm.add_node(pgm.StochasticNodePGM("mu0", n_sim, value=mu0, sampled_from="Log-normal"))
        
        # deterministic
        bisse_pgm.add_node(pgm.DeterministicNodePGM("l0r", value=l0rate))
        bisse_pgm.add_node(pgm.DeterministicNodePGM("mu0r", value=mu0rate))
        bisse_pgm.add_node(pgm.DeterministicNodePGM("q01r", value=q01rate))
        bisse_pgm.add_node(pgm.DeterministicNodePGM("l1r", value=l1rate))
        bisse_pgm.add_node(pgm.DeterministicNodePGM("mu1r", value=mu1rate))
        bisse_pgm.add_node(pgm.DeterministicNodePGM("q10r", value=q10rate))
        bisse_pgm.add_node(pgm.DeterministicNodePGM("meh", value=meh))

        # more rv
        bisse_pgm.add_node(pgm.StochasticNodePGM("trs", n_sim, value=trs, sampled_from="DnSSE", replicate_size=n_repl))

        sorted_node_pgm_list = bisse_pgm.get_sorted_node_pgm_list()
            
        ###################
        # Output handling #
        ###################
        # populating dataframe to be dumped
        data_df_names_list, data_df_list = pjwrite.prep_data_df(sorted_node_pgm_list)

        scalar_data_df, tree_data_df = data_df_list

        # writing to file handle
        scalar_data_outfile = io.StringIO()
        pjwrite.write_data_df(scalar_data_outfile, scalar_data_df)
        scalar_data_outfile.seek(0) # move to the start of the file handle

        # getting content of file handle
        csvstring_containg_data_dump = scalar_data_outfile.read()

        # checking file name suffixes are ok
        self.assertEqual(["_rvs.csv", "_trs.tsv"], data_df_names_list)
        # checking dump
        self.assertEqual(csvstring_containg_data_dump, ",simulation,l0,mu0\n0,0,1.0,0.23\n1,1,1.1,0.24\n2,2,0.9,0.25\n3,3,0.95,0.26\n4,4,1.05,0.27\n")
        self.assertEqual(tree_data_df.shape, (10, 3))
        self.assertEqual("\n".join(str(i) for i in tree_data_df["simulation"]), "0\n0\n1\n1\n2\n2\n3\n3\n4\n4")
        self.assertEqual("\n".join(str(i) for i in tree_data_df["replicate"]), "0\n1\n0\n1\n0\n1\n0\n1\n0\n1")

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
    # $ python3 tests/readwrite/test_data_dump.py
    # 
    # or
    #
    # $ python3 -m tests.readwrite.test_data_dump
    #
    # or 
    #
    # $ python3 -m unittest tests.readwrite.test_data_dump.TestDataDump.test_bisse_data_dump

    unittest.main()