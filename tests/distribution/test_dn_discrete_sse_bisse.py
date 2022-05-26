import unittest
import math
import statistics

# pj imports
import phylojunction.utility.helper_functions as pjh
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class TestBiSSETrees(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        total_n_states = 2

        rates_t0_s0 = [ sseobj.AtomicSSERateParameter(name="lambda0", val=0.5, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
                        sseobj.AtomicSSERateParameter(name="mu0", val=0.25, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]),
                        sseobj.AtomicSSERateParameter(name="q01", val=0.75, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[0,1]) ]

        rates_t0_s1 = [ sseobj.AtomicSSERateParameter(name="lambda1", val=1.5, event=sseobj.MacroevolEvent.W_SPECIATION, states=[1,1,1]),
                        sseobj.AtomicSSERateParameter(name="mu1", val=0.25, event=sseobj.MacroevolEvent.EXTINCTION, states=[1]),
                        sseobj.AtomicSSERateParameter(name="q10", val=0.75, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[1,0]) ]

        rates_t0 = rates_t0_s0 + rates_t0_s1

        # original implementation
        # matrix_atomic_rate_params = [ [rates_t0_s0, rates_t0_s1] ] # 1D: time slices, 2D: states, 3D: parameters of state, several parameters -> matrix
        matrix_atomic_rate_params = [ rates_t0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

        fig_rates_manager = sseobj.FIGRatesManager(matrix_atomic_rate_params, total_n_states)

        cls.event_handler = sseobj.MacroEvolEventHandler(fig_rates_manager)


    def test_tree_size_state_count_max_taxa_bisse(self):
        """
        Test if BiSSE trees simulated here have similar root ages and number of tips for both states
        as BiSSE trees simulated with diversitree
        """

        stop_condition = "size"
        stop_condition_value = 50 # 50 living taxa

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
                start_states_list=start_states_list, epsilon=1e-12, runtime_limit=3600,
                condition_on_speciation=True, condition_on_survival=True,
                debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)


        # "expectations" from diversitree
        n0_ci_width_maxtaxa_divtree = [
             0.85302550537514, 0.981511843246225, 0.789956054829554,  0.85477208870448, 0.791726103318692, 0.730577291932332, 0.737466005590569, 0.927140583481524, 0.797030131442858, 0.958460078841244,
            0.876910381891682, 0.756412043457255, 0.799218924720947, 0.902042951988357,  0.82727558972556, 0.757368465978368, 0.746478102650882, 0.833724065249121, 0.844229929529298, 0.816000773454261,
            0.752642550859668, 0.850323767081645, 0.816023182913051, 0.838842144124324, 0.775320058633599, 0.709721442761728, 0.853878131334347, 0.836709914561553, 0.875065809491018, 0.863024118588821,
            0.671794086099042, 0.702833666115489, 0.727150641189855,  0.91223569495931, 0.831530420444972, 0.858747041946161, 0.890070577301799, 0.879496633171041, 0.880572652219107, 0.701852364892549,
            0.853701185928911, 0.871388751442404, 0.718639593182936, 0.808251184467417, 0.850710215624417,  0.82061209718528, 0.913220689675005,  0.81224077937762, 0.895848048199738, 0.797077259282284,
            0.879117414634762, 0.796710284839118, 0.793619610267941, 0.786747993003085, 0.940504096301022, 0.817380850839527, 0.743646233127549, 0.892873883913765,  0.80462461978811, 0.765817917490149,
            0.790150922344406, 0.879392280275814, 0.847105731536116, 0.894743991686553, 0.889757647862509, 0.782299009440086, 0.781883890326635, 0.775360949306452, 0.895154318730108, 0.897915705424027,
             0.82893470785287, 0.813424308637474, 0.944264702856052, 0.821307970622894, 0.757241806201082, 0.766645019339678, 0.672673227431425, 0.966408996754125,  0.87739439909474,  0.79974091549382,
            0.829914496677141,  0.73324107979549, 0.862495683539176, 0.854827440585091, 0.771376144015068, 0.798372158525566, 0.977220400330101, 0.758524864509712, 0.891425868225381, 0.859256067051054,
            0.834490840331986, 0.735730564756913, 0.779375405354769, 0.820695557375095, 0.800912524505712, 0.858178494113665, 0.804120090637704, 0.971417774090328, 0.855250829648175, 0.871506660891987
            ]

        n0_mean_maxtaxa_divtree = [
            17.0549450549451, 17.1931818181818, 17.2717391304348, 16.7894736842105, 16.6629213483146, 16.8021978021978, 17.2747252747253, 18.1868131868132, 16.2553191489362, 16.6421052631579,
            16.9381443298969, 16.1978021978022, 17.2934782608696,  17.438202247191, 16.6559139784946,             16.3, 16.5913978494624, 16.9468085106383, 17.4130434782609, 17.4772727272727,
            17.1847826086957, 17.2555555555556, 16.4536082474227, 17.5212765957447, 16.5531914893617,          16.9375, 17.3777777777778, 16.2765957446809,  16.978021978022, 16.8817204301075,
            16.8823529411765, 17.3440860215054, 16.5979381443299,         16.46875, 16.5340909090909, 18.2637362637363, 16.2553191489362, 17.4130434782609, 16.7978723404255,         16.90625,
            17.4777777777778, 16.1145833333333, 16.7634408602151, 16.2395833333333, 17.2903225806452, 17.3586956521739, 17.4838709677419, 17.0215053763441, 17.2340425531915, 17.3406593406593,
            16.5111111111111, 16.5957446808511, 17.3707865168539,             17.0, 16.7553191489362,             17.4, 16.5789473684211, 17.0760869565217, 16.9148936170213, 16.8064516129032,
            16.5280898876404, 16.8478260869565, 16.0434782608696, 16.5161290322581, 17.4479166666667, 16.6813186813187, 16.6129032258065, 17.1978021978022, 17.3258426966292, 16.9574468085106,
            16.6989247311828, 16.6705882352941, 17.032967032967,  17.1521739130435, 16.9222222222222, 16.7659574468085, 17.1847826086957, 16.8181818181818, 16.5955056179775, 16.7340425531915,
            16.9340659340659, 17.2371134020619, 17.4591836734694, 17.5217391304348,             17.0, 16.3020833333333, 16.9574468085106, 17.1505376344086, 16.5483870967742, 17.1666666666667,
            17.1789473684211, 16.6326530612245, 16.4044943820225, 16.8085106382979, 17.2333333333333, 16.8369565217391, 17.6363636363636, 16.8494623655914, 16.6741573033708, 16.6477272727273
            ]

        n1_ci_width_maxtaxa_divtree = [
             0.85302550537514, 0.981511843246225, 0.789956054829554,  0.85477208870448, 0.791726103318692, 0.730577291932332, 0.737466005590569, 0.927140583481524, 0.797030131442858, 0.958460078841244,
            0.876910381891682, 0.756412043457255, 0.799218924720947, 0.902042951988357,  0.82727558972556, 0.757368465978368, 0.746478102650882, 0.833724065249121, 0.844229929529298, 0.816000773454261,
            0.752642550859668, 0.850323767081645, 0.816023182913051, 0.838842144124324, 0.775320058633599, 0.709721442761728, 0.853878131334347, 0.836709914561553, 0.875065809491018, 0.863024118588821,
            0.671794086099042, 0.702833666115489, 0.727150641189855,  0.91223569495931, 0.831530420444972, 0.858747041946161, 0.890070577301799, 0.879496633171041, 0.880572652219107, 0.701852364892549,
            0.853701185928911, 0.871388751442404, 0.718639593182936, 0.808251184467417, 0.850710215624417,  0.82061209718528, 0.913220689675005,  0.81224077937762, 0.895848048199738, 0.797077259282284,
            0.879117414634762, 0.796710284839118, 0.793619610267941, 0.786747993003085, 0.940504096301022, 0.817380850839527, 0.743646233127549, 0.892873883913765,  0.80462461978811, 0.765817917490149,
            0.790150922344406, 0.879392280275814, 0.847105731536116, 0.894743991686553, 0.889757647862509, 0.782299009440086, 0.781883890326635, 0.775360949306452, 0.895154318730108, 0.897915705424027,
             0.82893470785287, 0.813424308637474, 0.944264702856052, 0.821307970622894, 0.757241806201082, 0.766645019339678, 0.672673227431425, 0.966408996754125,  0.87739439909474,  0.79974091549382,
            0.829914496677141,  0.73324107979549, 0.862495683539176, 0.854827440585091, 0.771376144015068, 0.798372158525566, 0.977220400330101, 0.758524864509712, 0.891425868225381, 0.859256067051054,
            0.834490840331986, 0.735730564756913, 0.779375405354769, 0.820695557375095, 0.800912524505712, 0.858178494113665, 0.804120090637704, 0.971417774090328, 0.855250829648175, 0.871506660891987
            ]

        n1_mean_maxtaxa_divtree = [
            32.9450549450549, 32.8068181818182, 32.7282608695652, 33.2105263157895, 33.3370786516854, 33.1978021978022, 32.7252747252747, 31.8131868131868, 33.7446808510638, 33.3578947368421,
            33.0618556701031, 33.8021978021978, 32.7065217391304,  32.561797752809, 33.3440860215054,             33.7, 33.4086021505376, 33.0531914893617, 32.5869565217391, 32.5227272727273,
            32.8152173913044, 32.7444444444444, 33.5463917525773, 32.4787234042553, 33.4468085106383,          33.0625, 32.6222222222222, 33.7234042553191,  33.021978021978, 33.1182795698925,
            33.1176470588235, 32.6559139784946, 33.4020618556701,         33.53125, 33.4659090909091, 31.7362637362637, 33.7446808510638, 32.5869565217391, 33.2021276595745,         33.09375,
            32.5222222222222, 33.8854166666667, 33.2365591397849, 33.7604166666667, 32.7096774193548, 32.6413043478261, 32.5161290322581, 32.9784946236559, 32.7659574468085, 32.6593406593407,
            33.4888888888889, 33.4042553191489, 32.6292134831461,             33.0, 33.2446808510638,             32.6, 33.4210526315789, 32.9239130434783, 33.0851063829787, 33.1935483870968,
            33.4719101123596, 33.1521739130435, 33.9565217391304, 33.4838709677419, 32.5520833333333, 33.3186813186813, 33.3870967741936, 32.8021978021978, 32.6741573033708, 33.0425531914894,
            33.3010752688172, 33.3294117647059,  32.967032967033, 32.8478260869565, 33.0777777777778, 33.2340425531915, 32.8152173913044, 33.1818181818182, 33.4044943820225, 33.2659574468085,
            33.0659340659341, 32.7628865979381, 32.5408163265306, 32.4782608695652,             33.0, 33.6979166666667, 33.0425531914894, 32.8494623655914, 33.4516129032258, 32.8333333333333,
             32.821052631579, 33.3673469387755, 33.5955056179775, 33.1914893617021, 32.7666666666667, 33.1630434782609, 32.3636363636364, 33.1505376344086, 33.3258426966292, 33.3522727272727
            ]

        root_ages_ci_width_maxtaxa_divtree = [
            0.297091385854989, 0.296201258364643, 0.364677547007252, 0.403410726912187,  0.30510452309174, 0.332343398872262, 0.351974400796566, 0.371568927607893, 0.297531811009305, 0.269371236670627,
            0.363329989787168, 0.260701089242488, 0.300001568287857,  0.2530829110793,  0.270014664456289, 0.242136675205298, 0.318485244803122, 0.257715075020624,  0.25907057622389, 0.318105354093215,
            0.295649788296284, 0.311937665024648, 0.305725412362601, 0.332610918147052, 0.345035665033601, 0.346490614568642, 0.370157999432243, 0.316642772699646, 0.358319231492461, 0.285686778286896,
            0.235180477377879, 0.320436587439621, 0.245151314575713, 0.283951344103828, 0.344043808743368, 0.221872924471466, 0.362864547425411, 0.308158207188226, 0.269551937319812, 0.380830627500729,
            0.310512688872746, 0.341875659534818, 0.260369097203599, 0.386563845965079, 0.272009984464708, 0.344875889238449, 0.300178558228157,  0.29843094916527, 0.289186417553936, 0.286954801361071,
            0.321754021568278, 0.281748836173812, 0.256108929712384, 0.307933933303829, 0.325808767321394, 0.312068840058396, 0.347252738581748, 0.332694301661293, 0.275018954210629,  0.30645146122481,
            0.293173460034019, 0.351531190777319, 0.267837033908018, 0.294715004084522, 0.347381154902021, 0.325611688392348, 0.286164044316842, 0.376131282477176, 0.320331348113289, 0.294079255065471,
            0.232423419170186,  0.36972671927173, 0.293116813409044, 0.21141538360236,  0.263866004795028, 0.294293160168103, 0.314429317812875, 0.356704136939503,  0.27014180128031,  0.28571716587361,
             0.43447488617703, 0.279609392718712, 0.306019862379033, 0.264160807235181, 0.283032353344332, 0.306969205178319, 0.325236488829131, 0.310253519049271, 0.272322060747632, 0.255640803405894,
            0.249799371340539, 0.348720429600812, 0.365819330062694, 0.329868016278265, 0.294871195211231, 0.272180825217148, 0.293453252979965, 0.328912047324881, 0.359020707335335, 0.506617695105554
            ]

        root_ages_mean_maxtaxa_divtree = [
             4.1708473048743, 4.38718729944211, 4.41905690583559, 4.52602287188492,  4.2666219495327, 4.37984424449476, 4.26871971605143, 4.48983972903491, 4.28619418582605, 4.10116468196298,
            4.45480456210228, 4.25028477461168, 4.27446268322006, 4.26488454036112, 4.48656571755536, 4.01211346950374, 4.28286176631673, 3.96149423317711, 4.25224318956742, 4.43581870172241,
            4.29646091514782, 4.35719938293082, 4.21829284306175, 4.49034761595346, 4.50895295180056, 4.64584676142023, 4.71785737264991, 4.26107100310806, 4.63898098764023, 3.97510756979318,
            3.99690372012977, 4.39831333657496,  4.2200443518875, 4.20721376766939, 4.37649054423245, 4.07307891361905, 4.54008512360625, 4.50566228766296, 4.22702010754777, 4.47326920370557,
             4.2789224513745, 4.57915846442125, 4.10499158280693, 4.39743731849134, 4.24053408711583, 4.56322982762148, 4.37327926153719, 4.44016579094305, 4.15085289653408, 4.20785589277263,
            4.46968742200277, 4.21507109788349, 4.17899948964936, 4.17313295251671, 4.18242859482202,  4.5222140964834, 4.45805623699317, 4.36607934916459, 4.27093107633833,  4.4339775889521,
            4.35208103425127, 4.53733122889178, 4.27007730098223, 4.38020890972613, 4.34464143922207, 4.62886171599614, 4.34634918937408, 4.73836497112102, 4.45125280505314, 4.16295960077813,
            4.15717661360985, 4.58361703340409, 4.34695224653625, 4.14131065669417, 4.10604689055055, 4.25116267237091,  4.1708907812679, 4.21625227185613, 4.10110237869258, 4.22001633310081,
             4.6806168777176, 4.28886464136169, 4.29760351525763, 4.18350023467351,  4.3888435882561, 4.23094763858926, 4.30406684568677, 4.11357021011849, 3.88989256757541, 4.01962725397889,
             4.2809863046246, 4.45468669353114, 4.28450479861332, 4.46815087019072,  4.4859859281205, 4.28027397975184, 4.44877911001947, 4.32569690095601,   4.366935252278, 4.65710256038698
             ]


        # parsing simulations
        n0_ci_overlap_count = 0
        n1_ci_overlap_count = 0
        root_age_ci_overlap_count = 0
        for i, batch in enumerate(sim_batches):
            n0s = [ann_tr.state_count_dict[0] for ann_tr in batch]
            n1s = [ann_tr.state_count_dict[1] for ann_tr in batch]
            root_ages = [ann_tr.root_age for ann_tr in batch]

            mean_n0 = statistics.mean(n0s)
            mean_n1 = statistics.mean(n1s)
            mean_root_ages = statistics.mean(root_ages)

            stdevs_n0 = statistics.stdev(n0s)
            stdevs_n1 = statistics.stdev(n1s)
            stdevs_root_ages = statistics.stdev(root_ages)

            sterr_n0 = stdevs_n0 / math.sqrt(n_sim)
            sterr_n1 = stdevs_n1 / math.sqrt(n_sim)
            sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)

            n0_ci_width_maxtaxa = 1.96 * sterr_n0
            n1_ci_width_maxtaxa = 1.96 * sterr_n1
            root_ages_ci_width_maxtaxa = 1.96 * sterr_root_ages

            if abs(mean_n0 - n0_mean_maxtaxa_divtree[i]) <= (n0_ci_width_maxtaxa + n0_ci_width_maxtaxa_divtree[i]):
                n0_ci_overlap_count += 1

            if abs(mean_n1 - n1_mean_maxtaxa_divtree[i]) <= (n1_ci_width_maxtaxa + n1_ci_width_maxtaxa_divtree[i]):
                n1_ci_overlap_count += 1

            if abs(mean_root_ages - root_ages_mean_maxtaxa_divtree[i]) <= (root_ages_ci_width_maxtaxa + root_ages_ci_width_maxtaxa_divtree[i]):
                root_age_ci_overlap_count += 1

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n95% CIs of simulations here and from diversitree overlapped " + str(n0_ci_overlap_count) + " times for state 0 count.")
        print("\n95% CIs of simulations here and from diversitree overlapped " + str(n1_ci_overlap_count) + " times for state 1 count.")
        print("\n95% CIs of simulations here and from diversitree overlapped " + str(root_age_ci_overlap_count) + " times for root age.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(n0_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        self.assertAlmostEqual(n1_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        self.assertAlmostEqual(root_age_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)


    def test_tree_size_state_count_max_t_bisse(self):
        """
        Test if BiSSE trees simulated here have similar root ages and number of tips (for both states)
        as BiSSE trees simulated with diversitree
        """

        stop_condition = "age"
        stop_condition_value = 3.0 # 3.0 time units

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
            sse_sim = distsse.DnSSE(self.event_handler, stop_condition_value, n=100, stop=stop_condition, origin=start_at_origin,
                                    start_states_list=start_states_list, epsilon=1e-12, runtime_limit=3600,
                                    condition_on_speciation=True, condition_on_survival=True,
                                    debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)


        # "expectations" from diversitree
        n0_mean_maxt_divtree = [
            4.53, 6.42, 3.86, 4.12, 5.07, 5.33, 4.82, 4.48, 4.09, 6.16,
            5.39, 3.75, 5.42, 4.33, 4.99, 5.47, 4.43, 4.83, 4.62, 5.15,
            4.06, 4.12, 3.99, 4.95, 3.45,  5.3, 5.14, 4.48, 4.86,  5.3,
             4.4, 5.18, 5.43, 5.35, 3.75, 4.73, 4.42, 4.06, 4.96, 4.47,
            4.85, 5.14, 5.21, 4.23, 5.06, 4.31, 4.77, 3.92, 5.43, 5.26,
            4.43, 4.64,    5, 4.91, 4.82, 5.36, 5.23, 5.26,  5.6, 3.65,
            4.17, 4.89, 4.52, 6.19, 4.97, 4.07, 4.75, 3.84, 4.06, 4.76,
            4.17, 4.97, 5.28, 5.23, 4.92, 5.46,  4.4,  5.3, 4.76, 3.58,
            5.55,  4.6, 4.77, 4.46, 4.96, 4.61,  4.7, 4.97, 5.43, 4.34,
            5.19, 4.31, 4.42,  4.6, 3.74,  5.2, 4.17, 3.94, 3.14, 4.34
        ]

        n0_ci_width_maxt_divtree = [
             1.01842872915273,  1.60110418263034, 0.794335708995799, 0.961928077279474,  1.07651592978028,  1.29174936493818, 0.971201880416864,   1.0429102218835,    0.811767413572,  1.20427477249869,
             1.10940645082464,  0.79986349845559,  1.23181884802248, 0.978238376270365,  1.02679204193144,   1.0823396214505,  1.05319389257322,  1.26104411508956,  1.00474394426409, 0.956708055276166,
            0.752079260813283,  1.08227687863144, 0.916994915553797,  1.18859501135168, 0.714809206725964,  1.18798272272578,  1.10383983906318,  1.04662435428854,  1.12955405806938,  1.18176030316702,
            0.709701658758414,  1.11620097100691,  1.00685621414361,  1.38154758780642, 0.871822310688446,  1.00322693163792, 0.900204453544759, 0.980126698880602,   1.0082677085582, 0.913586331692055,
             0.93122139844899,  1.36903275093451,  1.38012001176299,  0.99740817333849,  1.03781523744388, 0.914537263573577, 0.979744572773326, 0.800496347034864,  1.16381223172257,  1.49103206016653,
            0.914774842088025,    1.146771315873,  0.99395717648542,  1.22039862273924,  1.00728969361782,   1.2119193252478,  1.25388483739718,  1.24355705128393,  1.57146082716469, 0.704970033033141,
            0.871234592258293,  1.03221967672105,  1.05143311223061,  1.54774828376548,   1.1460723554996,  1.05393051809391,  1.18237581405261, 0.754088807566436, 0.914165923537722, 0.924842786749197,
            0.894528924951987,  1.03693058396384,  1.50410692506224,  1.19038272158091,  1.11942938999828,  1.25127173061585,  1.11328390382016,  1.17748391226223,  1.08687079250294, 0.790496518050842,
             1.41512526482694, 0.930648258913206,  1.01896201572725, 0.938710642909257,  1.11951258080732, 0.905240213892691,  1.25840887190487,  1.24194431534761,  1.15745988066694,  1.08833361195412,
            0.944182261980378, 0.986391721724624, 0.913044622362041,  1.03488719924034, 0.816307848445881,  1.12196393467307,  1.03152773705957, 0.797164027834126, 0.633442419951442, 0.944398002042903
        ]

        n1_mean_maxt_divtree = [
            7.22, 11.29, 7.15,  8.06, 9.41, 10.09, 8.02,  8.21,  8.32, 11.87,
            8.48,  6.65,  8.7,  7.56,  8.2,  9.36, 8.74,  8.87,  9.75, 10.13,
            8.04,  7.64, 6.81,  7.93, 7.81, 10.21, 8.11,  8.18,  7.61,  9.43,
            8.03,  9.89, 9.67,  9.18, 7.35,   9.4, 6.84,  7.74,  9.09,  6.81,
            8.93,  8.35, 8.61,  7.89, 9.29,  8.86, 8.69,  6.62,  9.53,  8.91,
            7.96,  8.14, 8.53,  8.76, 8.12,  8.75, 8.91,  9.58,  9.74,  6.66,
            7.95,  9.67, 8.69, 10.86, 8.33,  7.97, 9.47,  7.98,  8.32,  9.14,
            8.27,  9.68, 9.54,  9.45, 8.81, 10.28, 7.92, 10.55,  9.62,  7.93,
            9.88,  8.72, 7.72,  9.09,  9.4,  8.43, 8.65,  8.61, 11.41,  7.08,
             9.4,  7.25, 7.87,  8.49, 6.31,  9.42, 8.73,  7.31,   5.7,  8.19
        ]

        n1_ci_width_maxt_divtree = [
            1.94660643991516, 2.79295849340061, 1.50701807379589, 1.78375772418394, 2.23111094092317, 2.22396868628489, 1.90130715099487, 1.96223490763142, 2.22358217716051, 2.66270661447097,
            1.98683655812119, 1.53710149635638, 2.03889108303628, 1.86261363625949, 1.86483728539787, 2.29851813624045, 2.29819735207639, 2.47164729954201, 2.23758537786434, 2.14664857172874,
            1.57907743632128, 2.08164903707021, 1.80191159504667, 1.67479967190411, 1.96251174413265, 2.81592730019793, 1.70798309976736, 2.00702623416996, 1.78399156541015, 2.53773845107423,
            1.54844009665511, 2.32927015195921, 2.23770329985497, 2.45885428719486, 1.88334285865116, 2.34572226175499, 1.55625540943443, 2.51305604900777,    1.98167671848, 1.44827819136931,
            2.08234516239588, 1.92410927050635, 2.72283280263609, 1.81375992263852, 2.06494961280762, 2.15483340338961, 1.61882308614698, 1.64051135462284, 2.28088201653074, 2.82841204518962,
            1.73414274062099, 2.14852139111548, 1.96759858392782, 2.55265892636353, 1.68842000423682, 2.31078492939066,  2.4161443595242, 2.56283272441518, 2.45699769950107, 1.63469170182025,
            1.83428496513855, 2.16436282587283, 2.23989413133656, 2.99071008871597, 2.32045067409525, 2.02321607807616, 2.61577551039994, 1.62702062517581, 2.19459890760306, 1.86400893221946,
            2.15354275030094, 2.39238433992095, 2.81566890959943, 1.99402934033103,  1.9243754601208, 2.35843872890092, 2.11844918921194, 2.30254186287742, 2.22600656374104, 1.88743445176169,
            3.36799857386183, 2.10860831734099, 1.84182281954206,  2.0233388224178, 2.29200588784126,  2.0294971665639, 2.37896454005273, 2.34388450055014, 2.80083889857264, 1.80470573935381,
            1.76212919876248, 2.00857140045035, 2.04595017801243, 2.11236188321428, 1.58979825363872, 2.37321682928212, 2.08059275586492, 1.74813315212994, 1.13759164354998, 2.13977570390691
        ]

        root_ages_mean_maxt_divtree = [
            2.00732372015192, 2.10257928188434, 1.96924293478096, 2.08250892711859, 2.02606497794669,  2.1910383608989,   2.144213932243, 2.05344591114125, 2.10888450947495, 2.13939132596943,
            2.02240137059822, 2.00589399615393, 2.08575353049747, 2.05458060433435, 2.11653820081295, 2.17023732987548, 1.91592946420162, 1.93555734662312, 2.09715568896766, 2.07068452260627,
            1.94676783383699, 2.02640356283249, 1.94946190589675, 2.04421288948422,  1.9994327048228, 2.06766611347784, 2.08509584760876, 2.03834627098419, 2.00618234613033, 2.15702796870466,
            2.03767685320768,  2.1500388013299, 2.16908041892826, 2.05761919687995, 1.96151561929508, 1.99680070282347, 1.99081726249908,  2.0020929100754, 2.09287663009162, 2.07952960521424,
            2.11666124170689,  2.0843141751505, 2.03625917344531, 1.94875032408262, 2.11529827531097, 2.00261143960927,  2.0075434821194, 1.91984485471279, 2.16721120245804, 2.00071729007383,
            2.05601090136388, 2.06827967482995, 2.01059498587899, 2.06315000226629, 2.06143741551062,  1.9953055160624, 1.93319966731066, 2.04575552601094, 1.95284027914759, 1.90526435624664,
            2.01664689527347,  2.0020219815412, 1.91649329880205, 2.16094830412944, 2.02470068767687, 1.93669068777118, 1.99753587691501, 1.99380890942962,  1.9595439348218, 2.13792713239682,
            1.95396613994661, 1.87821353516645, 2.09743180648574, 2.09572419802156, 2.15962745109445, 2.19517653866974, 2.01542282351472,  2.0439752160935,  2.1361131092323, 1.92635463251225,
            2.06472685735143, 1.93991676414349, 2.11393395125013,  2.1018255016506, 2.09268604725736, 2.04155509875933, 1.89302471893026, 2.04282570807596, 2.13749597206643, 1.95231362180622,
            2.12001314945457, 2.02511639478348,  1.9267504893811, 2.11211125949106,  2.0050087749629, 2.04797546625838, 1.97387690374822, 2.04193431544719, 2.04275906552701, 2.07206553135074
        ]

        root_ages_ci_width_maxt_divtree = [
            0.297091385854989, 0.296201258364643, 0.364677547007252, 0.403410726912187,  0.30510452309174, 0.332343398872262, 0.351974400796566, 0.371568927607893, 0.297531811009305, 0.269371236670627,
            0.363329989787168, 0.260701089242488, 0.300001568287857,   0.2530829110793, 0.270014664456289, 0.242136675205298, 0.318485244803122, 0.257715075020624,  0.25907057622389, 0.318105354093215,
            0.295649788296284, 0.311937665024648, 0.305725412362601, 0.332610918147052, 0.345035665033601, 0.346490614568642, 0.370157999432243, 0.316642772699646, 0.358319231492461, 0.285686778286896,
            0.235180477377879, 0.320436587439621, 0.245151314575713, 0.283951344103828, 0.344043808743368, 0.221872924471466, 0.362864547425411, 0.308158207188226, 0.269551937319812, 0.380830627500729,
            0.310512688872746, 0.341875659534818, 0.260369097203599, 0.386563845965079, 0.272009984464708, 0.344875889238449, 0.300178558228157,  0.29843094916527, 0.289186417553936, 0.286954801361071,
            0.321754021568278, 0.281748836173812, 0.256108929712384, 0.307933933303829, 0.325808767321394, 0.312068840058396, 0.347252738581748, 0.332694301661293, 0.275018954210629,  0.30645146122481,
            0.293173460034019, 0.351531190777319, 0.267837033908018, 0.294715004084522, 0.347381154902021, 0.325611688392348, 0.286164044316842, 0.376131282477176, 0.320331348113289, 0.294079255065471,
            0.232423419170186,  0.36972671927173, 0.293116813409044,  0.21141538360236, 0.263866004795028, 0.294293160168103, 0.314429317812875, 0.356704136939503,  0.27014180128031,  0.28571716587361,
             0.43447488617703, 0.279609392718712, 0.306019862379033, 0.264160807235181, 0.283032353344332, 0.306969205178319, 0.325236488829131, 0.310253519049271, 0.272322060747632, 0.255640803405894,
            0.249799371340539, 0.348720429600812, 0.365819330062694, 0.329868016278265, 0.294871195211231, 0.272180825217148, 0.293453252979965, 0.328912047324881, 0.359020707335335, 0.506617695105554
        ]

        # parsing simulations
        n0_ci_overlap_count = 0
        n1_ci_overlap_count = 0
        root_age_ci_overlap_count = 0
        for i, batch in enumerate(sim_batches):
            n0s = [ann_tr.state_count_dict[0] for ann_tr in batch]
            n1s = [ann_tr.state_count_dict[1] for ann_tr in batch]
            root_ages = [ann_tr.root_age for ann_tr in batch]

            mean_n0 = statistics.mean(n0s)
            mean_n1 = statistics.mean(n1s)
            mean_root_ages = statistics.mean(root_ages)

            stdevs_n0 = statistics.stdev(n0s)
            stdevs_n1 = statistics.stdev(n1s)
            stdevs_root_ages = statistics.stdev(root_ages)

            sterr_n0 = stdevs_n0 / math.sqrt(n_sim)
            sterr_n1 = stdevs_n1 / math.sqrt(n_sim)
            sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)

            n0_ci_width_maxt = 1.96 * sterr_n0
            n1_ci_width_maxt = 1.96 * sterr_n1
            root_ages_ci_width_maxt = 1.96 * sterr_root_ages

            if abs(mean_n0 - n0_mean_maxt_divtree[i]) <= (n0_ci_width_maxt + n0_ci_width_maxt_divtree[i]):
                n0_ci_overlap_count += 1

            if abs(mean_n1 - n1_mean_maxt_divtree[i]) <= (n1_ci_width_maxt + n1_ci_width_maxt_divtree[i]):
                n1_ci_overlap_count += 1

            if abs(mean_root_ages - root_ages_mean_maxt_divtree[i]) <= (root_ages_ci_width_maxt + root_ages_ci_width_maxt_divtree[i]):
                root_age_ci_overlap_count += 1

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n95% CIs of simulations here and from diversitree overlapped " + str(n0_ci_overlap_count) + " times for state 0 count.")
        print("\n95% CIs of simulations here and from diversitree overlapped " + str(n1_ci_overlap_count) + " times for state 1 count.")
        print("\n95% CIs of simulations here and from diversitree overlapped " + str(root_age_ci_overlap_count) + " times for root age.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(n0_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        self.assertAlmostEqual(n1_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        self.assertAlmostEqual(root_age_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)

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
    # $ python3 tests/distribution/test_dn_discrete_sse_bisse.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_bisse
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_bisse.TestBDTrees.test_tree_size_state_count_max_taxa_bisse

    unittest.main()