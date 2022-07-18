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

        rates_t0_s0 = [ sseobj.MacroevolStateDependentRateParameter(name="lambda0", val=0.5, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
                        sseobj.MacroevolStateDependentRateParameter(name="mu0", val=0.25, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]),
                        sseobj.MacroevolStateDependentRateParameter(name="q01", val=0.8, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[0,1]) ]

        rates_t0_s1 = [ sseobj.MacroevolStateDependentRateParameter(name="lambda1", val=1.5, event=sseobj.MacroevolEvent.W_SPECIATION, states=[1,1,1]),
                        sseobj.MacroevolStateDependentRateParameter(name="mu1", val=0.25, event=sseobj.MacroevolEvent.EXTINCTION, states=[1]),
                        sseobj.MacroevolStateDependentRateParameter(name="q10", val=0.4, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[1,0]) ]

        rates_t0 = rates_t0_s0 + rates_t0_s1

        # original implementation
        # matrix_atomic_rate_params = [ [rates_t0_s0, rates_t0_s1] ] # 1D: time slices, 2D: states, 3D: parameters of state, several parameters -> matrix
        matrix_atomic_rate_params = [ rates_t0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

        fig_rates_manager = sseobj.FIGRatesManager(matrix_atomic_rate_params, total_n_states)

        cls.event_handler = sseobj.MacroevolEventHandler(fig_rates_manager)


    def test_tree_size_state_count_max_taxa_bisse(self):
        """
        Test if BiSSE trees simulated here have similar root ages and number of tips for both states
        as BiSSE trees simulated with diversitree
        """
        stop_condition = "size"
        stop_condition_value = [ 50 ] # 50 living taxa

        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]

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
        n0_ci_width_maxtaxa_divtree = [ 0.69857566176268, 0.682181290516907, 0.76148707238108, 0.659864747076619, 0.660451229917595, 0.630322020627465, 0.679375574025315, 0.684939970477629, 0.626276111252214, 0.650732224569865, 0.691289573122707, 0.649051786282571, 0.615812527367672, 0.621703633507002, 0.626700228233863, 0.751828202108001, 0.797753519987023, 0.75695689206503, 0.657415436891266, 0.66893069513576, 0.736935133715545, 0.627674339863497, 0.772694255427808, 0.653883277313985, 0.694208227682342, 0.697388865373198, 0.701924381142374, 0.636653681618197, 0.690857929827726, 0.716310315594732, 0.593933750238816, 0.594928872754016, 0.690287604255378, 0.750791185023354, 0.680515481470936, 0.665889901848514, 0.68140133801784, 0.662532365313363, 0.780446942599069, 0.615167699286637, 0.68518592047271, 0.762936034251056, 0.667501463300847, 0.671750693886387, 0.747830984224012, 0.690662847114914, 0.700209952499606, 0.758270095982651, 0.660436766804131, 0.745837642837692, 0.717507842121821, 0.692333307688362, 0.679680761934844, 0.783311978805604, 0.722241819228017, 0.669684615144087, 0.619384423671881, 0.708978899757948, 0.69885619660113, 0.698526395169173, 0.698908376898035, 0.646233816520023, 0.524284895554321, 0.639301885443374, 0.671639907716296, 0.666481333551904, 0.640899772897716, 0.697880916199803, 0.735296979117167, 0.624453167632793, 0.687693785247505, 0.670464857205835, 0.702839809119529, 0.609142738817578, 0.580027453365198, 0.63124339537534, 0.64105061655259, 0.745641826989972, 0.737839532186865, 0.609461707309304, 0.62555104675997, 0.644266345010154, 0.639762510271821, 0.701369207891536, 0.725318739884268, 0.599698564191191, 0.674392208484714, 0.686387124811052, 0.723609362978651, 0.735629372876929, 0.682255794168838, 0.730028893652169, 0.609284860502141, 0.636888137305202, 0.625368743182459, 0.680609108546504, 0.724767723226951, 0.682673007331043, 0.754500441360665, 0.718855177954899 ]

        n0_mean_maxtaxa_divtree = [ 9.63157894736842, 10.2626262626263, 10.3541666666667, 9.75789473684211, 9.91578947368421, 9.47777777777778, 9.91011235955056, 10.09375, 9.48421052631579, 9.70652173913044, 10.0842105263158, 9.6, 9.69148936170213, 9.36170212765957, 9.41758241758242, 10.4782608695652, 9.86666666666667, 9.43617021276596, 9.34736842105263, 9.93548387096774, 9.1304347826087, 9.75, 10.2755102040816, 9.83333333333333, 9.66315789473684, 9.94623655913978, 9.78125, 9.71276595744681, 9.78409090909091, 10.3695652173913, 9.56989247311828, 10.030303030303, 10.445652173913, 9.64583333333333, 9.65957446808511, 9.81720430107527, 9.93258426966292, 9.42708333333333, 10.7444444444444, 9.45833333333333, 10.0104166666667, 9.56382978723404, 9.95833333333333, 9.46875, 9.74468085106383, 9.45360824742268, 9.8936170212766, 9.63736263736264, 10.1808510638298, 9.0625, 10.0322580645161, 9.43298969072165, 9.52040816326531, 9.83333333333333, 10.5161290322581, 9.75824175824176, 10.0531914893617, 9.76595744680851, 9.73404255319149, 10.7473684210526, 10.1458333333333, 10.0520833333333, 9.26881720430108, 10.0681818181818, 9.30851063829787, 10.0315789473684, 9.78350515463918, 10.1122448979592, 9.76344086021505, 9.69473684210526, 10.5670103092784, 10.1473684210526, 9.69318181818182, 10.4782608695652, 9.78494623655914, 9.40425531914894, 9.69791666666667, 9.97872340425532, 9.74489795918367, 9.72631578947368, 10.0851063829787, 9.76595744680851, 9.76595744680851, 9.61290322580645, 9.8494623655914, 10.3913043478261, 10.1397849462366, 9.6, 10.1368421052632, 9.69148936170213, 9.94623655913978, 9.40816326530612, 9.91836734693878, 9.68421052631579, 9.89473684210526, 10.3510638297872, 9.91578947368421, 9.77894736842105, 9.89473684210526, 10.0106382978723 ]

        n1_ci_width_maxtaxa_divtree = [ 0.69857566176268, 0.682181290516907, 0.76148707238108, 0.659864747076619, 0.660451229917595, 0.630322020627465, 0.679375574025315, 0.684939970477629, 0.626276111252214, 0.650732224569865, 0.691289573122707, 0.649051786282571, 0.615812527367672, 0.621703633507002, 0.626700228233863, 0.751828202108001, 0.797753519987023, 0.75695689206503, 0.657415436891266, 0.66893069513576, 0.736935133715545, 0.627674339863497, 0.772694255427808, 0.653883277313985, 0.694208227682342, 0.697388865373198, 0.701924381142374, 0.636653681618197, 0.690857929827726, 0.716310315594732, 0.593933750238816, 0.594928872754016, 0.690287604255378, 0.750791185023354, 0.680515481470936, 0.665889901848514, 0.68140133801784, 0.662532365313363, 0.780446942599069, 0.615167699286637, 0.68518592047271, 0.762936034251056, 0.667501463300847, 0.671750693886387, 0.747830984224012, 0.690662847114914, 0.700209952499606, 0.758270095982651, 0.660436766804131, 0.745837642837692, 0.717507842121821, 0.692333307688362, 0.679680761934844, 0.783311978805604, 0.722241819228017, 0.669684615144087, 0.619384423671881, 0.708978899757948, 0.69885619660113, 0.698526395169173, 0.698908376898035, 0.646233816520023, 0.524284895554321, 0.639301885443374, 0.671639907716296, 0.666481333551904, 0.640899772897716, 0.697880916199803, 0.735296979117167, 0.624453167632793, 0.687693785247505, 0.670464857205835, 0.702839809119529, 0.609142738817578, 0.580027453365198, 0.63124339537534, 0.64105061655259, 0.745641826989972, 0.737839532186865, 0.609461707309304, 0.62555104675997, 0.644266345010154, 0.639762510271821, 0.701369207891536, 0.725318739884268, 0.599698564191191, 0.674392208484714, 0.686387124811052, 0.723609362978651, 0.735629372876929, 0.682255794168838, 0.730028893652169, 0.609284860502141, 0.636888137305202, 0.625368743182459, 0.680609108546504, 0.724767723226951, 0.682673007331043, 0.754500441360665, 0.718855177954899 ]

        n1_mean_maxtaxa_divtree = [ 40.3684210526316, 39.7373737373737, 39.6458333333333, 40.2421052631579, 40.0842105263158, 40.5222222222222, 40.0898876404494, 39.90625, 40.5157894736842, 40.2934782608696, 39.9157894736842, 40.4, 40.3085106382979, 40.6382978723404, 40.5824175824176, 39.5217391304348, 40.1333333333333, 40.563829787234, 40.6526315789474, 40.0645161290323, 40.8695652173913, 40.25, 39.7244897959184, 40.1666666666667, 40.3368421052632, 40.0537634408602, 40.21875, 40.2872340425532, 40.2159090909091, 39.6304347826087, 40.4301075268817, 39.969696969697, 39.554347826087, 40.3541666666667, 40.3404255319149, 40.1827956989247, 40.0674157303371, 40.5729166666667, 39.2555555555556, 40.5416666666667, 39.9895833333333, 40.436170212766, 40.0416666666667, 40.53125, 40.2553191489362, 40.5463917525773, 40.1063829787234, 40.3626373626374, 39.8191489361702, 40.9375, 39.9677419354839, 40.5670103092784, 40.4795918367347, 40.1666666666667, 39.4838709677419, 40.2417582417582, 39.9468085106383, 40.2340425531915, 40.2659574468085, 39.2526315789474, 39.8541666666667, 39.9479166666667, 40.7311827956989, 39.9318181818182, 40.6914893617021, 39.9684210526316, 40.2164948453608, 39.8877551020408, 40.2365591397849, 40.3052631578947, 39.4329896907216, 39.8526315789474, 40.3068181818182, 39.5217391304348, 40.2150537634409, 40.5957446808511, 40.3020833333333, 40.0212765957447, 40.2551020408163, 40.2736842105263, 39.9148936170213, 40.2340425531915, 40.2340425531915, 40.3870967741936, 40.1505376344086, 39.6086956521739, 39.8602150537634, 40.4, 39.8631578947368, 40.3085106382979, 40.0537634408602, 40.5918367346939, 40.0816326530612, 40.3157894736842, 40.1052631578947, 39.6489361702128, 40.0842105263158, 40.2210526315789, 40.1052631578947, 39.9893617021277 ]

        root_ages_ci_width_maxtaxa_divtree = [ 0.237212993809004, 0.296979779440502, 0.28806325958551, 0.266263651627938, 0.208575345317851, 0.219387757908826, 0.228462264833061, 0.262163861788498, 0.261355252477221, 0.201814174135395, 0.249221113245812, 0.227648626977026, 0.237461927675334, 0.201605514450131, 0.306443230702792, 0.253880308576643, 0.264693219550334, 0.247337895174861, 0.241695689593438, 0.252664409659066, 0.235370579335703, 0.277314442215125, 0.295148293136944, 0.232284025442549, 0.250796849515626, 0.345442823254911, 0.361594500566786, 0.310413070539827, 0.284136942428714, 0.316424933388719, 0.252354970707, 0.220257230805152, 0.230854563320041, 0.256373129296999, 0.234181478640652, 0.230416800243785, 0.211788777972263, 0.274927260561929, 0.265629656173639, 0.268635899575133, 0.276282070089439, 0.288650303478448, 0.216524899252375, 0.249700782678979, 0.366316701087045, 0.253311517999019, 0.236061913023948, 0.275602678394168, 0.192792682036902, 0.279252545937737, 0.176414003947713, 0.236075981338354, 0.302828249032044, 0.290686064311163, 0.214037457947085, 0.212662394571809, 0.225272421773522, 0.258677463358675, 0.262750401898894, 0.354915786639729, 0.266055800698691, 0.241419983384779, 0.211596214411715, 0.289079919013829, 0.260373512149469, 0.228523191146005, 0.232202052193019, 0.25417826422539, 0.259505274369138, 0.264670809292306, 0.258234211907364, 0.245114095879749, 0.253049897615512, 0.20679941023967, 0.285348557125798, 0.266690612886022, 0.221490727774162, 0.278958598628794, 0.287669065400626, 0.276015513185732, 0.231450002235123, 0.240590876451501, 0.231710272433262, 0.245502953956807, 0.22548404426651, 0.241267274481501, 0.219972841791308, 0.264374424927063, 0.245319021977643, 0.226684738116317, 0.21946245162533, 0.271629194029264, 0.275005878348308, 0.287177471717258, 0.272486464169989, 0.222265234389857, 0.196359245942406, 0.261702987054514, 0.213987492259741, 0.22563362332201 ]

        root_ages_mean_maxtaxa_divtree = [ 3.55172673053542, 3.9127941530696, 3.80430236524155, 3.73141417706893, 3.74556420199403, 3.66778891066675, 3.78692902376881, 3.82355427080851, 3.84183817250497, 3.39796869853732, 3.74235553252069, 3.62481722797747, 3.68348713831474, 3.4905364833486, 3.81810100696981, 3.82890319642297, 3.80880424382715, 3.74328937501739, 3.64858811295637, 3.78360095544487, 3.71369231703396, 3.9522425883831, 3.90438932184811, 3.67288412830961, 3.79884982112664, 3.84202517772462, 4.2265311063898, 3.92442063318461, 3.66633567518946, 3.90519524370196, 3.52471406271283, 3.62370284152017, 3.48618096580436, 3.76015121106648, 3.73003110012213, 3.76391373389574, 3.62033392230768, 3.84409305346056, 3.71020837256848, 3.73183352622779, 3.74564098325164, 3.9788601071399, 3.58401971485518, 3.71106697654576, 4.17593358262528, 3.72643511100609, 3.93983044683273, 3.88412385660439, 3.56812188399591, 3.61818361857227, 3.48258903923192, 3.78035893877615, 4.11804122568628, 3.68410297548811, 3.62766205435023, 3.75966381109593, 3.55533251807556, 3.77592628858866, 3.77846953611231, 3.93273223188244, 3.71778065650567, 3.82054183789191, 3.5763812100404, 3.82858401874954, 3.72426522210724, 3.71750643708242, 3.57482672841752, 3.92436176028331, 3.65568302579977, 3.76106590011963, 3.86024937603752, 3.74342937935341, 3.76319533044903, 3.6844556766734, 3.9379632770938, 3.73877504348952, 3.67310551279245, 3.91181291298287, 3.89328236529354, 4.11901284215373, 3.59271350915588, 3.86283995250313, 3.88521704143744, 3.93570383469811, 3.71959252602608, 3.8555927333007, 3.68110003385647, 4.01937281069923, 3.61964428541228, 3.4476883192311, 3.65231859420909, 3.61942791920685, 3.82680907358801, 4.01316244739759, 3.87737899433156, 3.8972061633736, 3.59952140342085, 3.57187052195236, 3.79747604161167, 3.70491914421258 ]


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

            # print("mean n0 = " + str(mean_n0) + ", diversitree mean n0 = " + str(n0_mean_maxtaxa_divtree[i]))
            # print("mean n1 = " + str(mean_n1) + ", diversitree mean n1 = " + str(n1_mean_maxtaxa_divtree[i]))
            # print("mean root age = " + str(mean_root_ages) + ", diversitree mean root age = " + str(root_ages_mean_maxtaxa_divtree[i]))

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
        stop_condition_value = [ 3.0 ] # 3.0 time units

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
            4.27, 3.51, 3.39, 4.51, 3.67, 2.93, 4.09, 3.32, 4.49,  3.8,
            3.9,  2.85, 3.49, 3.45, 3.65,  3.7,  4.5, 3.68, 3.89, 3.23,
            4.56, 4.66,  4.1, 3.12, 4.06, 3.38, 3.99, 3.11,  3.2, 4.07,
            3.6,  3.34, 4.19, 3.16, 3.33, 3.55,  3.4, 4.11, 3.63, 3.18,
            3.29,  3.5,  3.6, 3.51, 3.71, 4.36, 3.19, 3.21,  3.6, 3.94,
            3.71, 3.49, 3.68, 4.08, 3.49, 3.79, 4.03, 3.49, 3.41, 3.71,
            3.62, 3.59, 3.62, 2.83, 4.07, 3.79, 3.69, 3.38, 3.84, 4.51,
            3.51, 3.41, 2.91, 2.94, 3.82, 3.05, 3.18, 4.16, 4.69, 3.67,
            4.44, 3.16, 3.99, 4.13,  3.1, 2.92, 3.73, 3.34, 3.23, 3.43,
            4.03, 3.81, 3.04, 3.42, 3.54,  3.9,  4.5, 3.23, 4.19, 4.25 ]

        n0_ci_width_maxt_divtree = [
             0.88659874061984, 0.751395307061967,  1.00285554233293,  1.12317205856063, 0.867664144864145,  1.10567860309103, 0.941877959882877, 0.838472469687731,  1.00347444814545, 0.953305106596323,
            0.739949329414528, 0.536135178190099, 0.745172421509128, 0.844694792767952, 0.832665113907862, 0.910201220954061,  1.01118072798411, 0.819751700273327, 0.719419495499705, 0.705234192493973,
             1.12621684690943,  1.19570547638732,  1.00887560174581, 0.844191615952471, 0.697468975741009, 0.901324507756684, 0.872320670349609, 0.679475233046233, 0.876981231730806,  1.02443113635389,
            0.905713759496572, 0.840478597664041, 0.903867460561905, 0.615862170071561, 0.741811281508806,  0.77772615990333, 0.788933969404796, 0.873831806852249, 0.828591435674258,  0.78961243632301,
            0.818527141402239, 0.791144216762767,  0.88139485413116, 0.755515422701214, 0.967966156015756,  1.09271032064701, 0.638734530461631, 0.648142469685563, 0.844062901576761, 0.893558871574614,
            0.821838942170963, 0.872320670349609, 0.921379006893165, 0.921294772687774, 0.736793488603479, 0.967966156015756,  1.01231407351392, 0.811480551370472, 0.836718959128844, 0.834489921261407,
            0.792751358729631, 0.988906237142119, 0.766369645758099, 0.667467653360692, 0.916893350263582, 0.803213605060033, 0.885652862200102, 0.830520981436388, 0.853080591359067,  1.16222406017971,
            0.863378095481682, 0.790454791991179, 0.644419863132725, 0.704664475212326, 0.897268470797268, 0.602263038610845, 0.636717435386571,  1.02687706934937,  1.18935218017408,  0.91214105855787,
            0.937842151651407, 0.709548547893787,  1.05916048249393, 0.864348350089841, 0.658362461956278, 0.663751122885162, 0.695259528318147, 0.804625460520093, 0.720475904447115, 0.834936204905876,
             1.00539063497174,  1.12909878355718, 0.646118710673842, 0.767583889798886, 0.922868679934497, 0.946154962002109, 0.946154962002109, 0.647279773086155,  1.08278409499015, 0.963175786634816
             ]

        n1_mean_maxt_divtree = [
            15.71, 13.13, 13.72, 18.28, 13.95,  11.7, 16.25, 13.72, 16.06, 13.68,
            14.34, 11.64, 13.75, 13.48, 12.74, 12.86, 18.05, 13.74, 15.27, 12.29,
            17.55, 18.13, 15.49, 12.27, 13.92, 12.36, 15.72, 11.89, 11.62, 15.76,
            14.24, 13.23, 15.58, 12.16, 11.84, 14.14, 12.57, 15.02, 13.49, 13.03,
            12.34, 14.18, 13.36, 12.39, 12.99, 16.88, 14.94, 13.28, 12.42, 14.94,
            13.27, 13.52, 12.63, 15.63, 13.46, 14.49, 13.25, 12.39, 12.54,  13.1,
            14.48, 13.03, 14.61,  13.2, 15.75,  13.5, 15.51, 13.99, 13.43, 17.94,
            13.33, 14.86, 12.02, 11.16, 15.63, 12.52, 13.67, 17.42, 18.16, 13.88,
            16.76, 10.43, 15.96, 16.71, 11.54, 11.56, 15.12, 13.55, 13.69, 13.65,
            15.05, 13.29, 12.15, 13.36, 12.96, 15.01, 17.72, 12.72, 16.97, 14.18
            ]

        n1_ci_width_maxt_divtree = [
            3.41911306568818, 2.62247218733266, 3.49130165671513, 4.59348964258826, 3.48265443289883, 3.09796761008428, 4.06733690949662, 3.15328549287818, 4.01069659966264,  3.8110706098858,
            3.34672510971547, 2.28955308735089, 3.27584498467744, 3.52166776285316,  2.9810411641952, 2.93041800141726, 4.13677082378089, 2.55607695747316, 3.10414348595778, 2.57237126277203,
            4.37412976395525, 4.36136135377267, 4.07040013400514, 3.32189411004984, 3.02254444580534, 3.40081751407517, 3.47347294920932, 2.71769747433746, 3.24339399638058, 4.13289918032699,
            3.49186178212348, 3.33332212628419, 3.22412252217247, 2.67289779388694, 2.53689196435084, 3.19933513396048, 3.16129569421505, 3.54949014884781, 3.36098799783875, 3.71345619534518,
            2.78668894988916, 3.26302735458465, 3.86843556278689,  2.9461171512842, 3.52541376587289, 4.61583835579938, 2.99055438636012, 3.03502309671208, 2.96649961695156, 3.60425816943226,
            3.25105090602029, 3.63659036578013, 3.66149262073777, 3.89892332840851, 3.11805872888367, 4.24885865106945, 3.33629165088995, 2.93767545814267, 2.98153576696591, 3.12465751820269,
             3.0123594285657, 3.63978742658236, 3.15547702878533,  2.6752602139084, 3.54165343455639, 2.97085492575773, 3.68955844698758, 3.25875710309463, 3.19523722607713, 4.83376925332017,
            3.29811935661902, 3.73674476257157, 2.82839626786644, 2.77685462465442, 3.74842303941594, 2.51040807049452, 2.74485612906996, 4.82857736097384, 3.98194546513042, 3.49287955785145,
             3.6852043620813, 2.23999114392189, 3.96370001290754, 3.39176572923606, 3.31564633222605,  3.2183571601836, 3.18185346049129, 3.06812241468869,  3.4196123919993, 3.21462392666685,
            3.90028857259241, 4.40638384623991, 2.77456961809513, 3.22055080401041, 3.87368820864274, 3.22211376734233, 3.86618396485488, 2.65049911548978, 3.93062256687567, 3.30027468119142
            ]

        root_ages_mean_maxt_divtree = [
            2.16019074739235, 2.23985325286739, 2.10912011582493, 2.08439540715148, 2.13755192815195, 2.04205603439848, 2.08936733990997, 2.08124115321389, 2.08807639506063, 2.11766684110982,
             2.0757680302778, 2.00991995256452, 2.09897928747077, 2.04872196530228, 1.94931979163084, 1.97173449766771, 2.13918824968117, 2.03011397192554, 2.09541540835165, 2.06150834429441,
            2.17898464720682, 2.02624688600771, 2.09667511152114,  1.9754571847295, 2.07751262302644,   2.158776386014, 2.11504982259296, 2.10446610835592, 1.98656140946664, 2.13145930566523,
            2.05392749413261, 1.99715972917505,   2.094930002603,  2.0021278866453, 1.99800273456582, 2.06648552997892, 1.98789166119901, 2.09592859576434, 2.00193105859188, 2.00967553866122,
            2.12256784098586, 2.09696623504926,  1.9873263882158, 2.00317034938936, 2.05783703333234, 2.16276211673741, 2.01938809176279,  2.0778645042394, 2.07649292566626, 2.19972519294765,
            2.10734467222416, 1.95657366059353, 2.05713384191161, 2.11220638537172, 2.12465050375679, 2.22008365896472, 2.10599092168607, 1.92049420120249, 2.01865311783623, 1.98159012747726,
            2.03027863148591, 1.95820335112912, 1.96147032624551,  2.0424904697804, 2.02598968913367, 2.05281521660041, 1.87214418547339, 2.04883756427977,  2.1634103222118, 2.05967692385619,
            2.05446496840829, 2.00338268809876, 2.04269943590095, 1.96483549376368, 1.94699243666359, 2.00933205297941, 2.09417603053312, 1.92748936307839, 2.04613235181204, 2.02311834630362,
            2.17709563733529, 1.92738468601729, 2.16764177144626, 2.10799911987284,  2.0678850656439, 1.92466950356008, 2.07922029495231, 1.98672184902005,  2.0120457850184, 2.05198182016563,
            1.99223230590896, 1.98443208217947,  2.1411956180265, 2.07942520081614, 1.98297995887198, 2.09256649157022, 2.13358578217028, 2.05623500706743, 2.15917821269505, 2.04823393358692
            ]

        root_ages_ci_width_maxt_divtree = [
            0.131509336256971, 0.113368190652701, 0.132905496094447, 0.150075373915352, 0.124348500246804, 0.134352098826473, 0.142229926075057, 0.147769194139888, 0.142024033059407, 0.142960742226984,
            0.135376437999344, 0.127962850153195, 0.141706314053254, 0.138853892488107, 0.167679239140398, 0.141929760304129, 0.130405591631781, 0.140592563673068,  0.13479959633956, 0.140810658626993, 
            0.133925886541875, 0.140655099303131, 0.125554463378381, 0.163574574741462, 0.137954565306727, 0.133925102730878, 0.133117618975698, 0.123658101466697, 0.145921842031752, 0.136109701567449, 
             0.12919792470185, 0.146469249860248, 0.138961220099624, 0.141154937526522, 0.144905559356234, 0.136989581099801, 0.153565856318093,  0.14619780010837, 0.142149733645523, 0.125054194902086, 
            0.126578945088363, 0.130274514759025, 0.153630694718841, 0.137166225385444, 0.151670094693458, 0.131659161345938, 0.129825151242959, 0.153244384030126, 0.125697814332391, 0.122163352340789, 
             0.14169065050112,  0.14423459316781, 0.133861971841433, 0.129031343436785, 0.133393151436868, 0.124243690494388, 0.144132490964085, 0.149497756055504, 0.141231842519813, 0.154111281742736, 
            0.150148408014678, 0.155317703312488, 0.160034800733246, 0.143822107792082, 0.136440771932919, 0.144500713812307, 0.151579473914853, 0.139653180814706, 0.140066921781386, 0.141832052213146, 
            0.135947837835367, 0.137563663037994, 0.131445063292238, 0.157989987160941, 0.157077862816222,  0.13224158335642, 0.140424798436539, 0.146143851988532, 0.146891598902084, 0.150635164827912, 
            0.127420233505672, 0.157284650622664, 0.138457705182025, 0.140304209487322, 0.138073021199882, 0.149866940133464, 0.136411943841081,   0.1514282969283, 0.137968831355123, 0.141772791438584, 
            0.163207912189686, 0.148880499430412, 0.141555578522862, 0.143785731101053, 0.143693129538231, 0.126868118721876, 0.119099619176478, 0.138228325030658, 0.138096034395849, 0.148795113183592
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

            # debugging
            # print("PJ's vs. diversitree mean state 0 count: " + str(mean_n0) + " <-> " + str(n0_mean_maxt_divtree[i]))
            # print("PJ's vs. diversitree mean state 1 count: " + str(mean_n1) + " <-> " + str(n1_mean_maxt_divtree[i]))
            # print("PJ's vs. diversitree mean root age: " + str(mean_root_ages) + " <-> " + str(root_ages_mean_maxt_divtree[i]))

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
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_bisse.TestBiSSETrees.test_tree_size_state_count_max_taxa_bisse

    total_n_states = 2

    rates_t0_s0 = [ sseobj.MacroevolStateDependentRateParameter(name="lambda0", val=0.5, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
                    sseobj.MacroevolStateDependentRateParameter(name="mu0", val=0.25, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]),
                    sseobj.MacroevolStateDependentRateParameter(name="q01", val=0.8, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[0,1]) ]

    rates_t0_s1 = [ sseobj.MacroevolStateDependentRateParameter(name="lambda1", val=1.5, event=sseobj.MacroevolEvent.W_SPECIATION, states=[1,1,1]),
                    sseobj.MacroevolStateDependentRateParameter(name="mu1", val=0.25, event=sseobj.MacroevolEvent.EXTINCTION, states=[1]),
                    sseobj.MacroevolStateDependentRateParameter(name="q10", val=0.4, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[1,0]) ]

    rates_t0 = rates_t0_s0 + rates_t0_s1

    matrix_atomic_rate_params = [ rates_t0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

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

    # unittest.main()