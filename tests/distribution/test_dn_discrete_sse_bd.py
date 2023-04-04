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
        rates_t0_s0 = [ sseobj.DiscreteStateDependentRate(name="lambda", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
                        sseobj.DiscreteStateDependentRate(name="mu", val=0.8, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]) ]

        matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

        state_dep_par_manager = sseobj.DiscreteStateDependentParameterManager(matrix_atomic_rate_params, total_n_states)

        event_handler = sseobj.MacroevolEventHandler(state_dep_par_manager)

        cls.sse_stash = sseobj.SSEStash(event_handler)


    def test_tree_size_total_count_max_taxa_bd(self):
        """
        Test if birth-death trees simulated here have similar total number of taxa and root ages
        as trees simulated with geiger, starting from the root

        Note that the method used by both PhyloJunction and geiger, when simulating for max taxa, is what
        Stadler (2011) calls the "simple sampling approach" (SSA).

        Trees simulated with TreeSim (which simulates using the GSA approach, see paper above) for the
        parameters below will have twice the number of tips, and larger root ages.
        """

        stop_condition = "size"
        stop_condition_value = [ 10 ] # 10 taxa

        start_at_origin = False

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]

        # simulations
        sim_batches = list()
        for i in range(n_batches):
            sse_sim = distsse.DnSSE(
                self.sse_stash,
                stop_condition_value,
                n=100,
                stop=stop_condition,
                origin=start_at_origin,
                start_states_list=start_states_list,
                epsilon=1e-12,
                runtime_limit=3600,
                condition_on_speciation=True,
                condition_on_survival=True,
                debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)

        n_total_mean_maxtaxa_geiger = [
            22.89, 23.03, 23.47, 22.71, 21.6, 23.14, 24.3, 24.36, 24.26, 26.05, 23.41, 24.5, 22.17, 23.03, 23.34, 23.19, 24.01, 23.67, 23.82, 21.03, 25.77, 22.82, 22.04, 23.03, 24.13, 22.77, 23.54, 23.74, 21.12, 23.83, 23, 22.83, 23.73, 22.94, 21.71, 21.7, 22.21, 25.24, 23.82, 24.15, 22.56, 23.78, 25.77, 25.72, 22.6, 24.23, 21.65, 25.78, 23.55, 24.99, 23.02, 23.01, 22.35, 23.19, 24.46, 22.37, 21.06, 22.26, 25.04, 24.45, 22.4, 22.98, 23.94, 23.56, 25.2, 22.06, 22.57, 23.64, 22.57, 23.47, 22.19, 24.1, 23.42, 25.42, 24.62, 22.93, 23.43, 22.76, 23.09, 22.9, 24.21, 23.12, 23.67, 24.73, 22.41, 22.7, 24.46, 23.38, 23.58, 21.7, 22.06, 21.82, 23.59, 23.86, 23.39, 23.97, 23.16, 23.75, 21.54, 25.02
        ]

        n_total_ci_width_maxtaxa_geiger = [
            2.16433413988943, 2.05008065512683, 2.53960324301242, 2.12149040768489, 1.94669016175112, 1.98166203236199, 2.29530489601198, 2.16587995190421, 2.26194735946036, 2.31581722382001, 2.00077448398303, 2.13629192572532, 1.87935798665933, 1.94737873693088, 2.12805128061551, 2.11826875718752, 2.31298706656587, 2.30501451113993, 2.13627012857041, 1.33558482086963, 2.48529970422039, 2.26121986422007, 1.953439849821, 2.13918806120736, 2.01595323458666, 1.91652757165904, 1.90425195005622, 2.23415685368523, 1.83975696945729, 2.38900417477046, 2.07389484747512, 1.90438134327294, 2.33142819921198, 2.00865736915649, 2.18101540952045, 1.856599824704, 2.084401426663, 2.64757530360948, 2.07040798243664, 2.05686837574783, 1.79254648232985, 2.40172494106325, 2.4116046616393, 2.42719828209523, 1.76762590310657, 2.06467899371082, 1.67413226128727, 2.7150811483306, 2.2131735155882, 2.40251891349412, 1.92966841596051, 1.83726960306413, 1.76529189505517, 2.09672627520411, 2.14588289249037, 2.24300337291523, 1.61116436420523, 2.15979786708541, 2.34570902774157, 2.29274646296409, 1.90722068662091, 1.83927179090206, 2.35394929324136, 2.26515310036227, 2.34158299214248, 1.79784228204565, 1.76746234843385, 2.18745096947072, 1.89584496276153, 2.31616572425684, 1.86227506740188, 2.23135528866257, 2.07186935979091, 2.42337433536235, 2.19210438297509, 2.01105817705882, 2.4510894385884, 2.53527008072557, 2.03919461978835, 2.06216749753234, 2.36192671188655, 2.35912966523117, 2.1225437039914, 2.40548248244159, 1.96752758516465, 2.04743754113489, 2.27423139961392, 2.42090716981228, 1.90888425224916, 1.3434194874245, 1.96016233671159, 2.23533065856739, 2.13730797126675, 2.60783658284569, 1.97044431394126, 2.25711044058772, 2.5028597603945, 2.00179825217834, 1.61995529755009, 2.16542483610216
        ]

        root_ages_mean_maxtaxa_geiger = [
            4.41540003385216, 4.71985308803584, 4.78379383074576, 4.45153696783428, 4.15493573002888, 4.4971376309519, 4.78609524630304, 4.85569472317847, 5.0639355631863, 5.28769069559235, 4.43338595924275, 4.99011701274351, 4.23555962232338, 4.52679664106864, 4.38418496772901, 4.58210611081857, 4.59798256048134, 4.76777716753414, 4.81256943543864, 4.09850103212462, 5.1476188129078, 4.36044023301376, 4.18169289420837, 4.37196755153014, 4.9122562415769, 4.48727782593609, 4.50827839153032, 4.58357231556283, 3.84098982778127, 4.88005050840536, 4.46224177624419, 4.4456937972613, 4.58177908352567, 4.39708876765071, 4.27278205928572, 4.15384372566421, 4.41753686559664, 5.36563825384304, 4.53917587273469, 4.57453608826879, 4.55620981249849, 4.45794365124995, 5.05272811015187, 5.01535850756605, 4.38910373238569, 4.862381656804, 4.10992932243869, 5.41385509206486, 4.56647109179781, 5.18159798645258, 4.66964068612041, 4.4541388448132, 4.40769619801235, 4.85158157115605, 5.19060502952889, 4.29951899167654, 3.97719811293962, 4.38312483410883, 5.16621405796241, 4.74725089637243, 4.18418324597874, 4.52221768601807, 4.9439341606343, 4.72384620343728, 5.11609079651273, 4.37979450012976, 4.25897648661528, 4.91755355338507, 4.35958543156291, 4.77674181275586, 4.2520448671119, 4.597822219599, 4.62740603262086, 5.3125759502359, 4.76164916282395, 4.08450520628267, 4.66405611675163, 4.52891752875248, 4.56063152895865, 4.56488191432983, 4.88138905204797, 4.36758713015436, 4.71921913628319, 4.88086094612638, 4.20295644209544, 4.43075410746347, 4.82556491591985, 4.69655259967226, 4.34264329420629, 4.13476454568959, 4.43580899043324, 4.00780146654569, 4.78997328319289, 4.59595115636419, 4.53576414993745, 4.46778315657736, 4.43009892547679, 4.41078796096448, 4.12916569383616, 4.93607105673977
        ]

        root_ages_ci_width_maxtaxa_geiger = [
            0.660375229352299, 0.598881818945529, 0.756839478185466, 0.611259204624358, 0.519048727326599, 0.578474277847004, 0.690616749550038, 0.569387891184332, 0.672662965200767, 0.632236562369346, 0.550946147735116, 0.631789795466899, 0.528984675670863, 0.575986080462836, 0.571070556065218, 0.613448835087499, 0.588891952157204, 0.61739594704027, 0.595035379814095, 0.452720914371002, 0.627198994353223, 0.598592897674005, 0.540861934038963, 0.582802793964172, 0.571700133184533, 0.561995890371707, 0.536844358548194, 0.605985747183232, 0.578833852211988, 0.645051344045209, 0.561503716284546, 0.548873314789491, 0.63257827206128, 0.518310279379211, 0.580700729754077, 0.57125559915378, 0.688543061490186, 0.780492656496002, 0.573593461594399, 0.562782897379988, 0.613782380741358, 0.605116465113965, 0.700136142705465, 0.695960552345689, 0.585087544753621, 0.52712255835094, 0.440768591214856, 0.839259844743783, 0.667130435428858, 0.724511105825431, 0.584789659709947, 0.504638070565292, 0.544683273517427, 0.610756357306458, 0.67476670574648, 0.564811884419875, 0.495773302291837, 0.620697898089063, 0.694447853989534, 0.584467286678396, 0.522019343033574, 0.479306433583065, 0.656327347649981, 0.628845133186212, 0.65420800918892, 0.562288132599757, 0.471001188757444, 0.672223362881446, 0.595934033992796, 0.655112577959579, 0.533079150771734, 0.591877980232354, 0.552705723616435, 0.68074031245756, 0.592410678608659, 0.488611338824823, 0.744861905313832, 0.760466849228129, 0.53836466360587, 0.609966786838998, 0.664359196652447, 0.581306538627766, 0.685454228706172, 0.669432498273439, 0.527543478889788, 0.606514744745762, 0.670207510533916, 0.691762283055467, 0.486598243073245, 0.394061161797239, 0.632166225960083, 0.580384916946076, 0.615571754635696, 0.68624668032033, 0.539826737325181, 0.576087281990491, 0.636250367953265, 0.502187046774736, 0.485283600695029, 0.636452298329218
        ]

        # parsing simulations
        n_total_ci_overlap_count = 0
        root_age_ci_overlap_count = 0
        global_mean_total = 0.0
        global_mean_root_age = 0.0
        for i, batch in enumerate(sim_batches):
            n_total = [len(ann_tr.tree.leaf_nodes()) for ann_tr in batch]
            root_ages = [ann_tr.root_age for ann_tr in batch]

            mean_total = statistics.mean(n_total)
            mean_root_ages = statistics.mean(root_ages)
            
            global_mean_total += mean_total
            global_mean_root_age += mean_root_ages

            # debugging
            # print("PJ's vs. geiger mean total taxon count: " + str(mean_total) + " <-> " + str(n_total_mean_maxtaxa_geiger[i]))
            # print("PJ's vs. geiger mean root age: " + str(mean_root_ages) + " <-> " + str(root_ages_mean_maxtaxa_geiger[i]))
            
            stdevs_n_total = statistics.stdev(n_total)
            stdevs_root_ages = statistics.stdev(root_ages)

            sterr_n_total = stdevs_n_total / math.sqrt(n_sim)
            sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)

            n_total_ci_width_maxtaxa = 1.96 * sterr_n_total
            root_ages_ci_width_maxtaxa = 1.96 * sterr_root_ages

            if abs(mean_total - n_total_mean_maxtaxa_geiger[i]) <= (n_total_ci_width_maxtaxa + n_total_ci_width_maxtaxa_geiger[i]):
                n_total_ci_overlap_count += 1

            if abs(mean_root_ages - root_ages_mean_maxtaxa_geiger[i]) <= (root_ages_ci_width_maxtaxa + root_ages_ci_width_maxtaxa_geiger[i]):
                root_age_ci_overlap_count += 1

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n\nPJ global mean total taxon count = " + str(global_mean_total / 100.0))
        print("geiger global mean total taxon count = " + str(statistics.mean(n_total_mean_maxtaxa_geiger)))
        print("\n\nPJ global mean root age = " + str(global_mean_root_age / 100.0))
        print("geiger global mean total root age = " + str(statistics.mean(root_ages_mean_maxtaxa_geiger)))

        print("\n95% CIs of simulations here and from geiger overlapped " + str(n_total_ci_overlap_count) + " times for the total taxon count.")
        print("\n95% CIs of simulations here and from geiger overlapped " + str(root_age_ci_overlap_count) + " times for root age.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(n_total_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_geiger) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        self.assertAlmostEqual(root_age_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_geiger) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)


    def test_tree_size_total_count_max_t_from_root_bd(self):
        """
        Test if birth-death trees simulated here have similar total number of taxa and
        number of extant taxa as trees simulated with geiger, starting from the root
        """

        stop_condition = "age"
        stop_condition_value = [ 3.0 ] # 3.0 time units

        start_at_origin = False

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]

        # simulations
        sim_batches = list()
        for i in range(n_batches):
            sse_sim = distsse.DnSSE(
                self.sse_stash,
                stop_condition_value,
                n=100,
                stop=stop_condition,
                origin=start_at_origin,
                start_states_list=start_states_list, epsilon=1e-12, runtime_limit=3600,
                condition_on_speciation=True, condition_on_survival=True,
                debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)


        # geiger from root
        n_total_mean_maxt_geiger = [
            17.22, 17.2, 16.46, 15.83, 15.63, 16.97, 16.6, 16.36, 16.22, 15.83, 16.31, 15.87, 15.82, 18.18, 17.28, 16.06, 17.76, 16.42, 16.13, 15.94, 18.57, 16.66, 16.21, 16.58, 16.32, 16.77, 18.09, 18.28, 18.35, 16.2, 16.8, 14.86, 17.2, 15.68, 16.79, 15.44, 17.71, 16.51, 14.83, 17.19, 16.41, 16.94, 15.69, 17.68, 16.09, 16.69, 17.85, 17.4, 16.13, 16.17, 17.22, 18.29, 16.82, 17.07, 15.39, 15.93, 17.19, 16.44, 18.1, 17.6, 17.66, 15.81, 16.43, 16.41, 16.64, 17.79, 16.32, 18.31, 18.22, 17.09, 16.53, 15.24, 17.3, 15.19, 15.57, 16.72, 17.71, 15.48, 17.85, 17.18, 17.6, 18.17, 16.28, 15.12, 16.02, 16.48, 17.59, 16.78, 16.58, 16.7, 19.04, 16.35, 16.67, 15.88, 15.33, 15.57, 17.73, 17.58, 16.3, 16.07
        ]

        n_total_ci_width_maxt_geiger = [
            2.02783782588273, 1.83631999148055, 1.78606217001299, 1.69341703968722, 1.83159209035597, 2.38435421780397, 2.00540538232916, 1.94978136210192, 2.32053846624729, 1.5375255529124, 1.84679139235288, 1.8844304258579, 1.84399157560687, 2.16907740048458, 2.05171255040635, 2.03933257610415, 1.78295264087299, 1.89521546910726, 1.9579508277892, 1.77656436751365, 1.87609284558428, 1.97404976511524, 2.0434147134914, 2.04813487403106, 1.66215358708836, 1.83671193740611, 2.1302155319416, 2.11136691460315, 2.31749222700362, 2.12655187302987, 1.91837827133211, 1.48226188283074, 1.9536543742795, 2.09194324661801, 1.72637749257475, 1.77120494921429, 2.03523271803636, 1.53814123637275, 1.70141826091818, 1.90695515469562, 1.78316482647674, 2.40014752543221, 1.81969780552551, 2.30834956161017, 1.61171219276992, 1.98551096708229, 2.43538509596667, 2.15176627471714, 2.31452329223601, 2.0486795002014, 1.91424335491028, 2.7547680637855, 2.55026510042079, 2.03636874317541, 1.58740445235407, 2.01683847041789, 2.22023165194149, 1.66085506224287, 2.09890988540973, 1.953852987036, 2.29745431152294, 1.97591139497906, 2.12641592556227, 1.90030990306233, 1.71175600262574, 2.27744966667596, 1.84767051853151, 1.96132502686579, 2.00532411184211, 1.98812808521793, 2.07135893167707, 1.78491031899888, 2.19258581717035, 1.7376690477752, 2.13642906105709, 1.98597702734229, 2.11305982705821, 1.62236325740874, 2.21492614282775, 2.26892902284471, 1.98967037521802, 2.074279316882, 1.8927364044065, 1.78124553222194, 1.72900195231559, 2.19098889015159, 1.81294675870501, 2.03242519849434, 1.89930600233481, 1.8763503366532, 2.91247337255786, 1.67644851885277, 1.8191347531453, 1.92509822118564, 1.48410242833093, 1.8669699541495, 1.8190920905604, 2.16241220936997, 1.93278684021207, 1.78276437301652
        ]

        n_extant_mean_maxt_geiger = [
            7.64, 7.41, 7.18, 6.96, 6.55, 7.68, 6.82, 6.89, 6.47, 6.79, 6.66, 6.61, 7.22, 8.04, 7.25, 6.88, 7.31, 7.17, 6.5, 6.67, 7.83, 7.17, 7.17, 7.27, 6.15, 6.95, 7.94, 7.64, 7.88, 6.64, 6.71, 6.34, 7.47, 6.62, 7.29, 6.69, 7.74, 7.04, 6.11, 7.14, 7.05, 7.76, 6.85, 7.48, 6.58, 6.89, 7.97, 7.8, 7.6, 6.94, 7.52, 8.02, 7.11, 7.3, 6.54, 7.01, 7.35, 6.98, 7.72, 7.73, 7.39, 6.88, 7.03, 6.35, 6.7, 7.47, 7.23, 8.07, 7.85, 7.43, 7.06, 6.45, 7.58, 6.73, 6.95, 7.27, 7.35, 6.12, 8.01, 7.5, 7.62, 7.63, 7.23, 6.45, 6.77, 6.68, 7.4, 7.51, 6.98, 6.97, 8.62, 6.99, 7.01, 7.06, 6.74, 6.56, 7.87, 7.65, 6.97, 6.86
        ]

        n_extant_ci_width_maxt_geiger = [
            1.02642350839498, 1.08435978896268, 1.03352895587219, 0.94345248974394, 0.970400458833904, 1.26712985215495, 1.06021692158165, 1.00440209180398, 1.06389799085697, 0.907124357739357, 1.00761323668731, 0.890110676582219, 1.05735826024345, 1.20594914813137, 0.987052402076224, 1.12447896333809, 0.891765733607411, 1.0490592310747, 0.987101542165473, 0.901012324941768, 0.945972439430249, 1.05201421889995, 1.13504926142374, 1.07961143103443, 0.815239965912002, 0.945693460489244, 1.25789072384342, 1.09093328888988, 1.2530876973994, 1.07082907859077, 0.955051775210817, 0.828855990464745, 1.08520399125061, 1.0400413014332, 0.945660634026538, 1.00006572632487, 1.09735249734729, 0.90224535532865, 0.773483558031749, 0.915099291626517, 1.02976570837323, 1.26014065208226, 0.950196300824735, 1.20677260274978, 0.8618868997025, 0.991570923007293, 1.28131339876651, 1.2311319255798, 1.357927833134, 1.1205173148726, 0.917411638659263, 1.51180106435398, 1.36246240478203, 1.13588482841983, 0.841032442953903, 1.13280436782243, 1.13413267711187, 0.941006206357875, 1.07573164465546, 1.06039075778315, 1.15733918375843, 1.02330362685142, 1.12350367632112, 0.940344123556409, 0.895155536397689, 1.15147693722875, 0.960949337842006, 1.13784915014463, 1.22872625186418, 1.10848266553789, 1.17164297247476, 0.890755643350979, 1.11299805221861, 0.977762253228454, 1.15649395269766, 1.12122182344697, 1.06859997145326, 0.734041813137606, 1.21478960160252, 1.21319299887745, 1.09811604030712, 1.16594417871474, 1.07565047907362, 0.957113568934421, 0.965381022595679, 1.06259326464109, 0.901849592154194, 1.06464178942872, 1.00211428691562, 1.04698577725668, 1.59517968640501, 0.945742698046932, 0.99493074648592, 1.08426137517527, 0.82152017158633, 1.03986966083022, 1.07304993887197, 1.08517538497775, 1.10855267622229, 0.954126988138869 , 0.985872305128916, 1.15826420530591, 0.987083852015162, 0.991758747953056, 1.01179262377987, 1.39449697038181, 0.955102561735829, 1.19947542136746, 0.966040005027292, 1.19562758703063, 1.15156623708054, 0.933969649756707, 1.02045563679499, 1.00651700673584, 1.03775354190548, 1.06965978467431, 0.894721942021852, 0.879356307279089, 1.11804196520506, 1.39607662375456, 1.17744271772998, 0.998934022553112, 1.05591684015993, 1.22293967342976, 0.901718349720256, 1.03315155781431, 1.09962565604589, 0.877819317365136, 1.04677820598504, 1.00685621414361, 1.13540475157496, 0.984937058954028, 1.24813562070665, 0.941548313506358, 1.0324452089554, 1.27123899731043, 1.21677642941993, 0.96923211592758, 1.45421412231327, 1.19602996077469, 1.22509542254611, 1.09341322694052, 1.25210876605038, 0.951235058099509, 1.28453163668964, 1.40470119760234, 1.06347481456409, 1.08411821190022, 1.09052594113246, 0.998871867901408, 0.707563317936203, 0.980346403424356, 1.07183419417276, 1.14289373367109, 1.11183820565564, 1.22686948257454, 0.908560527383395, 0.89644640394812, 0.979358335386445, 0.830574710606451, 1.24996856138242, 0.942999953083324, 1.09021454658051, 1.03954682511956, 0.997019048681791, 0.965879319061082, 1.07834011163623, 1.10588915439938, 1.27430454790371, 1.10131819712249, 1.03091249880968, 0.745372879094007, 1.03874956995943, 1.06277401457723, 0.983367789021406, 1.096883857281, 1.33250154115172, 1.26632572827887, 1.31497175819425, 1.07731765365127, 1.01481980288287, 1.19296331219465, 1.0298203463996, 1.08787000411749, 1.0022846499852, 0.985165537851227, 0.931889962724281, 1.08440451954988, 1.14676454832408, 1.03803955394559, 0.844235282442763, 0.879464413691364, 1.087406198625, 0.939931375100858, 0.924110339294999, 1.17629528279392, 1.28923558214262, 1.0254249681242, 0.807914998029127
        ]

        # parsing simulations
        n_total_ci_overlap_count = 0
        n_extant_ci_overlap_count = 0
        # root_age_ci_overlap_count = 0
        global_mean_total = 0.0
        global_mean_extant = 0.0
        for i, batch in enumerate(sim_batches):
            n_total = [len(ann_tr.tree.leaf_nodes()) for ann_tr in batch]
            n_extant = [ann_tr.n_extant_terminal_nodes for ann_tr in batch]
            # root_ages = [ann_tr.root_age for ann_tr in batch]

            mean_total = statistics.mean(n_total)
            mean_extant = statistics.mean(n_extant)
            # mean_root_ages = statistics.mean(root_ages)

            # debugging
            # print("PJ's vs. geiger mean total taxon count: " + str(mean_total) + " <-> " + str(n_total_mean_maxt_geiger[i]))
            # print("PJ's vs. geiger mean root age: " + str(mean_root_ages) + " <-> " + str(root_ages_mean_maxt_geiger[i]))

            global_mean_total += mean_total
            global_mean_extant += mean_extant

            stdevs_n_total = statistics.stdev(n_total)
            stdevs_n_extant = statistics.stdev(n_extant)

            sterr_n_total = stdevs_n_total / math.sqrt(n_sim)
            sterr_n_extant = stdevs_n_extant / math.sqrt(n_sim)

            n_total_ci_width_maxt = 1.96 * sterr_n_total
            n_extant_ci_width_maxt = 1.96 * sterr_n_extant

            if abs(mean_total - n_total_mean_maxt_geiger[i]) <= (n_total_ci_width_maxt + n_total_ci_width_maxt_geiger[i]):
                n_total_ci_overlap_count += 1

            if abs(mean_extant - n_extant_mean_maxt_geiger[i]) <= (n_extant_ci_width_maxt + n_extant_ci_width_maxt_geiger[i]):
                n_extant_ci_overlap_count += 1

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n\nPJ global mean total taxon count = " + str(global_mean_total / 100.0))
        print("geiger global mean total taxon count = " + str(statistics.mean(n_total_mean_maxt_geiger)))
        print("\n\nPJ global mean extant taxon count = " + str(global_mean_extant / 100.0))
        print("geiger global mean extant taxon count = " + str(statistics.mean(n_extant_mean_maxt_geiger)))

        print("\n\n95% CIs of simulations here and from geiger overlapped " + str(n_total_ci_overlap_count) + " times for the sampled ancestor count.")
        print("\n95% CIs of simulations here and from geiger overlapped " + str(n_extant_ci_overlap_count) + " times for the sampled ancestor count.")

        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(n_total_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_treesim) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        self.assertAlmostEqual(n_extant_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_treesim) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)


    def test_tree_size_total_count_max_t_from_origin_bd(self):
        """
        Test if birth-death trees simulated here have similar total number of taxa and
        number of extant taxa as trees simulated with TreeSim, starting from the origin
        """

        stop_condition = "age"
        stop_condition_value = [ 3.0 ] # 3.0 time units

        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]

        # simulations
        sim_batches = list()
        for i in range(n_batches):
            sse_sim = distsse.DnSSE(
                self.sse_stash,
                stop_condition_value,
                n=100,
                stop=stop_condition,
                origin=start_at_origin,
                start_states_list=start_states_list,
                epsilon=1e-12,
                runtime_limit=3600,
                condition_on_speciation=True, condition_on_survival=True,
                debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)


        # treesim from origin
        n_total_mean_maxt_treesim = [
            12.86, 11.36, 11.82, 12.37, 10.83, 11.74, 12.82, 11.53, 11.21, 10.57, 12.51, 10.12, 9.89, 12.38, 10.78, 12.45, 10.51, 11.57, 10.66, 12.68, 12.31, 11.16, 11.01, 12.15, 11.18, 10.53, 12.43, 12.59, 11.27, 12.36, 10.71, 12, 12.07, 12.22, 11.3, 9.85, 11.14, 11.18, 12.53, 11.29, 11.5, 11.02, 11.6, 10.7, 11.61, 10.94, 10.93, 12.57, 11.72, 11.35, 12.33, 12.31, 11.45, 13.19, 12.18, 12.7, 10.52, 13.22, 11.35, 12.27, 10.53, 11.73, 10.48, 11.39, 11.83, 11.68, 11.78, 12.35, 11.3, 10.47, 11.54, 10.5, 10.9, 11.3, 13.14, 11.88, 12.34, 12, 11.34, 12.46, 10.53, 13.04, 10.54, 11.33, 10.7, 13.26, 12.52, 9.68, 13.07, 10.53, 10.63, 9.87, 11.21, 11.07, 11.19, 10.61, 10.45, 12.46, 11.6, 12.7
        ]

        n_total_ci_width_maxt_treesim = [
            2.00598965887543, 1.53334677411764, 1.84672520462957, 1.81520559806109, 1.38622464458418, 1.74856040019535, 1.8082929395293, 1.58998374478699, 1.40740440743256, 1.46520900070863, 1.60626173647404, 1.56311217974419, 1.4053626464656, 1.72940587964933, 1.6328771241665, 1.79428933626164, 1.43692696734905, 1.55613946091544, 1.56483407659494, 1.68602813449332, 1.7642628535923, 1.35498704395759, 1.59875522235585, 1.82325125293064, 1.52849538534806, 1.45276179531931, 1.82641540010409, 1.99630876298525, 1.49135343371128, 1.8181671330316, 1.51136850636111, 1.91918720169543, 1.65451954074501, 1.75254154999088, 1.52165603567596, 1.26421573085507, 1.48644461121907, 1.71376331168713, 1.84253373681817, 1.26756310161098, 1.78586227935081, 1.6572671656018, 1.34038318055214, 1.51501121774691, 1.41279805701275, 1.34264515886598, 1.38488034767403, 1.81148214376287, 1.78468420815232, 1.55940764459464, 2.02658883861727, 1.73565808623796, 1.55666801635137, 1.693105371914, 1.57451484963374, 1.86119227563105, 1.57108544913536, 2.25110677730495, 1.34049173850871, 1.66191427704392, 1.30293585568914, 1.3732311647786, 1.32743208983786, 1.65724726317865, 1.62202238697201, 1.45794108373293, 1.65201399754385, 1.80011899718891, 1.61493639365628, 1.40662116313016, 1.24505394568941, 1.17451421292534, 1.59608435006589, 1.53359476032259, 2.04032178250432, 1.92469504067423, 2.0369336606227, 1.66567970777364, 1.66129891690602, 2.0020008109061, 1.55718642468143, 1.71998777726453, 1.43393305790244, 1.42378571251137, 1.77299614019075, 1.92677857449981, 1.72754253952213, 1.42618886136543, 1.8242895618998, 1.47581701583984, 1.44262988302403, 1.08777547519902, 1.79556916341478, 1.47681582173036, 1.32544718674434, 1.56276831761636, 1.26175780864505, 1.63734818076379, 1.96810043268409, 1.74919941283394
        ]

        root_ages_mean_maxt_treesim = [
            2.41134926469317, 2.38920913840397, 2.43975399368484, 2.36167155882785, 2.3587007364241, 2.33778283109888, 2.4577397815274, 2.42626338160149, 2.39713681308128, 2.4084376193809, 2.4974528389355, 2.40511012325899, 2.30647267157316, 2.42565744285912, 2.26765674607618, 2.4834094327036, 2.24709301128232, 2.45818899585454, 2.40231025551501, 2.42716561100311, 2.51587695830744, 2.43829360327129, 2.36654427609215, 2.29366342006442, 2.44129237874595, 2.44364909306366, 2.45982058094011, 2.3500470152416, 2.418030790723, 2.36690045733178, 2.33344159000591, 2.3488046990303, 2.43983650212711, 2.45543345103501, 2.42322554561376, 2.37134138005089, 2.39025548763228, 2.38353641595131, 2.41053305434273, 2.45320698395875, 2.56030978698553, 2.38549730708735, 2.40636439184482, 2.33156546281527, 2.45372703782439, 2.43273393188503, 2.43442793638736, 2.45165518183406, 2.35961928551049, 2.40281359050419, 2.3125562566843, 2.39925423527952, 2.39144320514972, 2.4277004420588, 2.44469274408377, 2.41753460371012, 2.31053533778282, 2.34655967858245, 2.36056707080135, 2.4257683669853, 2.36392051223435, 2.4433124600778, 2.3708687917585, 2.37442558977751, 2.40860605171841, 2.4345556612665, 2.35992545428336, 2.44112488756847, 2.38450054327198, 2.44956570399421, 2.46092076374808, 2.43811520578555, 2.29448856565299, 2.37649818075683, 2.38469577201613, 2.4643347062554, 2.44166991306082, 2.45571811529939, 2.48262532651773, 2.4040188589081, 2.49845301372491, 2.52427137935501, 2.30222139001482, 2.31649414607209, 2.40718261057246, 2.47965263780969, 2.52176339091588, 2.34595828935434, 2.41761244731798, 2.32563882473638, 2.35136820012941, 2.43721756305005, 2.37140250618802, 2.34120221673433, 2.42689883102106, 2.30395083398858, 2.32984965307698, 2.48754373987955, 2.35823136681131, 2.43976798550804
        ]

        root_ages_ci_width_maxt_treesim = [
            0.118614537742372, 0.0973025403085302, 0.0924813356752944, 0.110194862975739, 0.114440741277904, 0.110178798503782, 0.0957698644204067, 0.116719073615701, 0.110946947854556, 0.106754991316681, 0.0871260862329165, 0.117324236823729, 0.12459877289986, 0.119319780757935, 0.135223244008495, 0.0845561395512595, 0.125747814556681, 0.109656611487338, 0.122301468085595, 0.0910437169307571, 0.110916473580753, 0.117261807568291, 0.105129241091335, 0.130308337844548, 0.117198030247258, 0.111154148813296, 0.0899255607893612, 0.12119151921046, 0.116914638975389, 0.121049494696862, 0.120881228390186, 0.127553739150464, 0.112069672234155, 0.0981225418422778, 0.10700017134518, 0.124603154132274, 0.13080762953373, 0.111977296852409, 0.105382701132416, 0.101228706099052, 0.0880208207778674, 0.122068533115286, 0.0958017554961159, 0.122574672405576, 0.0934483208090163, 0.104521531144681, 0.105018649124541, 0.110618124060909, 0.129357799956822, 0.0991629983683035, 0.113002372302917, 0.120914396167643, 0.0998846417472681, 0.112128717510314, 0.0962641661378887, 0.107585124592922, 0.1271274597584, 0.120671574088347, 0.117222788974941, 0.113849859692386, 0.0988910493234009, 0.109599399979178, 0.109263277034012, 0.12104537332385, 0.118114588466149, 0.0992280983891192, 0.111306108533451, 0.101338565822798, 0.0992505816563806, 0.107205384355046, 0.0846739513590895, 0.101923577275351, 0.120585421036524, 0.13533080458462, 0.119862388629098, 0.107469642266994, 0.10187092799781, 0.102821551426194, 0.0848564790648721, 0.109156482672756, 0.110050492222876, 0.0881731057795677, 0.113992856348529, 0.121231280141388, 0.104053085645244, 0.10507398152035, 0.0917024098969028, 0.119101188258717, 0.104633851513134, 0.123410247752883, 0.123102428889098, 0.091508461525494, 0.115638713124, 0.117966708499831, 0.103408689657657, 0.111697405092111, 0.118632053538949, 0.0972386360773996, 0.113507182558513, 0.0976551271056787
        ]

        # parsing simulations
        n_total_ci_overlap_count = 0
        root_age_ci_overlap_count = 0
        global_mean_total = 0.0
        global_mean_root_age = 0.0
        for i, batch in enumerate(sim_batches):
            n_total = [len(ann_tr.tree.leaf_nodes()) for ann_tr in batch]
            root_ages = [ann_tr.root_age for ann_tr in batch]

            mean_total = statistics.mean(n_total)
            mean_root_ages = statistics.mean(root_ages)

            # debugging
            # print("PJ's vs. TreeSim mean total taxon count: " + str(mean_total) + " <-> " + str(n_total_mean_maxt_treesim[i]))
            # print("PJ's vs. TreeSim mean root age: " + str(mean_root_ages) + " <-> " + str(root_ages_mean_maxt_treesim[i]))

            global_mean_total += mean_total
            global_mean_root_age += mean_root_ages

            stdevs_n_total = statistics.stdev(n_total)
            stdevs_root_ages = statistics.stdev(root_ages)

            sterr_n_total = stdevs_n_total / math.sqrt(n_sim)
            sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)

            n_total_ci_width_maxt = 1.96 * sterr_n_total
            root_ages_ci_width_maxt = 1.96 * sterr_root_ages

            if abs(mean_total - n_total_mean_maxt_treesim[i]) <= (n_total_ci_width_maxt + n_total_ci_width_maxt_treesim[i]):
                n_total_ci_overlap_count += 1

            if abs(mean_root_ages - root_ages_mean_maxt_treesim[i]) <= (root_ages_ci_width_maxt + root_ages_ci_width_maxt_treesim[i]):
                root_age_ci_overlap_count += 1

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n\nPJ global mean total taxon count = " + str(global_mean_total / 100.0))
        print("TreeSim global mean total taxon count = " + str(statistics.mean(n_total_mean_maxt_treesim)))
        print("\n\nPJ global mean root age = " + str(global_mean_root_age / 100.0))
        print("TreeSim global mean root age = " + str(statistics.mean(root_ages_mean_maxt_treesim)))

        print("\n\n95% CIs of simulations here and from TreeSim overlapped " + str(n_total_ci_overlap_count) + " times for the sampled ancestor count.")
        print("\n95% CIs of simulations here and from TreeSim overlapped " + str(root_age_ci_overlap_count) + " times for the sampled ancestor count.")

        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(n_total_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_treesim) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        self.assertAlmostEqual(root_age_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_treesim) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)


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
            sse_sim = distsse.DnSSE(
                self.sse_stash,
                stop_condition_value,
                n=n_sim,
                stop=stop_condition,
                origin=start_at_origin,
                start_states_list=start_states_list,
                epsilon=1e-12,
                runtime_limit=3000,
                condition_on_speciation=False,
                condition_on_survival=False,
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
            n_extant_obs_nodes = [ann_tr.n_extant_terminal_nodes for ann_tr in batch]
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

    #############
    # Debugging #
    #############
    # total_n_states = 1

    # # not state-dependent (just state 0, and no transition)
    # rates_t0_s0 = [ sseobj.DiscreteStateDependentRate(name="lambda", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
    #                 sseobj.DiscreteStateDependentRate(name="mu", val=0.8, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]) ]

    # matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

    # state_dep_par_manager = sseobj.DiscreteStateDependentParameterManager(matrix_atomic_rate_params, total_n_states)

    # event_handler = sseobj.MacroevolEventHandler(state_dep_par_manager)

    # # stop_condition = "size"
    # stop_condition = "age"
        
    # # stop_condition_value = [ 10 ] # 10 living taxa
    # stop_condition_value = [ 3.0 ] # origin age 

    # start_at_origin = False

    # # simulation initialization
    # n_sim = 100
    # start_states_list = [0 for i in range(n_sim)]

    # sse_sim = distsse.DnSSE(
    #     event_handler,
    #     stop_condition_value,
    #     n=n_sim,
    #     stop=stop_condition,
    #     origin=start_at_origin,
    #     start_states_list=start_states_list,
    #     epsilon=1e-12,
    #     runtime_limit=3600,
    #     condition_on_speciation=True,
    #     condition_on_survival=True,
    #     debug=False)

    # trs = sse_sim.generate()

    # print(trs[0].tree.as_string(schema="newick"))

    # n_leaves = float()
    # root_age = float()
    # origin_age = float()
    # for tr in trs:
        # print(tr.tree.as_string(schema="newick"))
        # print("root age = " + str(tr.root_age) + "\n")
        # print("total taxa = " + str(tr.n_extant_terminal_nodes + tr.n_extinct_terminal_nodes) + "\n")
        # n_leaves += tr.n_extant_terminal_nodes + tr.n_extinct_terminal_nodes
        # n_leaves += len(tr.tree.leaf_nodes())
        # root_age += tr.root_age
        # origin_age += tr.origin_age

    # print("mean extant and fossil leaves = " + str(n_leaves / n_sim))
    # print("mean root age = " + str(root_age / n_sim))

    unittest.main()