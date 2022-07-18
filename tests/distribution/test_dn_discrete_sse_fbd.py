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
                        sseobj.MacroevolStateDependentRateParameter(name="mu", val=0.5, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]),
                        sseobj.MacroevolStateDependentRateParameter(name="psi", val=0.75, event=sseobj.MacroevolEvent.ANCESTOR_SAMPLING, states=[0]) ]

        matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

        fig_rates_manager = sseobj.FIGRatesManager(matrix_atomic_rate_params, total_n_states)

        cls.event_handler = sseobj.MacroevolEventHandler(fig_rates_manager)

    
    def test_tree_size_sa_count_max_taxa_fbd(self):
        """
        Test if FBD trees simulated here have similar root ages and number of sampled ancestors
        as trees simulated with FossilSim
        """

        stop_condition = "size"
        stop_condition_value = [ 10 ] # 10 taxa

        start_at_origin = True

        # simulation initialization
        n_batches = i = 100
        n_sim = 100
        start_states_list = [0 for i in range(n_sim)]

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

        n_sa_mean_maxtaxa_fossilsim = [
            11.79, 11.43, 12.31, 11.86, 13.68, 11.66, 11.75, 11.7, 11.82, 12.38, 11.98, 12.36, 11.95, 12.93, 11.88, 11.39, 11.44, 11.86, 13.44, 10.6, 12.51, 10.84, 11.73, 12.21, 11.81, 11.73, 12.26, 12.11, 12, 12.6, 11.77, 12.97, 12.35, 11.2, 11.58, 11.06, 11.66, 12.77, 12.48, 11.33, 11.79, 10.97, 13.07, 11.86, 12.33, 12.28, 13.29, 11.29, 12.66, 11.6, 12.35, 11.26, 11.77, 11.69, 11.92, 11.3, 12.24, 11.4, 12.84, 11.8, 11.34, 12.29, 12.21, 12.32, 11.44, 11.65, 13.44, 12.05, 12.22, 11.61, 11.89, 10.98, 11.78, 11.9, 12.22, 13.22, 9.65, 12.55, 11.21, 12.72, 11.86, 10.57, 12.34, 11.67, 12.22, 11.75, 12.85, 11.86, 11.55, 12.26, 11.82, 12.35, 11.59, 12.7, 11.93, 12.52, 11.47, 12.49, 11.33, 12.63
        ]

        n_sa_ci_width_maxtaxa_fossilsim = [
            1.33379389560741, 1.24593720418313, 1.56990316319771, 1.11572807519475, 1.12088433773778, 1.20861690315247, 1.23313962442756, 1.28887585503685, 1.20931019723283, 1.27989681149724, 1.2932279844629, 1.32948845877661, 1.17248886622183, 1.41481811716081, 1.28519152744726, 1.18057597926956, 1.32504457801893, 1.31732008581242, 1.39633926241354, 1.14048736702492, 1.39471956502033, 1.11144201101019, 1.30216129336553, 1.20928292233977, 1.30062271891956, 1.14213804186126, 1.51504195308115, 1.22572874485328, 1.4264935607243, 1.37016038285577, 1.37718152878708, 1.19350802242631, 1.3215418451509, 1.47464518588114, 1.19504325503897, 1.25170581886068, 1.34910349222262, 1.68878193914425, 1.42961298485452, 1.17568814046769, 1.19637868333542, 1.18698769304249, 1.53782837886036, 1.21750969813407, 1.30905657199636, 1.22356935449223, 1.30705118584707, 1.29960793404096, 1.25677968722534, 1.44218465033316, 1.30350754473916, 1.20327547850512, 1.25295608184417, 1.08421663868249, 1.4198625262804, 1.2987731273643, 1.26536470504265, 1.40318103254748, 1.44640815676613, 1.40179763377622, 1.23025538549333, 1.34162600518094, 1.2137669947629, 1.32900970255862, 1.17643547492698, 1.21347283580323, 1.59701279932633, 1.10009136443343, 1.27223370955071, 1.25760534285754, 1.44538162672309, 1.10838639362623, 1.2651990962419, 1.30800234763327, 1.27680062023731, 1.49751380007469, 0.910574177831402, 1.51625933924193, 1.3502750690212, 1.33721810796999, 1.00406398765872, 1.2727444940249, 1.42246736498397, 1.20083338805123, 1.23477643075149, 1.25311711545669, 1.33410806280163, 1.18717238389836, 1.43094101975595, 1.33650406265561, 1.14265604181617, 1.40466666666667, 1.25920880658684, 1.27403807246919, 1.31680007057206, 1.37707586317141, 1.31980242843428, 1.29847282418675, 1.15739282838279, 1.27894901992048
        ]

        root_ages_mean_maxtaxa_fossilsim = [
            3.06158475192417, 3.32709888885735, 3.1214987606465, 3.2038958126961, 3.03141896233894, 3.46974124531352, 3.50598379407897, 3.17238009650997, 3.2730016625159, 3.28443930284859, 3.42614433739572, 3.1799702838571, 3.19976981581782, 3.47522343843126, 3.183949056896, 3.16728712301645, 3.28761041840292, 3.11426523229198, 3.58531310345, 2.94317680635375, 3.29449840365967, 3.28878994083599, 3.23802550455375, 3.22190113550461, 3.36282308204269, 3.54382048378497, 3.20440286112108, 3.3090674701662, 3.34275889007505, 3.12691938172564, 3.33182999768509, 3.43457794037564, 3.27309332225342, 3.20495251252992, 3.43589418010246, 3.24701970302327, 3.28278784848371, 3.44945343469407, 3.19295764870201, 3.22819435797767, 3.28178198407458, 3.39723970249208, 3.60982901913729, 3.33668733806059, 3.44476782122808, 3.32025918256223, 3.51571270353605, 2.98238129728412, 3.46001154879427, 3.34857910267609, 3.58346899764418, 3.42406958205434, 3.50028638328715, 3.23462859507321, 3.26271504845709, 3.40853305622832, 2.97400800063136, 3.16754353212038, 3.64684379585319, 3.58176312946842, 3.1133520652582, 3.19872729959642, 3.12828002043412, 3.32382912239342, 3.26820188009603, 3.37489591781346, 3.62748409916923, 3.29571394856256, 3.30896200621846, 3.06386390733467, 3.32127490587135, 3.10097762840498, 3.20307109196218, 3.33146391677673, 3.26885893364074, 3.48889690295416, 2.64522072444859, 3.66406943536744, 3.05953473561372, 3.20892751670681, 3.51632147125951, 3.168446136035, 3.15701055536265, 3.23219303029446, 3.25490487429717, 3.31696905343935, 3.25640838754555, 3.04585612546581, 3.01576569359621, 3.90536299984387, 3.52072820067763, 3.2701692005532, 3.56062136501843, 3.36883225323135, 3.12855041294481, 3.53681433150419, 3.05565680557827, 3.54836816485227, 2.85258012676722, 3.53368111295637
        ]

        root_ages_ci_width_maxtaxa_fossilsim = [
            0.408647553505562, 0.509592298750689, 0.38409696476414, 0.409577052775503, 0.343181507630708, 0.425925871023221, 0.440410622940826, 0.373512788035978, 0.423466332276909, 0.446359844966612, 0.39550978128479, 0.362226738238393, 0.433430020093021, 0.441539771069205, 0.386224196961461, 0.416661593610685, 0.476251662617177, 0.394567481807831, 0.401871801250315, 0.377139877102328, 0.361005260062301, 0.47457693912755, 0.494768499892774, 0.364915372943971, 0.420463635827272, 0.470143983528569, 0.395008062847314, 0.450118684986919, 0.473402295787168, 0.374673879265125, 0.435951339487528, 0.422889389606616, 0.400228202075088, 0.427262159720166, 0.469923951198503, 0.418199242842914, 0.391065817814109, 0.473126296469968, 0.482313867304151, 0.467750693540001, 0.417134096093507, 0.460851228415337, 0.548406647618856, 0.402711017600549, 0.473654696487989, 0.449154692120125, 0.403322820470335, 0.454830338414632, 0.453998441357947, 0.485982655331107, 0.421788588327507, 0.441012302705562, 0.493085692361023, 0.3967052365066, 0.400791691564759, 0.462200455705074, 0.372325916558093, 0.437067833644792, 0.469982650720636, 0.549914384317812, 0.376250577728838, 0.399247826476163, 0.381805266058534, 0.399376308077217, 0.419076282996469, 0.406647571777729, 0.503686696621813, 0.402734834181398, 0.409564416715385, 0.414961924085472, 0.482648165564663, 0.379673919691292, 0.40496599740129, 0.438897720965336, 0.420990053026926, 0.501251156477642, 0.301854144258771, 0.483862166647989, 0.418853230258794, 0.395362614193838, 0.365611355058353, 0.532996409564814, 0.387043185904377, 0.403102021917312, 0.358345071423213, 0.360662145626226, 0.394681747797261, 0.346148061927539, 0.414274968052868, 0.525565802883034, 0.473154548027078, 0.456245337946746, 0.430655387821212, 0.362858201816026, 0.424946940315028, 0.444499780064382, 0.431266924028455, 0.399629131019165, 0.355152164986533, 0.473725519301307
        ]

        # parsing simulations
        n_sa_ci_overlap_count = 0
        root_age_ci_overlap_count = 0
        for i, batch in enumerate(sim_batches):
            n_sas = [ann_tr.n_sa for ann_tr in batch]
            root_ages = [ann_tr.root_age for ann_tr in batch]

            mean_sa = statistics.mean(n_sas)
            mean_root_ages = statistics.mean(root_ages)

            # debugging
            # print("PJ's vs. FossilSim mean sa count: " + str(mean_sa) + " <-> " + str(n_sa_mean_maxtaxa_fossilsim[i]))
            # print("PJ's vs. FossilSim mean root age: " + str(mean_root_ages) + " <-> " + str(root_ages_mean_maxtaxa_fossilsim[i]))

            stdevs_n_sa = statistics.stdev(n_sas)
            stdevs_root_ages = statistics.stdev(root_ages)

            sterr_n_sa = stdevs_n_sa / math.sqrt(n_sim)
            sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)

            n_sa_ci_width_maxtaxa = 1.96 * sterr_n_sa
            root_ages_ci_width_maxtaxa = 1.96 * sterr_root_ages

            if abs(mean_sa - n_sa_mean_maxtaxa_fossilsim[i]) <= (n_sa_ci_width_maxtaxa + n_sa_ci_width_maxtaxa_fossilsim[i]):
                n_sa_ci_overlap_count += 1

            if abs(mean_root_ages - root_ages_mean_maxtaxa_fossilsim[i]) <= (root_ages_ci_width_maxtaxa + root_ages_ci_width_maxtaxa_fossilsim[i]):
                root_age_ci_overlap_count += 1

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n95% CIs of simulations here and from FossilSim overlapped " + str(n_sa_ci_overlap_count) + " times for the sampled ancestor count.")
        print("\n95% CIs of simulations here and from FossilSim overlapped " + str(root_age_ci_overlap_count) + " times for root age.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(n_sa_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        self.assertAlmostEqual(root_age_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)

    def test_tree_size_sa_count_max_t_fbd(self):
        """
        Test if FBD trees simulated here have similar root ages and number of sampled ancestors
        as trees simulated with FossilSim

        Note: condition_on_speciation=False to match FossilSim!
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
            # print("Doing batch " + str(n_batches - i))
            sse_sim = distsse.DnSSE(self.event_handler, stop_condition_value, n=100, stop=stop_condition, origin=start_at_origin,
                                    start_states_list=start_states_list, epsilon=1e-12, runtime_limit=3600,
                                    condition_on_speciation=True, condition_on_survival=True,
                                    debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)


        # "expectations" from FossilSim
        n_sa_mean_maxt_fossilsim = [
            9.56, 8.73, 8.35, 9.66, 9.86, 9.6, 9.49, 8.24, 8.15, 8.41, 9.31, 9.07, 8.43, 8.8, 6.8, 9.44, 8.42, 9.07, 9.23, 10.48, 9.77, 7.74, 10.14, 9.47, 9.18, 9.19, 8.49, 9.78, 10, 8.84, 7.85, 9.2, 9, 8.37, 8.69, 9.59, 8.72, 9.04, 9.41, 8.7, 7.57, 7.82, 9.13, 9.79, 10.82, 9.63, 8.17, 8.86, 10.01, 8.02, 7.87, 9.29, 9.34, 9.22, 8.65, 9.6, 9.82, 8.42, 8.25, 10.01, 9.47, 9.4, 9.95, 7.54, 8.73, 9.75, 9.93, 9.48, 8.5, 10.7, 9.47, 7.48, 8.4, 9.49, 8.81, 8.12, 7.91, 7.54, 8.96, 8.33, 8.65, 9.75, 9.68, 8.55, 9.57, 8.53, 8.82, 10.95, 8.79, 9.13, 9.8, 8.11, 8.95, 9.14, 9.83, 8.81, 8.2, 9.5, 7.66, 9.24
        ]

        n_sa_ci_width_maxt_fossilsim = [
            1.57622429004417, 1.22349799634928, 1.28763334702222, 1.67618115459824, 1.52891167683276, 1.51385819225667, 1.59681232915626, 1.12131353277657, 1.21379257043385, 1.24245662095282, 1.32456860993466, 1.64605467759007, 1.29899868259122, 1.30339590642632, 0.848647741939487, 1.24595744791714, 1.19990399346587, 1.33610768996811, 1.24331831723811, 1.33733417686308, 1.50146668819444, 1.21099042794138, 1.3965337777237, 1.44499497940479, 1.50186075885546, 1.15492259673298, 1.25377342305683, 1.4591914827857, 1.67983029445888, 1.29494915870715, 1.44201647533983, 1.24429948650685, 1.38312660010542, 1.41025790663705, 1.54699596233305, 1.52027707150421, 1.25642765414528, 1.50921671248781, 1.33140613388002, 1.33006623971147, 1.16780643961405, 1.14469179902167, 1.27712728837304, 1.47918921771969, 1.57155465766234, 1.50296489794155, 1.17568814046769, 1.27023128666246, 1.34226216349697, 1.06736643410146, 1.10230431259964, 1.64818438656639, 1.29119801484333, 1.80403906696931, 1.29514541908067, 1.3808803516071, 1.40770352480108, 1.14156882022642, 1.33759384453177, 1.75230129737088, 1.44392041583763, 1.78966072062542, 1.46603505683827, 1.00344544546288, 1.31402711778054, 1.39746568595953, 1.62994555328289, 1.60480797605196, 1.25840887190487, 1.52751004890146, 1.56836003442314, 1.27585813213434, 1.18028176742041, 1.45277247949282, 1.23701036733282, 1.29301792842649, 1.13813557907819, 1.29275981285673, 1.6672910268099, 1.32438110478331, 1.46471102400384, 1.5246812947318, 1.64643771000585, 1.22080554627598, 1.46732616081845, 1.28191894731092, 1.33846532108889, 1.66832754819934, 1.51239514851527, 1.42067534431319, 1.68029223047599, 1.24239415596979, 1.56784532065346, 1.43922729212027, 1.48148547400809, 1.3460727997483, 1.20146173598259, 1.44420122272389, 1.06890677202117, 1.49142238338325 
        ]

        root_ages_mean_maxt_fossilsim = [
            2.34928018306087, 2.37694044114638, 2.2741160991596, 2.27493120089464, 2.2898058530331, 2.30978142364105, 2.11564313816687, 2.44469964652236, 2.27601601175793, 2.15770303824444, 2.28981788343307, 2.09549344753009, 2.08705626932286, 2.49560102538168, 2.14952356423281, 2.3738096689444, 2.32452342108283, 2.46653387020446, 2.22863764667063, 2.42924317846784, 2.25241621228355, 2.41661765185379, 2.30872857034547, 2.42706618690718, 2.43354168323478, 2.43912168366146, 2.34459168262921, 2.22489768062937, 2.32677644645318, 2.33095560290669, 2.23966946267758, 2.42862320073951, 2.28479653862915, 2.28491090777854, 2.135661652917, 2.40692053301125, 2.19788790262409, 2.15084166337265, 2.36896620909097, 2.33121486476817, 2.17868291089316, 2.28667812030766, 2.40913987068157, 2.30812645328225, 2.42957031905339, 2.29011334084868, 2.11056825231319, 2.30774731902097, 2.49938152091942, 2.36100042161558, 2.30322056656281, 2.17322194712571, 2.34453669106415, 2.31349880731835, 2.2245359230374, 2.37034574828528, 2.39534584153939, 2.34398715683377, 2.15089338709287, 2.25144000801653, 2.35413553229237, 2.28929758965049, 2.25607461322042, 2.24464202176275, 2.26316075839537, 2.34102193076088, 2.39589343682254, 2.3056034848673, 2.26397435253567, 2.47219279938485, 2.32405432178518, 2.15979124361155, 2.18697551131298, 2.23993302313543, 2.16437391678019, 2.24898977648034, 2.07822030198584, 2.16721341468534, 2.07585584668402, 2.21190130281405, 2.02626864846111, 2.19174445107652, 2.36535097063222, 2.34513416994933, 2.31534865060971, 2.11341394795276, 2.30155049176814, 2.45509551716975, 2.26071287765277, 2.42743601018436, 2.24543026059735, 2.31096663087081, 2.29946248994171, 2.24976224634868, 2.4072150797085, 2.34768682555717, 2.16864709288917, 2.44151438050047, 2.30999422221719, 2.29281418422062
        ]

        root_ages_ci_width_maxt_fossilsim = [
            0.160068528160401, 0.159525643483163, 0.186148448915979, 0.175320261381863, 0.164665691122548, 0.155224633297259, 0.201864167965336, 0.129999588032008, 0.181493928204013, 0.202786341627259, 0.171086876128891, 0.204136497051271, 0.210764038711101, 0.135366891678624, 0.19117166107702, 0.160851487108604, 0.160190146652771, 0.143669186367426, 0.181853594013861, 0.161071200982838, 0.170028123271311, 0.138759346975339, 0.170844966012847, 0.154928093331149, 0.138996173451172, 0.135564120948094, 0.170271632362141, 0.189905227020525, 0.167822373963756, 0.17255311112302, 0.178076860800695, 0.13634617749235, 0.167728389241753, 0.159963886274955, 0.205033964911946, 0.137180367225992, 0.193773250306032, 0.200451397977918, 0.161443558860366, 0.150547656345217, 0.187948577722585, 0.169360902130034, 0.165539313326966, 0.174756913225232, 0.154644464376885, 0.171024408329313, 0.207535609713255, 0.173841179264627, 0.128262750618909, 0.148896776858856, 0.177585604030253, 0.191542041401052, 0.17149455919184, 0.160027660931104, 0.176742999958746, 0.159545127147759, 0.134625713830695, 0.16082661360874, 0.191910083203907, 0.18746724410156, 0.168270113690349, 0.182463808064623, 0.187420789740848, 0.185539991994409, 0.183264494580736, 0.168933224726715, 0.153030398541846, 0.175170877307819, 0.173451931440737, 0.126919532548284, 0.161934720776669, 0.185577540281363, 0.183389757967752, 0.175668593801917, 0.196428306483782, 0.175722161795536, 0.187216233705957, 0.187608193036746, 0.210367030696603, 0.182559250749097, 0.214683518096769, 0.184837659631984, 0.152385150912183, 0.168830005584426, 0.169698527073066, 0.20029865372015, 0.166722291738935, 0.142396029704996, 0.185357837515451, 0.143839899685273, 0.180696150477176, 0.172362695453918, 0.179037078416514, 0.196060596745871, 0.158465912336148, 0.163406159340541, 0.195787087242981, 0.147432226118837, 0.167243144645565, 0.168850353377794
        ]

        # parsing simulations
        n_sa_ci_overlap_count = 0
        root_age_ci_overlap_count = 0
        for i, batch in enumerate(sim_batches):
            n_sas = [ann_tr.n_sa for ann_tr in batch]
            root_ages = [ann_tr.root_age for ann_tr in batch]

            mean_sa = statistics.mean(n_sas)
            mean_root_ages = statistics.mean(root_ages)

            # debugging
            # print("PJ's vs. FossilSim mean sa count: " + str(mean_sa) + " <-> " + str(n_sa_mean_maxt_fossilsim[i]))
            # print("PJ's vs. FossilSim mean root age: " + str(mean_root_ages) + " <-> " + str(root_ages_mean_maxt_fossilsim[i]))

            stdevs_n_sa = statistics.stdev(n_sas)
            stdevs_root_ages = statistics.stdev(root_ages)

            sterr_n_sa = stdevs_n_sa / math.sqrt(n_sim)
            sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)

            n_sa_ci_width_maxt = 1.96 * sterr_n_sa
            root_ages_ci_width_maxt = 1.96 * sterr_root_ages

            if abs(mean_sa - n_sa_mean_maxt_fossilsim[i]) <= (n_sa_ci_width_maxt + n_sa_ci_width_maxt_fossilsim[i]):
                n_sa_ci_overlap_count += 1

            if abs(mean_root_ages - root_ages_mean_maxt_fossilsim[i]) <= (root_ages_ci_width_maxt + root_ages_ci_width_maxt_fossilsim[i]):
                root_age_ci_overlap_count += 1

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n95% CIs of simulations here and from FossilSim overlapped " + str(n_sa_ci_overlap_count) + " times for the sampled ancestor count.")
        print("\n95% CIs of simulations here and from FossilSim overlapped " + str(root_age_ci_overlap_count) + " times for root age.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(n_sa_ci_overlap_count, exp_count,
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
    # $ python3 tests/distribution/test_dn_discrete_sse_fbd.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_discrete_sse_fbd
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_discrete_sse_fbd.TestFBDTrees.test_tree_size_sa_count_max_t_fbd

    #############
    # Debugging #
    #############
    # total_n_states = 1

    # # not state-dependent (just state 0, and no transition)
    # rates_t0_s0 = [ sseobj.MacroevolStateDependentRateParameter(name="lambda", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
    #                 sseobj.MacroevolStateDependentRateParameter(name="mu", val=0.5, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]),
    #                 sseobj.MacroevolStateDependentRateParameter(name="psi", val=0.75, event=sseobj.MacroevolEvent.ANCESTOR_SAMPLING, states=[0]) ]

    # matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

    # fig_rates_manager = sseobj.FIGRatesManager(matrix_atomic_rate_params, total_n_states)

    # event_handler = sseobj.MacroevolEventHandler(fig_rates_manager)

    # stop_condition = "size"
    # # stop_condition = "age"
        
    # stop_condition_value = [ 10 ] # 10 living taxa
    # # stop_condition_value = [ 3.0 ] # origin age 

    # start_at_origin = True

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
    #     condition_on_speciation=False,
    #     condition_on_survival=True,
    #     debug=False)

    # trs = sse_sim.generate()

    # # print(trs[0].tree.as_string(schema="newick"))

    # n_sa = float()
    # n_leaves = float()
    # root_age = float()
    # origin_age = float()
    # for tr in trs:
    #     # print(tr.tree.as_string(schema="newick"))
    #     # print(tr.n_sa)
    #     # print(str(tr.root_age) + "\n")
    #     n_sa += tr.n_sa
    #     n_leaves += tr.n_extant_terminal_nodes
    #     root_age += tr.root_age
    #     origin_age += tr.origin_age
    
    # print("mean sa = " + str(n_sa / n_sim))
    # print("mean leaves = " + str(n_leaves / n_sim))
    # print("mean root age = " + str(root_age / n_sim))
    # print("mean origin age = " + str(origin_age / n_sim))
    # print("mean leaves / mean sa = " + str(n_leaves / n_sa))
    # print("mean sa / mean root age = " + str(n_sa / root_age))

    unittest.main()