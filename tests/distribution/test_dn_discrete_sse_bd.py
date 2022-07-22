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
        rates_t0_s0 = [ sseobj.MacroevolStateDependentRateParameter(name="lambda", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
                        sseobj.MacroevolStateDependentRateParameter(name="mu", val=0.8, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]) ]

        # original implementation
        # matrix_atomic_rate_params = [ [rates_t0_s0] ] # 1D: time slices, 2D: states, 3D: parameters of state, several parameters -> matrix
        matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

        fig_rates_manager = sseobj.FIGRatesManager(matrix_atomic_rate_params, total_n_states)

        cls.event_handler = sseobj.MacroevolEventHandler(fig_rates_manager)


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
            sse_sim = distsse.DnSSE(self.event_handler, stop_condition_value, n=100, stop=stop_condition, origin=start_at_origin,
                                    start_states_list=start_states_list, epsilon=1e-12, runtime_limit=3600,
                                    condition_on_speciation=True, condition_on_survival=True,
                                    debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)

        n_total_mean_maxtaxa_geiger = [
            25.23, 21.84, 23.09, 22.63, 21.5, 23.11, 22.74, 23.17, 22.69, 22.49, 23.91, 25.15, 23.71, 23.52, 25.05, 22.68, 24.09, 23.96, 23.51, 23.25, 22.21, 23.49, 23.22, 23.88, 23.25, 24.24, 21.91, 22.79, 24.26, 25.06, 24.78, 22.59, 24.41, 21.74, 22.28, 23.47, 25.21, 24.54, 21.61, 23.04, 23.17, 25.48, 22.42, 25.99, 22.42, 23.5, 22.75, 23.1, 23.84, 22.52, 21.77, 24.41, 24.09, 22.27, 22.88, 24.28, 22.74, 23.83, 24.16, 23.4, 22.45, 23.32, 23.92, 23.93, 23.99, 23.35, 22.71, 22.24, 22.96, 23.28, 24.04, 23.19, 23.9, 24.15, 22.57, 23.63, 24.06, 22.73, 24.48, 22.45, 23.73, 22.22, 23.91, 24.71, 21.92, 24.48, 22.64, 20.97, 23.83, 23.22, 22.45, 23.72, 24.7, 22.7, 24.15, 24.46, 22.39, 22.72, 25.29, 23.6
        ]

        n_total_ci_width_maxtaxa_geiger = [
            2.63710970372773, 1.7170071890568, 1.8976534660395, 2.62720287743294, 2.12445238576767, 1.87063274726802, 2.03662883424665, 2.14943056062518, 2.0842897252997, 2.13665427084761, 2.28388109121028, 3.02635434912448, 2.151639134251, 2.15022926510762, 2.18671554823413, 1.86834641997958, 2.69262337447473, 2.22166870908542, 2.13629101751497, 2.00683191708582, 2.06964222116683, 2.3876133890045, 1.93179881253748, 2.2177108430619, 1.97879223443764, 2.08213364692, 1.63821889929893, 2.05609851550322, 2.27221713431057, 2.41240342506234, 2.57181079656598, 1.81978310129465, 2.58513163418434, 2.01324989752862, 2.04128389624043, 2.02570786370861, 2.44141957678294, 2.35223426442765, 1.70024101868416, 2.22062049168263, 2.03480559158745, 2.06611155633773, 2.33895987796009, 2.49972213372978, 2.31461130921096, 2.30087705719433, 2.05498095053718, 2.24349550496721, 2.09821093353747, 2.06141091076726, 1.80818027661162, 2.06058531215936, 2.33911167425723, 2.08766790617297, 1.97641898637333, 2.4477348559124, 1.97314533256786, 2.61438068406499, 2.40820152316238, 2.1479758853489, 2.04172012051237, 2.20201266099719, 2.029960774264, 3.00213662661253, 2.36721093152753, 2.42452615385091, 1.74826633155838, 2.12237276174521, 2.30010451549972, 2.21942491996317, 2.5355455680537, 2.07795017060447, 2.20652281871159, 2.19821976640511, 1.91052504328263, 2.17145095316938, 2.40224836233547, 1.95084283096079, 2.07023929559124, 1.89423140245869, 2.27856710990557, 1.79951639564129, 2.07241523821747, 2.14911280130366, 1.59125841110396, 2.23065608720707, 1.96069676504156, 1.62116092135522, 3.69298127870268, 2.08221937384323, 2.06928220610067, 2.1347619519551, 2.27714380731475, 2.23726019311952, 2.25158852117497, 2.34694939205959, 2.01079574338523, 1.76148606499504, 2.16745151621703, 2.13092674864267
        ]

        root_ages_mean_maxtaxa_geiger = [
            5.06659734743728, 4.23803293194664, 4.98279507195696, 4.3155903334607, 4.09430571259895, 4.60255582668252, 4.38215291960212, 4.34193016088942, 4.29046521064659, 4.37358409460006, 4.8496227481968, 4.92755985432126, 4.80410465532037, 4.36470691722471, 5.1024650945641, 4.56608895668551, 4.60126185356087, 4.69078233506166, 4.68146837210135, 4.58463378859495, 4.28327681215642, 4.61366929046375, 4.50586550974126, 4.64262101324975, 4.43507928338528, 5.11208300208071, 4.22942116869671, 4.22841445002169, 4.7596647970385, 5.18005203373315, 4.88519115561629, 4.63317656881595, 4.77547403365298, 4.12214707587524, 4.4642409620764, 4.70383201067341, 5.02039329287341, 4.96452047729738, 4.18188468312562, 4.19669714421154, 4.68312489919795, 4.99188418110751, 4.45279408918267, 5.04400772907278, 4.15571708973781, 4.66574290297132, 4.44352325560619, 4.55470212395163, 4.81210010116322, 4.52121024270688, 4.17523488823645, 4.78262229842709, 4.67687750561767, 4.28088858409093, 4.48603265546268, 4.60337922469557, 4.57962101809517, 4.86502882097972, 4.50685077739471, 4.25392899210395, 4.06595950314072, 4.2162352729703, 4.56076145309679, 4.86110226917551, 4.87321560375039, 4.42503618639634, 4.62032411301255, 4.31221733692068, 4.54743762333076, 4.57887648534795, 5.029688314397, 4.61104938167209, 4.90061936330257, 4.87571755247112, 4.38844111299254, 4.2965191578961, 4.63190153444414, 4.35331294706652, 5.0438855935475, 4.31053219610282, 4.82956399153212, 4.33136232493745, 4.67791916630722, 4.66232082345522, 4.22739962148843, 4.56228247455518, 4.44151585935606, 4.14159283267207, 4.55442468236585, 4.44182533888816, 4.16430081481615, 4.70477149383738, 4.86070433211574, 4.62275589866935, 4.49604585852457, 4.83996701890875, 4.35056924907805, 4.23273410429556, 4.93928962704673, 4.42127167578148
        ]

        root_ages_ci_width_maxtaxa_geiger = [
            0.820777733516174, 0.519051816334858, 0.650115479894051, 0.729174571091599, 0.63003935616219, 0.560508790800951, 0.633323455420104, 0.593571433952757, 0.558863434822771, 0.541641997336418, 0.644917025637329, 0.75826651788224, 0.645093205274147, 0.587498602023317, 0.612000655554466, 0.524589984692465, 0.695883176078397, 0.643926201041047, 0.56693867039975, 0.561522377123289, 0.566111924337786, 0.59973043787338, 0.548443390787214, 0.621840073870249, 0.523157853990445, 0.576999298365843, 0.516486965679786, 0.53646158051435, 0.625945307981113, 0.715335242667575, 0.688103781984694, 0.553142496458673, 0.663554653984096, 0.521324441464281, 0.699772807627341, 0.628288816192392, 0.608721537968437, 0.671767028550187, 0.518968863112461, 0.566616541416065, 0.546744744655922, 0.586124184835439, 0.67187197135834, 0.618481258119774, 0.703861756813563, 0.703714276243353, 0.562325051228904, 0.567748946528891, 0.606336778085757, 0.619280367027287, 0.534769319762792, 0.538334079371305, 0.653524344612808, 0.530269129029457, 0.587050330734391, 0.590163979560353, 0.65153138343531, 0.71131501970154, 0.566213258161102, 0.512791101283866, 0.501252999912154, 0.518786773140467, 0.506937212622307, 0.788959603132124, 0.637917810664974, 0.604268265862527, 0.529522506384781, 0.604298264368049, 0.602324917232518, 0.584741249565512, 0.71331567463174, 0.563257059756396, 0.655174840083523, 0.640855275420817, 0.550517443971482, 0.566488138883138, 0.59051170537572, 0.549347331263769, 0.676277031118791, 0.489354324998641, 0.649039739012612, 0.54730016374655, 0.590174084180873, 0.619878923828303, 0.485324606262793, 0.569736052092472, 0.592837049661994, 0.510201955596929, 1.01261229895688, 0.572320184424498, 0.552352472738077, 0.635452389868142, 0.633159479334567, 0.593738131736159, 0.527228454869711, 0.64796121462667, 0.575003163541575, 0.511069188135332, 0.619037931594416, 0.606935961514465
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

        print("\nPJ global mean total taxon count = " + str(global_mean_total / 100.0))
        print("geiger global mean total taxon count = " + str(statistics.mean(n_total_mean_maxtaxa_geiger)))
        print("PJ global mean root age = " + str(global_mean_root_age / 100.0))
        print("geiger global mean total root age = " + str(statistics.mean(root_ages_mean_maxtaxa_geiger)))

        print("\n95% CIs of simulations here and from geiger overlapped " + str(n_total_ci_overlap_count) + " times for the total taxon count.")
        print("\n95% CIs of simulations here and from geiger overlapped " + str(root_age_ci_overlap_count) + " times for root age.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(n_total_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_geiger) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        self.assertAlmostEqual(root_age_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_geiger) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)


    def test_tree_size_total_count_max_t_bd(self):
        """
        Test if birth-death trees simulated here have similar total number of taxa and
        number of extant taxa as trees simulated with phytools, starting from the root
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
            sse_sim = distsse.DnSSE(self.event_handler, stop_condition_value, n=100, stop=stop_condition, origin=start_at_origin,
                                    start_states_list=start_states_list, epsilon=1e-12, runtime_limit=3600,
                                    condition_on_speciation=True, condition_on_survival=True,
                                    debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)


        # phytools from root
        n_total_mean_maxt_phytools = [
            15.56, 15.12, 15.54, 15.09, 14.12, 15.41, 16.01, 14.33, 15.28, 14.87, 16.61, 15.49, 15.02, 15.18, 15.12, 14.51, 14.87, 14.51, 15.05, 16.22, 17.31, 15.43, 13.03, 14.64, 15.57, 13.88, 13.67, 15.28, 13.53, 15.15, 13.54, 14.49, 14.23, 13.65, 13.74, 14.56, 15.88, 16.32, 14.09, 16.82, 15.94, 15.28, 15.78, 15.41, 15.31, 18.09, 16.55, 14.44, 13.84, 15.74, 13.78, 13.42, 13.71, 14.32, 15, 15.88, 15.71, 14.16, 14.29, 15.94, 13.91, 16.13, 15.82, 16.38, 15.08, 15.42, 14.44, 15.78, 14, 15.43, 16, 15.48, 14.08, 14.66, 14.13, 15.29, 15.95, 15.94, 17.01, 14.77, 15.45, 14.43, 16.11, 15.62, 14.77, 14.56, 14.66, 14.28, 14.33, 14.94, 15.2, 12.27, 13.76, 14.63, 15.01, 14.07, 14.86, 17.61, 17.02, 14.42
        ]

        # tree sim with origin
        # n_total_mean_maxt_treesim = [
        #     12.86, 11.36, 11.82, 12.37, 10.83, 11.74, 12.82, 11.53, 11.21, 10.57, 12.51, 10.12, 9.89, 12.38, 10.78, 12.45, 10.51, 11.57, 10.66, 12.68, 12.31, 11.16, 11.01, 12.15, 11.18, 10.53, 12.43, 12.59, 11.27, 12.36, 10.71, 12, 12.07, 12.22, 11.3, 9.85, 11.14, 11.18, 12.53, 11.29, 11.5, 11.02, 11.6, 10.7, 11.61, 10.94, 10.93, 12.57, 11.72, 11.35, 12.33, 12.31, 11.45, 13.19, 12.18, 12.7, 10.52, 13.22, 11.35, 12.27, 10.53, 11.73, 10.48, 11.39, 11.83, 11.68, 11.78, 12.35, 11.3, 10.47, 11.54, 10.5, 10.9, 11.3, 13.14, 11.88, 12.34, 12, 11.34, 12.46, 10.53, 13.04, 10.54, 11.33, 10.7, 13.26, 12.52, 9.68, 13.07, 10.53, 10.63, 9.87, 11.21, 11.07, 11.19, 10.61, 10.45, 12.46, 11.6, 12.7
        # ]

        n_total_ci_width_maxt_phytools = [
            1.77764740956335, 1.91520396951565, 1.63925723312836, 2.56344435082917, 2.13268419328386, 1.7878341350767, 1.76066969581968, 1.87775505722911, 2.23040209411087, 2.04925522573495, 2.22882878007915, 2.04848345146088, 1.74239680690917, 1.73958511528851, 2.1495388769302, 1.45505708784532, 2.16465683528486, 2.31529856477311, 1.91612258775714, 1.61233805463502, 2.10608472015844, 2.23595551029988, 1.69291857292235, 1.84877067493101, 2.27308365784358, 1.98699572584982, 1.79185363038594, 2.13155580909174, 1.86980281097211, 2.22976785492766, 1.81442836989917, 2.02048404801494, 1.96787466587841, 2.00567142338566, 2.01973009858028, 2.06034519575157, 2.27540414312391, 1.93410242215709, 1.94507389343158, 2.14868393227298, 1.62806848440076, 1.78642278570694, 1.8984967765327, 1.86202500780655, 2.12198511976047, 1.833359206429, 1.82254452429297, 1.99422393198685, 1.68772464171814, 2.50582854331354, 1.83349042809036, 1.71331040003517, 2.2409264069157, 2.04128674768211, 2.07920563075168, 1.56139832486425, 2.02600379909303, 1.77724134613852, 1.54245546025442, 1.6419299587665, 1.58923310442553, 2.12667868867125, 1.77674674036395, 1.85033579960858, 2.25462398038549, 1.97722281511368, 1.63175387874767, 1.9580776633621, 1.48408804770247, 2.36773542713904, 1.9372377456658, 2.20016861943338, 2.10952365191287, 1.79083551811139, 2.07872687183208, 2.03227245290298, 1.79553134377115, 1.66941572601393, 2.00076672616406, 1.95382418943781, 2.00181375979108, 1.79984305443726, 2.17959605709893, 2.13278699216293, 1.87401611368297, 1.79796314655324, 1.80339366548939, 2.19818534379077, 1.85513516031321, 1.81956132397816, 2.28072294191508, 1.81299491672779, 1.9281626567038, 1.8903156728521, 1.73913335017674, 2.05131912227575, 2.47118254569374, 2.01542671983635, 1.79049637990828, 1.47685917563133
        ]

        # n_total_ci_width_maxt_treesim = [
        #     2.00598965887543, 1.53334677411764, 1.84672520462957, 1.81520559806109, 1.38622464458418, 1.74856040019535, 1.8082929395293, 1.58998374478699, 1.40740440743256, 1.46520900070863, 1.60626173647404, 1.56311217974419, 1.4053626464656, 1.72940587964933, 1.6328771241665, 1.79428933626164, 1.43692696734905, 1.55613946091544, 1.56483407659494, 1.68602813449332, 1.7642628535923, 1.35498704395759, 1.59875522235585, 1.82325125293064, 1.52849538534806, 1.45276179531931, 1.82641540010409, 1.99630876298525, 1.49135343371128, 1.8181671330316, 1.51136850636111, 1.91918720169543, 1.65451954074501, 1.75254154999088, 1.52165603567596, 1.26421573085507, 1.48644461121907, 1.71376331168713, 1.84253373681817, 1.26756310161098, 1.78586227935081, 1.6572671656018, 1.34038318055214, 1.51501121774691, 1.41279805701275, 1.34264515886598, 1.38488034767403, 1.81148214376287, 1.78468420815232, 1.55940764459464, 2.02658883861727, 1.73565808623796, 1.55666801635137, 1.693105371914, 1.57451484963374, 1.86119227563105, 1.57108544913536, 2.25110677730495, 1.34049173850871, 1.66191427704392, 1.30293585568914, 1.3732311647786, 1.32743208983786, 1.65724726317865, 1.62202238697201, 1.45794108373293, 1.65201399754385, 1.80011899718891, 1.61493639365628, 1.40662116313016, 1.24505394568941, 1.17451421292534, 1.59608435006589, 1.53359476032259, 2.04032178250432, 1.92469504067423, 2.0369336606227, 1.66567970777364, 1.66129891690602, 2.0020008109061, 1.55718642468143, 1.71998777726453, 1.43393305790244, 1.42378571251137, 1.77299614019075, 1.92677857449981, 1.72754253952213, 1.42618886136543, 1.8242895618998, 1.47581701583984, 1.44262988302403, 1.08777547519902, 1.79556916341478, 1.47681582173036, 1.32544718674434, 1.56276831761636, 1.26175780864505, 1.63734818076379, 1.96810043268409, 1.74919941283394
        # ]

        n_extant_mean_maxt_phytools = [
            6.66, 6.55, 6.63, 6.03, 5.65, 5.91, 6.69, 5.54, 6.23, 6.01, 6.98, 6.16, 5.98, 6.38, 6.25, 5.87, 6.21, 5.7, 6.25, 6.92, 7.05, 6.15, 5.62, 6.13, 6.59, 5.69, 5.35, 6.17, 5.11, 6.61, 5.43, 6.09, 5.6, 5.44, 5.29, 5.99, 6.56, 6.84, 5.47, 7.11, 6.66, 6.61, 6.5, 6.76, 6.04, 7.59, 7.1, 5.79, 5.46, 6.65, 5.74, 5.41, 5.65, 5.71, 6.28, 6.73, 6.5, 5.63, 5.48, 6.68, 5.61, 6.34, 5.94, 6.9, 6.03, 6.23, 6.09, 6.44, 5.73, 6.55, 6.27, 6.54, 5.32, 5.56, 5.75, 6.14, 6.79, 6.73, 7.07, 6.17, 6.52, 5.8, 6.62, 6.36, 5.96, 6.46, 6.22, 5.98, 5.66, 6.1, 6.46, 4.65, 5.74, 5.74, 5.85, 6.15, 6.61, 7.19, 6.68, 5.67
        ]

        n_extant_ci_width_maxt_phytools = [
            1.29419979972806, 0.985872305128916, 1.15826420530591, 0.987083852015162, 0.991758747953056, 1.01179262377987, 1.39449697038181, 0.955102561735829, 1.19947542136746, 0.966040005027292, 1.19562758703063, 1.15156623708054, 0.933969649756707, 1.02045563679499, 1.00651700673584, 1.03775354190548, 1.06965978467431, 0.894721942021852, 0.879356307279089, 1.11804196520506, 1.39607662375456, 1.17744271772998, 0.998934022553112, 1.05591684015993, 1.22293967342976, 0.901718349720256, 1.03315155781431, 1.09962565604589, 0.877819317365136, 1.04677820598504, 1.00685621414361, 1.13540475157496, 0.984937058954028, 1.24813562070665, 0.941548313506358, 1.0324452089554, 1.27123899731043, 1.21677642941993, 0.96923211592758, 1.45421412231327, 1.19602996077469, 1.22509542254611, 1.09341322694052, 1.25210876605038, 0.951235058099509, 1.28453163668964, 1.40470119760234, 1.06347481456409, 1.08411821190022, 1.09052594113246, 0.998871867901408, 0.707563317936203, 0.980346403424356, 1.07183419417276, 1.14289373367109, 1.11183820565564, 1.22686948257454, 0.908560527383395, 0.89644640394812, 0.979358335386445, 0.830574710606451, 1.24996856138242, 0.942999953083324, 1.09021454658051, 1.03954682511956, 0.997019048681791, 0.965879319061082, 1.07834011163623, 1.10588915439938, 1.27430454790371, 1.10131819712249, 1.03091249880968, 0.745372879094007, 1.03874956995943, 1.06277401457723, 0.983367789021406, 1.096883857281, 1.33250154115172, 1.26632572827887, 1.31497175819425, 1.07731765365127, 1.01481980288287, 1.19296331219465, 1.0298203463996, 1.08787000411749, 1.0022846499852, 0.985165537851227, 0.931889962724281, 1.08440451954988, 1.14676454832408, 1.03803955394559, 0.844235282442763, 0.879464413691364, 1.087406198625, 0.939931375100858, 0.924110339294999, 1.17629528279392, 1.28923558214262, 1.0254249681242, 0.807914998029127
        ]

        # root_ages_mean_maxt_treesim = [
        #     2.41134926469317, 2.38920913840397, 2.43975399368484, 2.36167155882785, 2.3587007364241, 2.33778283109888, 2.4577397815274, 2.42626338160149, 2.39713681308128, 2.4084376193809, 2.4974528389355, 2.40511012325899, 2.30647267157316, 2.42565744285912, 2.26765674607618, 2.4834094327036, 2.24709301128232, 2.45818899585454, 2.40231025551501, 2.42716561100311, 2.51587695830744, 2.43829360327129, 2.36654427609215, 2.29366342006442, 2.44129237874595, 2.44364909306366, 2.45982058094011, 2.3500470152416, 2.418030790723, 2.36690045733178, 2.33344159000591, 2.3488046990303, 2.43983650212711, 2.45543345103501, 2.42322554561376, 2.37134138005089, 2.39025548763228, 2.38353641595131, 2.41053305434273, 2.45320698395875, 2.56030978698553, 2.38549730708735, 2.40636439184482, 2.33156546281527, 2.45372703782439, 2.43273393188503, 2.43442793638736, 2.45165518183406, 2.35961928551049, 2.40281359050419, 2.3125562566843, 2.39925423527952, 2.39144320514972, 2.4277004420588, 2.44469274408377, 2.41753460371012, 2.31053533778282, 2.34655967858245, 2.36056707080135, 2.4257683669853, 2.36392051223435, 2.4433124600778, 2.3708687917585, 2.37442558977751, 2.40860605171841, 2.4345556612665, 2.35992545428336, 2.44112488756847, 2.38450054327198, 2.44956570399421, 2.46092076374808, 2.43811520578555, 2.29448856565299, 2.37649818075683, 2.38469577201613, 2.4643347062554, 2.44166991306082, 2.45571811529939, 2.48262532651773, 2.4040188589081, 2.49845301372491, 2.52427137935501, 2.30222139001482, 2.31649414607209, 2.40718261057246, 2.47965263780969, 2.52176339091588, 2.34595828935434, 2.41761244731798, 2.32563882473638, 2.35136820012941, 2.43721756305005, 2.37140250618802, 2.34120221673433, 2.42689883102106, 2.30395083398858, 2.32984965307698, 2.48754373987955, 2.35823136681131, 2.43976798550804
        # ]

        # root_ages_ci_width_maxt_treesim = [
        #     0.118614537742372, 0.0973025403085302, 0.0924813356752944, 0.110194862975739, 0.114440741277904, 0.110178798503782, 0.0957698644204067, 0.116719073615701, 0.110946947854556, 0.106754991316681, 0.0871260862329165, 0.117324236823729, 0.12459877289986, 0.119319780757935, 0.135223244008495, 0.0845561395512595, 0.125747814556681, 0.109656611487338, 0.122301468085595, 0.0910437169307571, 0.110916473580753, 0.117261807568291, 0.105129241091335, 0.130308337844548, 0.117198030247258, 0.111154148813296, 0.0899255607893612, 0.12119151921046, 0.116914638975389, 0.121049494696862, 0.120881228390186, 0.127553739150464, 0.112069672234155, 0.0981225418422778, 0.10700017134518, 0.124603154132274, 0.13080762953373, 0.111977296852409, 0.105382701132416, 0.101228706099052, 0.0880208207778674, 0.122068533115286, 0.0958017554961159, 0.122574672405576, 0.0934483208090163, 0.104521531144681, 0.105018649124541, 0.110618124060909, 0.129357799956822, 0.0991629983683035, 0.113002372302917, 0.120914396167643, 0.0998846417472681, 0.112128717510314, 0.0962641661378887, 0.107585124592922, 0.1271274597584, 0.120671574088347, 0.117222788974941, 0.113849859692386, 0.0988910493234009, 0.109599399979178, 0.109263277034012, 0.12104537332385, 0.118114588466149, 0.0992280983891192, 0.111306108533451, 0.101338565822798, 0.0992505816563806, 0.107205384355046, 0.0846739513590895, 0.101923577275351, 0.120585421036524, 0.13533080458462, 0.119862388629098, 0.107469642266994, 0.10187092799781, 0.102821551426194, 0.0848564790648721, 0.109156482672756, 0.110050492222876, 0.0881731057795677, 0.113992856348529, 0.121231280141388, 0.104053085645244, 0.10507398152035, 0.0917024098969028, 0.119101188258717, 0.104633851513134, 0.123410247752883, 0.123102428889098, 0.091508461525494, 0.115638713124, 0.117966708499831, 0.103408689657657, 0.111697405092111, 0.118632053538949, 0.0972386360773996, 0.113507182558513, 0.0976551271056787
        # ]

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
            # print("PJ's vs. TreeSim mean total taxon count: " + str(mean_total) + " <-> " + str(n_total_mean_maxt_treesim[i]))
            # print("PJ's vs. TreeSim mean root age: " + str(mean_root_ages) + " <-> " + str(root_ages_mean_maxt_treesim[i]))
            global_mean_total += mean_total
            global_mean_extant += mean_extant

            stdevs_n_total = statistics.stdev(n_total)
            stdevs_n_extant = statistics.stdev(n_extant)
            # stdevs_root_ages = statistics.stdev(root_ages)

            sterr_n_total = stdevs_n_total / math.sqrt(n_sim)
            sterr_n_extant = stdevs_n_extant / math.sqrt(n_sim)
            # sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)

            n_total_ci_width_maxt = 1.96 * sterr_n_total
            n_extant_ci_width_maxt = 1.96 * sterr_n_extant
            # root_ages_ci_width_maxt = 1.96 * sterr_root_ages

            if abs(mean_total - n_total_mean_maxt_phytools[i]) <= (n_total_ci_width_maxt + n_total_ci_width_maxt_phytools[i]):
                n_total_ci_overlap_count += 1

            if abs(mean_extant - n_extant_mean_maxt_phytools[i]) <= (n_extant_ci_width_maxt + n_extant_ci_width_maxt_phytools[i]):
                n_extant_ci_overlap_count += 1

            # if abs(mean_root_ages - root_ages_mean_maxt_treesim[i]) <= (root_ages_ci_width_maxt + root_ages_ci_width_maxt_treesim[i]):
            #     root_age_ci_overlap_count += 1

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time

        print("\n\nPJ global mean total taxon count = " + str(global_mean_total / 100.0))
        print("phytools global mean total taxon count = " + str(statistics.mean(n_total_mean_maxt_phytools)))
        print("PJ global mean extant taxon count = " + str(global_mean_extant / 100.0))
        print("phytools global mean extant taxon count = " + str(statistics.mean(n_extant_mean_maxt_phytools)))

        print("\n95% CIs of simulations here and from phytools overlapped " + str(n_total_ci_overlap_count) + " times for the sampled ancestor count.")
        print("\n95% CIs of simulations here and from phytools overlapped " + str(n_extant_ci_overlap_count) + " times for the sampled ancestor count.")
        # print("\n95% CIs of simulations here and from TreeSim overlapped " + str(root_age_ci_overlap_count) + " times for root age.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(n_total_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_treesim) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        # self.assertAlmostEqual(root_age_ci_overlap_count, exp_count,
        #                         msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_treesim) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)


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
            sse_sim = distsse.DnSSE(self.event_handler, stop_condition_value, n=n_sim, stop=stop_condition, origin=start_at_origin,
                    start_states_list=start_states_list, epsilon=1e-12, runtime_limit=3000,
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
    # rates_t0_s0 = [ sseobj.MacroevolStateDependentRateParameter(name="lambda", val=1.0, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
    #                 sseobj.MacroevolStateDependentRateParameter(name="mu", val=0.8, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]) ]

    # matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

    # fig_rates_manager = sseobj.FIGRatesManager(matrix_atomic_rate_params, total_n_states)

    # event_handler = sseobj.MacroevolEventHandler(fig_rates_manager)

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