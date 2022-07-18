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
            3.82808365734507, 4.00096314682251, 3.86032478422541, 4.05144732474238, 3.94573176174099, 4.22591940164165, 4.1858163386353, 3.85794143977763, 4.03797066101873, 4.05021841612791, 4.07014874709922, 3.96857489451474, 3.90608941309632, 4.32673571009221, 3.92753144479358, 4.03920522876101, 4.05028270249215, 3.82048173999274, 4.54022036490761, 3.57027605298589, 4.04380390175176, 3.82582120164118, 4.01474796262255, 4.03757673378799, 4.18353722826596, 4.27080992265196, 4.01707705621916, 4.19615567623368, 4.1045525724648, 4.0682372819021, 4.01928660677234, 4.37300674696915, 4.13921068303988, 3.795556414833, 4.1871457201025, 4.29653373720795, 4.07751116627911, 4.18038790733958, 4.17124014111614, 4.09148604256178, 4.0734121335006, 3.96331555575914, 4.47655706470497, 4.00022462160528, 4.20645520041668, 4.04880727047, 4.2459505989754, 3.88101504000803, 4.33608656882123, 4.19712909175085, 4.35953417991503, 4.02753978924281, 4.19721834119377, 3.93571430602648, 4.071225105414, 4.00450032805849, 3.80178147110162, 3.99193490912177, 4.29510191880387, 4.44211957297744, 3.8321177047849, 3.99575432256419, 3.85983102987576, 4.09734719995337, 3.98331501930357, 4.12040067352513, 4.52864800502805, 3.93538514352664, 3.94067286481623, 3.79763464790409, 3.98295470261725, 3.91509309156491, 4.02834885910831, 3.88048575370872, 4.01636606103248, 4.33101812478601, 3.32559072708332, 4.27389768864807, 3.81432643351652, 4.06658090480438, 4.22527775492514, 3.84336081551627, 3.8231414612518, 4.06636936266563, 3.94713932590288, 4.02397261320875, 4.07373563886444, 3.89002554079445, 3.938875203392, 4.57718489604075, 4.31790350658806, 3.9697582858912, 4.28141164653263, 4.2270783697956, 3.90703635106247, 4.28625954978261, 3.81801322965817, 4.33897038304726, 3.77574014321606, 4.17931249155088
        ]

        root_ages_ci_width_maxtaxa_fossilsim = [
            0.411272864550513, 0.495456124442128, 0.426652955380913, 0.42874565535312, 0.298752275494362, 0.465915355508509, 0.428956017225259, 0.379511818408814, 0.445776150381127, 0.47388810594271, 0.425032402303995, 0.413243614909602, 0.491651502882015, 0.486979141180629, 0.431630962956839, 0.431060583780089, 0.464831791268453, 0.443672234098104, 0.428301688593243, 0.361360644967735, 0.415602165842102, 0.467680116859945, 0.488803706003131, 0.38624279965205, 0.441932032445108, 0.438184743975959, 0.413224950150033, 0.443198015604551, 0.502637939733523, 0.477170241144126, 0.523782582281256, 0.414164133138349, 0.412169858700465, 0.468444792904064, 0.505730499572329, 0.526299584008709, 0.451390069237767, 0.522428048515174, 0.527481494295083, 0.442930644752631, 0.442648337409031, 0.465746978666139, 0.540527499552143, 0.410227907374619, 0.45914413549892, 0.480053697090059, 0.472367492691222, 0.514638786081069, 0.453481051487473, 0.523854315054489, 0.479785432214185, 0.420785467307656, 0.520227427361694, 0.374921518328542, 0.421453268903336, 0.490437349568258, 0.37325955245563, 0.51389351273983, 0.495377530719814, 0.583763514261502, 0.404308960676997, 0.428899562181744, 0.386146966569531, 0.400889359631582, 0.414430344755516, 0.451526883663354, 0.612853839616445, 0.392441426270844, 0.388644140559108, 0.429859643127227, 0.504598752159916, 0.398713601489738, 0.416956212820413, 0.458306304084857, 0.41956156849734, 0.503139655672981, 0.31439664095705, 0.490814487298686, 0.48903242883895, 0.433026024096485, 0.411947796313868, 0.526749395810516, 0.405621349300893, 0.36977947211838, 0.383589654305315, 0.4104445712992, 0.453663876413896, 0.364774054451824, 0.435626104102587, 0.551865932732791, 0.464958345594286, 0.493926295619062, 0.464239178321293, 0.39593158047593, 0.464686660093282, 0.441881834071931, 0.460238995172455, 0.396052689126035, 0.391689589300043, 0.466519347970454
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

    def test_tree_size_sa_count_max_t_fbd(self):
        """
        Test if FBD trees simulated here have similar root ages and number of sampled ancestors
        as trees simulated with FossilSim
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
            9.56, 8.73,  8.35, 9.66,  9.86,  9.6, 9.49, 8.24,  8.15,  8.41,
            9.31, 9.07,  8.43,  8.8,   6.8, 9.44, 8.42, 9.07,  9.23, 10.48,
            9.77, 7.74, 10.14, 9.47,  9.18, 9.19, 8.49, 9.78,  10.0 , 8.84,
            7.85,  9.2,   9.0, 8.37,  8.69, 9.59, 8.72, 9.04,  9.41,   8.7,
            7.57, 7.82,  9.13, 9.79, 10.82, 9.63, 8.17, 8.86, 10.01,  8.02,
            7.87, 9.29,  9.34, 9.22,  8.65,  9.6, 9.82, 8.42,  8.25, 10.01,
            9.47,  9.4,  9.95, 7.54,  8.73, 9.75, 9.93, 9.48,   8.5,  10.7,
            9.47, 7.48,   8.4, 9.49,  8.81, 8.12, 7.91, 7.54,  8.96,  8.33,
            8.65, 9.75,  9.68, 8.55,  9.57, 8.53, 8.82, 10.95, 8.79,  9.13,
             9.8, 8.11,  8.95, 9.14,  9.83, 8.81,  8.2,  9.5,  7.66,  9.24
        ]

        n_sa_ci_width_maxt_fossilsim = [
            1.57622429004417, 1.22349799634928, 1.28763334702222, 1.67618115459824,  1.52891167683276, 1.51385819225667, 1.59681232915626, 1.12131353277657, 1.21379257043385, 1.24245662095282,
            1.32456860993466, 1.64605467759007, 1.29899868259122, 1.30339590642632, 0.848647741939487, 1.24595744791714, 1.19990399346587, 1.33610768996811, 1.24331831723811, 1.33733417686308,
            1.50146668819444, 1.21099042794138,  1.3965337777237, 1.44499497940479,  1.50186075885546, 1.15492259673298, 1.25377342305683,  1.4591914827857, 1.67983029445888, 1.29494915870715,
            1.44201647533983, 1.24429948650685, 1.38312660010542, 1.41025790663705,  1.54699596233305, 1.52027707150421, 1.25642765414528, 1.50921671248781, 1.33140613388002, 1.33006623971147,
            1.16780643961405, 1.14469179902167, 1.27712728837304, 1.47918921771969,  1.57155465766234, 1.50296489794155, 1.17568814046769, 1.27023128666246, 1.34226216349697, 1.06736643410146,
            1.10230431259964, 1.64818438656639, 1.29119801484333, 1.80403906696931,  1.29514541908067,  1.3808803516071, 1.40770352480108, 1.14156882022642, 1.33759384453177, 1.75230129737088,
            1.44392041583763, 1.78966072062542, 1.46603505683827, 1.00344544546288,  1.31402711778054, 1.39746568595953, 1.62994555328289, 1.60480797605196, 1.25840887190487, 1.52751004890146,
            1.56836003442314, 1.27585813213434, 1.18028176742041, 1.45277247949282,  1.23701036733282, 1.29301792842649, 1.13813557907819, 1.29275981285673,  1.6672910268099, 1.32438110478331,
            1.46471102400384,  1.5246812947318, 1.64643771000585, 1.22080554627598,  1.46732616081845, 1.28191894731092, 1.33846532108889, 1.66832754819934, 1.51239514851527, 1.42067534431319,
            1.68029223047599, 1.24239415596979, 1.56784532065346, 1.43922729212027,  1.48148547400809,  1.3460727997483, 1.20146173598259, 1.44420122272389, 1.06890677202117, 1.49142238338325
        ]

        root_ages_mean_maxt_fossilsim = [
            2.56346380658813, 2.54575290712475, 2.56346944822295, 2.50736951884695, 2.53371320448893, 2.50663132679238, 2.52451105622246, 2.56238717276736, 2.50117846574106, 2.54348146246823,
            2.56439836873595,  2.5335756333469, 2.58650508477994, 2.61987224021307, 2.48662656353868, 2.58589449023329, 2.51526970489459, 2.63211196305548, 2.54148770760231,  2.6286950967592,
            2.50748682094593, 2.51810858763158, 2.58379144765191,  2.5995616947533, 2.59910254171525, 2.58125005952959,  2.5506504799352, 2.56562681186241, 2.55020412697261, 2.58662683848908,
            2.54740911663707, 2.56432699538845, 2.48823635460408, 2.49308503880128, 2.57036141890845, 2.53510169586257, 2.59490821977655, 2.55918390556403,  2.5679582206598, 2.50778787807515,
            2.53760142981706, 2.50940849634466, 2.62283536731434, 2.59968961205811, 2.60610912066357, 2.53458771204322, 2.59109859643892, 2.58537956201953, 2.61696768161528, 2.49940230417064,
            2.59052545588793, 2.52230018284912, 2.56792515254857, 2.51123422846728,  2.4513116655411, 2.55590202988248, 2.51462513061637, 2.57543484450858, 2.49834991805379, 2.56664348846925,
            2.58311131867186, 2.57243176738569,  2.5434887484404, 2.55814945720304, 2.56139385603142, 2.60596851717318, 2.57852725818135, 2.58377353970592, 2.52673120059733, 2.55728507308721,
            2.51418224339663, 2.47217752164821, 2.48179659797191, 2.51688361190748,  2.5494291388368, 2.54451130113411, 2.43042491700309, 2.48850561898614, 2.45710525478957, 2.51325938882337,
            2.51187911418291, 2.51307205834424, 2.52714508334476, 2.59052723854795, 2.56990177734216, 2.55054072227071, 2.55180110243501, 2.59963870734984, 2.58229472717322, 2.59177586674776,
            2.46618271356585, 2.55270693566421, 2.55106715320238, 2.62201598688621, 2.59891491330703, 2.56418323475208, 2.52452513612158, 2.60573934881315, 2.55427035501125, 2.53869540493527
        ]

        root_ages_ci_width_maxt_fossilsim = [
            0.0853142091235561,  0.098168246061348, 0.0942173793706048, 0.0981965715081872, 0.0766534466891648, 0.0735590159077762, 0.0788767039702253, 0.0746817093171096,  0.110406398815396, 0.0984156669395894,
            0.0833398545621408, 0.0736946034959934, 0.0742981878347941, 0.0764099379670089, 0.0918628520383513, 0.0841222431260386, 0.0910869784653398, 0.0755643289218411, 0.0872595345465329,  0.068076901561684,
             0.095206721694312,  0.100481760932599, 0.0801556985965358, 0.0869192312456692, 0.0695828001832828, 0.0798332979819988, 0.0923403035085774, 0.0751232415663398, 0.0930716445890904, 0.0812184184808355,
            0.0787715984442822, 0.0815017557269185, 0.0956307408700714, 0.0939062397125453, 0.0835571745122473, 0.0845821938271109, 0.0726154026981018, 0.0853040761436524, 0.0869252622656047, 0.0831415776771962,
            0.0806582985563312, 0.0967475544109498,  0.082062070397678, 0.0762447831963488, 0.0858767803245255, 0.0889908784841554, 0.0814915066863414, 0.0753148625160793, 0.0807551019021167, 0.0947672969955067,
            0.0816404799303369, 0.0894298561700682, 0.0947075215427731,  0.090385883414368, 0.0976410577503781, 0.0914555680850767, 0.0810708585890818, 0.0723167481077285,   0.09431152248645, 0.0855132995423713,
            0.0871991804742121, 0.0836512475559614, 0.0972628346708553, 0.0860582883132254, 0.0862338452933443, 0.0714805853827436, 0.0690858801235323, 0.0751519457283291, 0.0856333273280159, 0.0943181253655539,
            0.0872614443802423, 0.0947567084093306,  0.105596528217165, 0.0755661405995961, 0.0849260982349754, 0.0665649274722478,  0.104276828620311, 0.0952400203746332,  0.104052483608126, 0.0935856953174392, 
            0.0882439655638875, 0.0845361178012816, 0.0794157797493597, 0.0724138305248923, 0.0801120897138487, 0.0749223618669095, 0.0738595165003429, 0.0781043851627495, 0.0726309039112226, 0.0811696358387731,
             0.113884342818594, 0.0890966795011292, 0.0985459052178513, 0.0730176708305223,  0.079349324787874, 0.0800475349187237, 0.0897966476713933, 0.0714615879561518, 0.0889547717486681, 0.0837896632562448
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
    #     condition_on_speciation=True,
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

    # unittest.main()