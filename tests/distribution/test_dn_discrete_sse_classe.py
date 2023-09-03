import unittest
import math
import statistics

# pj imports
import phylojunction.utility.helper_functions as pjh
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.distribution.dn_discrete_sse as distsse

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class TestClaSSETrees(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        total_n_states = 3

        # calling state 0 "1" to match R unit test
        rates_t0_s1 = [ sseobj.DiscreteStateDependentRate(name="lambda1", val=0.9, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
                        sseobj.DiscreteStateDependentRate(name="mu1", val=0.6, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]),
                        sseobj.DiscreteStateDependentRate(name="q13", val=0.9, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[0,2]) ]
        
        rates_t0_s2 = [ sseobj.DiscreteStateDependentRate(name="lambda2", val=0.7, event=sseobj.MacroevolEvent.W_SPECIATION, states=[1,1,1]),
                        sseobj.DiscreteStateDependentRate(name="mu2", val=0.4, event=sseobj.MacroevolEvent.EXTINCTION, states=[1]),
                        sseobj.DiscreteStateDependentRate(name="q23", val=0.4, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[1,2]) ]
        
        rates_t0_s3 = [ sseobj.DiscreteStateDependentRate(name="lambda312", val=1.2, event=sseobj.MacroevolEvent.BW_SPECIATION, states=[2,0,1]),
                        sseobj.DiscreteStateDependentRate(name="lambda313", val=0.9, event=sseobj.MacroevolEvent.ASYM_SPECIATION, states=[2,0,2]),
                        sseobj.DiscreteStateDependentRate(name="lambda323", val=0.7, event=sseobj.MacroevolEvent.ASYM_SPECIATION, states=[2,1,2]),
                        sseobj.DiscreteStateDependentRate(name="q31", val=0.4, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[2,0]),
                        sseobj.DiscreteStateDependentRate(name="q32", val=0.6, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[2,1]) ]
        
        rates_t0 = rates_t0_s1 + rates_t0_s2 + rates_t0_s3

        # original implementation
        # matrix_atomic_rate_params = [ [ rates_t0_s1, rates_t0_s2, rates_t0_s3 ] ] # 1D: time slices, 2D: states, 3D: parameters of state, several parameters -> matrix
        matrix_atomic_rate_params = [ rates_t0 ]

        state_dep_par_manager = sseobj.DiscreteStateDependentParameterManager(matrix_atomic_rate_params, total_n_states)
        
        event_handler = sseobj.MacroevolEventHandler(state_dep_par_manager)

        cls.sse_stash = sseobj.SSEStash(event_handler)


    def test_tree_size_state_count_max_taxa_classe(self):
        """
        Test if ClaSSE (GeoSSE) trees simulated here have similar root ages and number of tips for the three states 
        (A, B, AB) as ClaSSE trees simulated with diversitree
        """
        
        stop_condition = "size"
        stop_condition_value = [50] # 50 living taxa

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
                runtime_limit=3600,
                condition_on_speciation=True,
                condition_on_survival=True,
                debug=False)

            trs = sse_sim.generate()

            # n0, n1, n2 = float(), float(), float()
            # for tr in trs:
            #     n0 += float(tr.state_count_dict[0])
            #     n1 += float(tr.state_count_dict[1])
            #     n2 += float(tr.state_count_dict[2])
            
            # print(str(n0/100) + " " + str(n1/100) + " " + str(n2/100))

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i, n_batches)


        # "expectations" from diversitree
        n1_ci_width_maxtaxa_divtree = [
            1.33148919498569, 1.31801361100251, 1.40783445479972, 1.27859702088585, 1.42576435017229, 1.45010567584266, 1.50838729265612,  1.3304731619794, 1.45984020379152, 1.39555536897638,
            1.45940420997091, 1.21229068487985, 1.23097432051923, 1.35532349800423, 1.43056810075896, 1.39452062274917, 1.28250000492252, 1.42137577058712, 1.33575913331043, 1.25284458491043,
              1.149974650577, 1.24997942673979, 1.16104988782576, 1.38224397910708, 1.20310616104626, 1.34280555088023,  1.2758018648614, 1.40721828850934,  1.3316931827299, 1.35186190436925,
            1.24986766426763, 1.37476751180523, 1.58350925772252, 1.39291134685665, 1.29779277672251, 1.23516605124227, 1.47226973593616, 1.36936150335428, 1.41180893147405, 1.37424805724013,
            1.35316874652264, 1.38575428900471, 1.51776979183607, 1.25562903805707, 1.33817392537198,  1.4498848935155, 1.14519170182961, 1.24070641266346, 1.20733197396306, 1.35463618374198,
            1.51585109260425, 1.27327640678715, 1.44198418349667, 1.42743309404672, 1.15158645495348, 1.34184580306081,   1.208809524727, 1.36003367456824, 1.39078553488322, 1.48894860823999,
            1.32952786090751, 1.32807942989914, 1.41570783951657, 1.38393856289164, 1.27415837399446, 1.37408568741678,  1.2044600342649, 1.51346211760128, 1.46169438141031, 1.36402226054883,
            1.32271583180862,  1.3198391796738, 1.27176848855633, 1.32769369413828, 1.09642563437616, 1.36987148022113,  1.3400632453038, 1.37660801929521, 1.21631711441888, 1.45814734038184,
            1.10079660778782, 1.57714966834592, 1.46557972450798, 1.17665974759748,  1.3232863064692, 1.25938136552769, 1.33053149188927, 1.27550528011293, 1.23523987675727, 1.23536552706256,
            1.27882765188759,  1.2947618596218, 1.48697575860088, 1.42722647697051, 1.47477543519087, 1.44228554591648,  1.3054945955422, 1.35463761600919, 1.50933369472102, 1.30562537297092
        ]

        n1_mean_maxtaxa_divtree = [
            28.18, 29.65, 29.77, 29.01, 29.44, 29.36, 29.69, 29.39, 29.14, 28.9,
            28.18, 29.08, 28.3, 29.89, 30, 28.62, 28.85, 29.66, 29.17, 30.51,
            27.8, 29.57, 29.02, 29.27, 29.41, 27.85, 29.29, 28.74, 29.72, 29.94,
            29.11, 28.79, 29.02, 29, 29.34, 28.94, 28.98, 29.42, 28.79, 29.97,
            29.75, 29.55, 29.88, 29.51, 29.75, 31.19, 28.77, 30.01, 29.66, 29.51,
            29.62, 29.6, 29.43, 29.47, 29.62, 29.33, 28.94, 28.95, 29.55, 29.74,
            28.63, 28.69, 28.7, 29.89, 28.89, 29.68, 29.29, 30.47, 30.4, 28.45,
            28.85, 29.78, 29.33, 30.85, 29.2, 28.98, 28.11, 29.06, 28.88, 30.13,
            30.45, 30.28, 29.87, 28.6, 28.56, 29.37, 27.41, 28.56, 28.17, 27.53,
            31.07, 29.59, 29.17, 29.31, 29.99, 29.05, 27.67, 28.7, 29.55, 30.99
        ]

        n2_ci_width_maxtaxa_divtree = [
            1.65879776545143, 1.66918326898353, 1.63472968191359, 1.65339926359258, 1.61573393349226, 1.75063633844809, 1.70857483289318, 1.80834444023959, 1.77321389390433, 1.43824420173988,
            1.60323551443557, 1.41262913045216, 1.89178796983919, 1.75472666886448, 1.87436808864139, 1.51010477929329, 1.64944583580442, 1.55246077925051, 1.77783840194612, 1.72791087686369,
            1.62711482897387, 1.71578636517169, 1.69486462616048, 1.81752247842935,  1.7089370405269, 1.56119949533762, 1.59109257931762, 1.82856423226398, 1.83053248575767,  1.8084259800314,
             1.4991906832176,  1.8735677674944, 1.70242146783381, 1.67881821060007, 1.77682208655586,  1.3373399800433,  1.7739632437858, 1.60526732774805, 1.74917168275661, 1.46710532503548,
            1.55779433764019, 1.63065960552499, 1.64628215066433, 1.67018725134203, 1.86435759304579, 1.67088410638491, 1.74586754747699, 1.41967120734686, 1.55535128072878, 1.93495290950974,
             1.7298770089988, 1.59618037967031,  1.6091303509861, 1.60044968226263, 1.50444098145014, 1.59695812813934, 1.44931873428013, 1.58676875607397, 1.57578109744682, 1.68828555150909,
            1.56187044354653, 1.73393238878824,  1.5461943374922, 1.50859308274065, 1.67453667899301, 1.57516534530522, 1.50408112611634, 1.61885304902682, 1.57901231422837, 1.68211448414719,
              1.419764136906, 1.61440648455516, 1.60449360701422, 1.57996184789636, 1.48005073131725, 1.63716687091827, 1.66626201148231, 1.48859282803484, 1.55524774004957, 1.72240007506556,
            1.72272333662299, 1.54959243178536, 1.53735392506381, 1.91992606381211, 1.72206210613141, 1.63902523384227, 1.58681399672489, 1.56614283957096,  1.6368030056155, 1.58348475244068,
            1.66667066714671, 1.81837520947252, 1.40356122816645, 1.64169360986966, 1.86345822750513, 1.52760022852244, 1.56780695781391, 1.58461650078038,  1.7346349537082, 1.61298291962125
        ]

        n2_mean_maxtaxa_divtree = [
            39.36, 37.67, 37.82, 37.48, 37.94, 39.98,  37.5, 39.74, 38.51, 37.95,
            39.02, 38.12, 39.03, 39.03, 38.96, 38.55, 38.13, 37.64, 38.37, 38.24,
            38.65, 38.56, 38.05,  39.1, 36.41, 39.78,  38.4, 38.55, 38.37,  37.8,
            36.67, 40.17, 38.03, 37.74,  38.6,  37.7, 40.04, 38.18, 41.15, 37.46,
            37.11, 38.43, 38.34, 38.45, 38.19, 36.95, 39.49, 38.02, 39.09, 38.71,
            37.68, 38.61, 39.18, 37.52, 37.95, 36.41, 37.22, 38.71, 38.64, 37.31,
            37.88, 40.02,  38.3, 36.49, 37.76, 39.14, 37.98, 38.06, 35.87, 37.89,
            38.44, 38.79, 38.58, 39.36, 38.78, 38.37,  38.7, 37.43, 37.92, 38.74,
            37.83, 37.67, 38.05, 38.63, 37.76,  37.9, 39.48,  39.3, 37.76, 39.68,
            37.07, 38.01, 37.85, 38.38, 38.35, 38.23, 38.34,  37.7, 39.76, 36.85
        ]

        n3_ci_width_maxtaxa_divtree = [
            0.477928973210815, 0.565086125195876, 0.612551670547807, 0.538100240746005, 0.590055397501974, 0.465424188403992, 0.543096116276992, 0.500813134368055, 0.497856212254912, 0.558555713289628,
            0.504422415573826, 0.466069885079242, 0.564838861798637, 0.474486208461534, 0.478882029714875, 0.505989306390776, 0.527558621760596, 0.462245138515061, 0.453717218923279, 0.460339814756996,
            0.492748172205764, 0.456501060881416, 0.545352809906328, 0.577932329756502, 0.528484589379975, 0.540086893121109, 0.501490644242369, 0.510113727070204,            0.5096, 0.548782322947631,
            0.521518359787994, 0.563986346994286,  0.49072802481969,  0.51555034201474, 0.498385937833626, 0.468167429258147, 0.492531561415001, 0.557860559674552, 0.489588028836284, 0.512076345318101,
            0.567597273646117, 0.503394383614424, 0.534337205226146, 0.494497260872192, 0.575378615444744, 0.532901012658945, 0.586016857757315,  0.44800671712136, 0.535342059629244, 0.522874504131448,
            0.564429952888019, 0.474064848310904, 0.502881508625215, 0.528572692171182, 0.534308156082766, 0.514627492621232, 0.494728698226093, 0.535787654304437, 0.490854527607444, 0.525166370945167,
            0.505977802842869, 0.527102389762112, 0.547508075305974, 0.557342106644011, 0.462392022230296, 0.431087496801972, 0.526583129546734, 0.525077696607605,  0.50059226175127, 0.504826124281991,
            0.482916558449017, 0.499412624687212, 0.486523017404774, 0.529581157161206, 0.481971479939163, 0.571911826219477, 0.499350461331842,  0.53346504282773, 0.558013567772489, 0.453396387771721,
            0.509474343482247,  0.49561030634721,  0.44455917286883, 0.528660780279813,  0.51595662490938,  0.44855638430301,  0.43609916581466, 0.509295324617032,  0.53660180899261, 0.490510521927116,
            0.476302363859536, 0.546106521595643, 0.493377771954297, 0.500115309531995, 0.494979626557569, 0.496767727769765, 0.499913534341872, 0.550945301873576, 0.493192910399863, 0.556366525509102
        ]

        n3_mean_maxtaxa_divtree = [
            8.56, 8.53, 8.48, 8.59, 8.74, 8.24, 8.67, 8.58, 8.35,  8.6,
            8.23, 8.61, 8.41, 8.41, 8.51, 8.61, 8.74, 8.44, 8.93, 8.17,
            8.73, 8.36, 8.66, 8.15, 8.32, 8.73, 8.83, 8.79, 8.74, 8.83,
            9.03, 8.23, 9.21, 8.52, 8.33, 8.46, 8.22,  8.8, 8.27, 8.68,
            8.76, 8.64, 8.61, 8.72, 8.78, 8.96,  8.7, 8.26, 9.12, 8.88,
             8.7, 8.78, 8.73,  9.2, 8.73, 8.57, 8.75, 8.61, 8.53, 8.95,
            8.32,  8.8, 8.57, 9.07, 8.49, 8.53, 8.71, 8.57, 9.11, 9.18,
            8.01, 8.15,  8.2, 8.35, 8.44, 8.97, 8.79, 8.69, 8.66, 8.68,
            8.97,  8.3, 8.37, 8.76, 8.86, 8.93, 8.17, 8.66, 8.86, 8.86,
            8.44, 8.88, 8.87, 9.12, 8.81, 8.98, 9.14, 8.76, 8.54, 8.77
        ]

        root_ages_ci_width_maxtaxa_divtree = [
            0.402583958506349,  0.35349828513856, 0.279136899947445, 0.335964859631695,  0.28499397804662, 0.350514303436294, 0.275736198280265, 0.455849583827377, 0.259865149105153, 0.316546126933438,
            0.305578753591548, 0.303276427111832, 0.356803245109667, 0.349309953728566, 0.297484943682441, 0.265405747167443, 0.363236300308106, 0.295812793909584, 0.285125509359068, 0.310244618451754,
            0.288206675600422, 0.336788789812009, 0.337736473094787, 0.272428124414821, 0.321512199607884, 0.267271492064018, 0.371021397380774, 0.352442900925604, 0.309540218809714, 0.355374796132488,
            0.305623973266641, 0.302257792816758, 0.300967278683057, 0.266956700043462, 0.293744911020102, 0.256051629776374, 0.309583989206934, 0.324916435756151,  0.28976540917237, 0.298235293668219,
            0.319584333981025, 0.271805481279532, 0.314020616539998, 0.313798716586409, 0.283627264545878, 0.296468492871938,  0.45298408408228, 0.306207593301471, 0.280898905689945, 0.341150709685206,
            0.298711003193309, 0.333910311099175, 0.319905289354395,  0.33376613540649, 0.261567434000333,  0.37741849941295, 0.294467075077753, 0.307454704085834, 0.346391306063803, 0.290871500055408,
            0.343205617473415, 0.350740045134609, 0.284444919739236, 0.246979986168107, 0.306765833389505, 0.309553536719499, 0.305167803383524, 0.254944666242712, 0.343337884707819, 0.283345167670977,
            0.297222585093481, 0.333984832861718, 0.286609263763584, 0.308404383388671, 0.349809192412382, 0.414744183506743, 0.296773462063394,  0.27413513569572, 0.351391384158802,  0.26791640108123,
            0.320410244823906, 0.299108070492746,  0.29918894547005, 0.309930051009676, 0.261255787898873,  0.28158359341994, 0.330409421230667, 0.306513711035209, 0.320042174273291, 0.317514715802498,
            0.286335101220367,  0.36542235825688,  0.32304600165077, 0.294266905897255, 0.343132428237696,  0.30837395518183, 0.305561551256897, 0.262667767787128, 0.375723289860407, 0.323511383964665
        ]

        root_ages_mean_maxtaxa_divtree = [
            5.08809704868853, 4.95753861326743, 4.73636284670131, 4.96319404849764, 4.73567730667835, 5.07163432338756, 4.62637335091147, 5.16307332294441, 4.76080447916564, 4.91593770775134,
            4.67908259075116,  4.9157819421004, 4.90971719868318, 4.82118785285073, 4.76865206960806, 4.57380150527143, 4.84420190154827, 4.82351228903679, 4.75515935691569, 4.77691149408655,
            4.77539392318328, 4.90003057814264, 4.80059210833242, 4.79140720380567, 4.62017727971943, 4.79818435964553, 4.89002658190016,  4.8802668875685, 4.82645552580452, 4.98041134655664,
            4.64552965351449, 4.77548921944079, 4.90378277256892, 4.79361323760426, 4.77139635474947, 4.69246861709882, 4.74756168773244, 5.02304518811842, 5.09361908167418, 4.74091267199378,
            4.69201783819571, 4.69175774995683, 4.83038879119169, 4.76759826602195, 4.79497437782467, 4.90607038574299, 5.09813837617416, 5.07873260639862,  4.8204596281595, 4.81986973522004,
            4.85284012770134, 4.80575537576375, 5.00501395688113, 5.04013324930276, 4.75502249254939, 5.02559335276541, 4.73056759503486,  4.8574636391649, 4.95127186276422, 4.81535239685493,
            4.74811431504154, 5.07250958068098, 4.74875653600774, 4.39826003370807, 4.81985408466801, 4.78999193814117, 4.95907610123399, 4.78060163162943, 4.79770442919843, 4.75592090276486,
            4.93917652978045, 4.79632734760973, 4.71875674806222, 5.03755596934894, 5.01843915972029,  5.1153635330767, 4.80821484679762, 4.75670039780254, 4.75572626809253, 4.74002377788157,
             4.8222321094791, 4.81827122346801, 4.92881100627846, 5.01322490498682, 4.76512822262139, 4.68956597661793, 4.99604570258193, 4.87552788348062, 4.82144322511717, 4.92249089073036,
            4.86686357554131, 4.82568772352076,  4.8844156477878,   4.948674394934, 4.86740785611967,  4.8069211295702, 4.77484137489311, 4.77531869830951, 5.14516991704117,  4.8200173812382
        ]

        
        # parsing simulations
        n1_ci_overlap_count = 0
        n2_ci_overlap_count = 0
        n3_ci_overlap_count = 0
        root_age_ci_overlap_count = 0
        for i, batch in enumerate(sim_batches):
            n1s = [ann_tr.state_count_dict[0] for ann_tr in batch]
            n2s = [ann_tr.state_count_dict[1] for ann_tr in batch]
            n3s = [ann_tr.state_count_dict[2] for ann_tr in batch]
            root_ages = [ann_tr.root_age for ann_tr in batch]
            
            mean_n1 = statistics.mean(n1s)
            mean_n2 = statistics.mean(n2s)
            mean_n3 = statistics.mean(n3s)
            mean_root_ages = statistics.mean(root_ages)
            
            stdevs_n1 = statistics.stdev(n1s)
            stdevs_n2 = statistics.stdev(n2s)
            stdevs_n3 = statistics.stdev(n3s)
            stdevs_root_ages = statistics.stdev(root_ages)
            
            sterr_n1 = stdevs_n1 / math.sqrt(n_sim)
            sterr_n2 = stdevs_n2 / math.sqrt(n_sim)
            sterr_n3 = stdevs_n3 / math.sqrt(n_sim)
            sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)
            
            n1_ci_width_maxtaxa = 1.96 * sterr_n1
            n2_ci_width_maxtaxa = 1.96 * sterr_n2
            n3_ci_width_maxtaxa = 1.96 * sterr_n3
            root_ages_ci_width_maxtaxa = 1.96 * sterr_root_ages

            if abs(mean_n1 - n1_mean_maxtaxa_divtree[i]) <= (n1_ci_width_maxtaxa + n1_ci_width_maxtaxa_divtree[i]):
                n1_ci_overlap_count += 1

            if abs(mean_n2 - n2_mean_maxtaxa_divtree[i]) <= (n2_ci_width_maxtaxa + n2_ci_width_maxtaxa_divtree[i]):
                n2_ci_overlap_count += 1

            if abs(mean_n3 - n3_mean_maxtaxa_divtree[i]) <= (n3_ci_width_maxtaxa + n3_ci_width_maxtaxa_divtree[i]):
                n3_ci_overlap_count += 1

            if abs(mean_root_ages - root_ages_mean_maxtaxa_divtree[i]) <= (root_ages_ci_width_maxtaxa + root_ages_ci_width_maxtaxa_divtree[i]):
                root_age_ci_overlap_count += 1

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time
    
        print("\n95% CIs of simulations here and from diversitree overlapped " + str(n1_ci_overlap_count) + " times for state 1 count.")
        print("\n95% CIs of simulations here and from diversitree overlapped " + str(n2_ci_overlap_count) + " times for state 2 count.")
        print("\n95% CIs of simulations here and from diversitree overlapped " + str(n3_ci_overlap_count) + " times for state 3 count.")
        print("\n95% CIs of simulations here and from diversitree overlapped " + str(root_age_ci_overlap_count) + " times for root age.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(n3_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        self.assertAlmostEqual(n2_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        self.assertAlmostEqual(n1_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        self.assertAlmostEqual(root_age_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)


    def test_tree_size_state_count_max_t_classe(self):
        """
        Test if ClaSSE (GeoSSE) trees simulated here have similar root ages and number of tips for the three states 
        (A, B, AB) as ClaSSE trees simulated with diversitree

        Note: condition_on_speciation=False to match diversitree!
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
            sse_sim = distsse.DnSSE(
                self.sse_stash,
                stop_condition_value,
                n=n_sim,
                stop=stop_condition,
                origin=start_at_origin,
                start_states_list=start_states_list,
                epsilon=1e-12,
                runtime_limit=3600,
                condition_on_speciation=False,
                condition_on_survival=True,
                debug=False)

            trs = sse_sim.generate()

            sim_batches.append(trs)

            # printing progress
            pjh.print_progress(i , n_batches)


        # "expectations" from diversitree
        n1_ci_width_maxt_divtree = [
            1.17669107640575, 1.5045647827377, 1.39848715432359, 1.31061633369244, 1.56348450829446, 1.46516662630169, 1.30855996168817, 1.7236330993022, 1.7506806692417, 1.36519524633148, 1.30283608215681, 1.60375338742488, 1.3919972971064, 1.63952706854289, 1.21629797252967, 1.98914275592134, 1.32397084468494, 1.47891636605245, 1.61526794900125, 1.37040250431703, 1.59596278531838, 1.57340543029469, 1.4452272487066, 1.74412414116598, 1.81470102565829, 1.40304828545319, 1.63137929174529, 1.69758238121456, 1.53606105255348, 1.4842762921758, 1.38907695592971, 1.62535505988909, 1.72964706903875, 1.48738410410066, 1.23482042631814, 1.43827792652859, 1.59315327548605, 1.42983418350634, 1.52433130893318, 1.8662475529492, 1.57484505904549, 1.41449992962499, 1.55595616982749, 1.41051240263239, 1.53001658405219, 1.32693504579447, 1.74649532539202, 1.32392541521436, 1.66040057233369, 1.26046702007108, 1.5462244529789, 1.81217712693708, 1.54353056086468, 1.51510214130462, 1.59073159155265, 1.50111516742007, 1.42016311824531, 1.77718348546516, 1.11735281605926, 1.15014166821731, 1.50524679677185, 1.4085633055914, 1.55524150243912, 1.42972155288639, 1.55286065067154, 1.34773224628513, 1.43393305790244, 1.35091433548474, 1.46442487569887, 1.2590316009224, 1.54753767063041, 1.34168385021934, 1.31645819231801, 1.22271119511791, 1.53001658405219, 1.32596966371441, 1.23682998121174, 1.45113688624687, 1.66152313522529, 1.48317132189812, 1.44931873428013, 1.62528224196925, 1.64183424156267, 1.52979465183, 1.29763579208069, 1.35795497998894, 1.94642203983918, 1.47145640984047, 1.52915530706705, 1.7495498831229, 1.46042486913358, 1.43199280569267, 1.66183138601832, 1.47624816232491, 1.35171119920662, 1.54960244833273, 1.72147162789764, 1.05936562759479, 1.53665838372842, 1.59447408181374
        ]

        n1_mean_maxt_divtree = [
            8.09, 8.27, 8.83, 7.44, 9.62, 8.09, 7.82, 9.59, 9.42, 8.3, 8.76, 9.24, 7.84, 8.74, 8.34, 9.21, 9.37, 8.57, 9.32, 8.73, 9, 9.32, 8.56, 10.13, 9.21, 8.36, 8.88, 8.93, 9.07, 9.16, 9.07, 10, 9.77, 8.26, 7.84, 8.1, 9.47, 9.21, 8.6, 9.38, 9.84, 8.59, 8.64, 8.22, 10.25, 8.38, 9.44, 7.9, 9.75, 7.58, 9.26, 10.5, 8.61, 9.27, 10.36, 9.51, 9.38, 8.87, 7.19, 6.99, 8.1, 8.1, 8.63, 8.32, 10.24, 8.03, 8.46, 8.64, 8.71, 7.36, 9.23, 8.49, 8.09, 7.65, 8.55, 8.52, 7.76, 8.75, 9.58, 8.49, 7.78, 9.19, 8.95, 8.5, 8.81, 8.59, 9.13, 9.11, 8.98, 9.28, 8.66, 7.93, 9.7, 8.91, 7.21, 8.59, 8.3, 7.33, 8.74, 10.32
        ]

        n2_ci_width_maxt_divtree = [
            1.31797828096327, 1.38161640015299, 1.56018009631238, 1.34284167256448, 1.71301028012034, 1.72060470037411, 1.45108474144989, 1.69735949704894, 1.84463118713895, 1.65820816034366, 1.25115853285802, 1.72923085605576, 1.59647207974776, 1.50538986446714, 1.5708977264834, 2.16505117545134, 1.67365007698802, 1.49538618369431, 2.0012011177129, 1.64573404039571, 1.72572103507754, 1.71487583388902, 1.66759355753472, 1.37901740746845, 1.75024617904817, 1.62097899786368, 1.66215708893078, 1.79089185433945, 1.86613111105205, 1.68720149431662, 1.42354721890229, 1.85597270094866, 1.84922713159915, 1.43535577457157, 1.44478013062839, 1.47465044870253, 1.59749261056447, 1.78559717243902, 1.78507879643868, 1.61093193230038, 1.59140471891783, 1.3729556267889, 1.29409035701577, 1.84124124317145, 1.68611443876329, 1.50833712706705, 1.76103991821178, 1.34052213316116, 1.42154911716945, 1.5781372059514, 1.84653083014248, 1.95092239298091, 1.60978493776609, 1.57714351735771, 1.64138039516017, 1.55025964276956, 1.44064478044879, 1.6662981076129, 1.54438507732712, 1.45676286541791, 1.5503935510806, 1.31768382713109, 1.92800668258199, 1.4345282821622, 1.67844372361224, 1.68092025973749, 1.62493840181619, 1.61181451372876, 1.54618304403346, 1.49538618369431, 1.74501385181617, 1.73279963065555, 1.37301639120043, 1.38778155880149, 1.6243173949681, 1.73561337173969, 1.23110671012596, 1.8863990147443, 1.53869614960462, 1.6114918797139, 1.74764474072907, 1.78848515145351, 1.54231960476171, 1.98147306126661, 1.5611845821023, 1.30963299662867, 1.98507802904927, 1.51571669269862, 1.87572050684476, 1.75816093079839, 1.77970686225444, 1.53801130731024, 1.87274742007403, 1.73221729292228, 1.85127402966965, 1.67856393841686, 1.75850961506118, 1.16684907664906, 1.4707771952202, 1.62018403992643
        ]

        n2_mean_maxt_divtree = [
             8.57, 8.26, 9.52, 8.7, 10.67, 9.37, 8.58, 11.12, 10.54, 9, 9.17, 10.4, 8.91, 8.83, 9.16, 10.61, 11.21, 10.25, 10.79, 9.89, 9.25, 10.79, 9.66, 8.95, 9.34, 9.19, 10.39, 9.92, 10.34, 10.02, 9.42, 10.7, 10.79, 8.92, 8.87, 8.86, 9.79, 10.88, 9.11, 9.23, 11.12, 8.68, 9.27, 9.56, 11.93, 9.3, 10.33, 8.48, 9.23, 8.41, 9.53, 11.57, 9.09, 10.33, 10.53, 9.84, 9.12, 9.37, 9.21, 8.97, 8.93, 8.93, 9.84, 9.26, 10.8, 9.84, 8.57, 9.64, 9.53, 8.25, 10.13, 10.04, 8.91, 9.26, 9.37, 10.51, 8.04, 10.66, 10.69, 9.58, 9.01, 9.78, 10.28, 10.67, 9.36, 9, 10.48, 10.57, 9.53, 9.4, 10.66, 9.02, 10.91, 9.56, 9.33, 10.64, 9.22, 8.25, 9.44, 11.65
        ]

        n3_ci_width_maxt_divtree = [
            0.435088072356702, 0.468879698124296, 0.488159290114664, 0.443440497414788, 0.566403040281745, 0.545093035891835, 0.399432963743921, 0.492937136566176, 0.5181929752476, 0.422945322540143, 0.420516221476952, 0.589446772299024, 0.517724742590992, 0.47399526741873, 0.505686292342027, 0.6081649377555, 0.499494202350556, 0.481919144901115, 0.519300066620369, 0.431829476258461, 0.53602298156672, 0.519927367110741, 0.497033241170329, 0.499261087948374, 0.653365999183374, 0.458604294561116, 0.472096268217242, 0.543827982310063, 0.525520918669843, 0.511773144187062, 0.437023576554683, 0.532008261498851, 0.531395222432647, 0.51744359962017, 0.417297270200958, 0.407225916501771, 0.441686894905558, 0.536862077730445, 0.585854604639904, 0.511602514584607, 0.551884764724805, 0.443440497414788, 0.479630975992286, 0.588, 0.558555713289628, 0.48580067242989, 0.529669077554252, 0.403531008909495, 0.51850739004893, 0.419264001591244, 0.549817239056431, 0.661889428954426, 0.501738190813907, 0.473737319704505, 0.597086359928043, 0.471203606509839, 0.43186541861786, 0.658288782535538, 0.416622554755876, 0.432493926389304, 0.516974688241918, 0.458705819329632, 0.563187664867214, 0.535033911509024, 0.479529835391191, 0.536135178190099, 0.463397970638542, 0.490379973885434, 0.453580358721784, 0.476269775007672, 0.460609477658776, 0.386396720817589, 0.500964201628351, 0.380975895205028, 0.514717967297224, 0.525373219528727, 0.353129812168874, 0.489461198718044, 0.556174692416734, 0.525077696607605, 0.487204469966723, 0.479335585614467, 0.453580358721784, 0.539054891603458, 0.485900507752712, 0.519221600770876, 0.644058469019041, 0.476758373886298, 0.543185421026892, 0.536963259281472, 0.55218350924198, 0.492531561415001, 0.572210286979579, 0.447660123601183, 0.455445801807238, 0.493979076152794, 0.509169592909312, 0.423866378792547, 0.522280463375993, 0.501386173875522
        ]

        n3_mean_maxt_divtree = [
            2.04, 2.12, 2.67, 2.25, 2.85, 2.27, 2.22, 2.59, 2.4, 2.49, 2.27, 2.81, 2.25, 2.51, 2.7, 2.78, 2.52, 2.57, 2.52, 2.12, 2.66, 2.56, 2.56, 2.58, 2.83, 2.4, 2.58, 2.28, 2.77, 2.48, 2.59, 2.81, 2.77, 2.2, 2.18, 1.92, 2.35, 2.82, 2.57, 2.43, 2.97, 2.45, 2.54, 2.5, 3, 2.59, 2.51, 1.94, 2.46, 2.1, 2.36, 3.1, 2.35, 2.42, 2.75, 2.59, 2.44, 2.85, 2.37, 2.14, 2.55, 2.24, 2.81, 2.77, 2.79, 2.55, 2.31, 2.23, 2.41, 2.12, 2.55, 2.18, 2.65, 2.14, 2.45, 2.63, 2.08, 2.31, 2.78, 2.43, 2.27, 2.33, 2.41, 2.54, 2.34, 2.55, 2.49, 2.32, 2.42, 2.64, 2.68, 2.22, 2.89, 2.34, 2.12, 2.46, 2.17, 2.1, 2.52, 2.96
        ]

        root_ages_ci_width_maxt_divtree = [
             0.0890150077020263, 0.100518538833934, 0.0874104901556335, 0.0869892588888224, 0.0938425765733717, 0.0868870971728242, 0.0997608042705438, 0.087853544594149, 0.102064901359859, 0.086851205328808, 0.075871793272718, 0.0813377391313974, 0.0926452862813964, 0.101199671593645, 0.0791886507707239, 0.107318411129113, 0.0745963752336609, 0.0901036209996327, 0.101116590243158, 0.0941207762367366, 0.0932618540448074, 0.0834792446607382, 0.105953800329746, 0.0886848805030738, 0.11023192406295, 0.096426369488969, 0.0996202102355397, 0.0800661606556725, 0.0995556703938077, 0.10418251020336, 0.0937627371585178, 0.0932883565934166, 0.0855686964547382, 0.0961577720432865, 0.0903651131338468, 0.113037603906392, 0.084308362365931, 0.0930913596185125, 0.106376131559568, 0.0959829834832216, 0.0853503011026563, 0.0940251386135392, 0.0896651518454209, 0.103780370719503, 0.0940381225787545, 0.100332873568977, 0.0853784325218767, 0.101012884225996, 0.0989205797003703, 0.105354708885621, 0.0990832036351655, 0.0934319132489726, 0.104307625075867, 0.083590375221642, 0.0894372067048722, 0.0890972129416705, 0.0926121036044804, 0.104531890959801, 0.0908701573284517, 0.0936802099138899, 0.0994553168518311, 0.0827393914400919, 0.0839525060976674, 0.0922856358708771, 0.0991299605436072, 0.090857355087559, 0.0968351931811871, 0.105560551761973, 0.0888630570752517, 0.0930000545144569, 0.0933621479799875, 0.0854634777621203, 0.0896437949079858, 0.0820277163674571, 0.106850057729717, 0.0903119768796217, 0.119035788606466, 0.0921128788396269, 0.0856496922485218, 0.123353641714617, 0.0912349458265095, 0.103776011720994, 0.0957357648932792, 0.0990316039317644, 0.103616427465083, 0.0960027476643842, 0.10351435335682, 0.0922342756437979, 0.0878019137425582, 0.10211010796414, 0.0921520529249162, 0.0751719530715145, 0.095069845372495, 0.0769698552138577, 0.0980111042916881, 0.0798281355669776, 0.0996935742695426, 0.0836533593668759, 0.0967188400811394, 0.111828594034809 
        ]

        root_ages_mean_maxt_divtree = [
            2.43945327182789, 2.42825587887783, 2.44113374202337, 2.46143845618943, 2.43182290624716, 2.45982299013164, 2.45606087973221, 2.52309738651583, 2.43904079163098, 2.51104711467526, 2.51610166330576, 2.495135812705, 2.41097982425691, 2.37443020658649, 2.47515420027852, 2.35062097590499, 2.55912987191545, 2.41814323973488, 2.42881679039035, 2.46732694382876, 2.48622759842651, 2.5014218394532, 2.38412906769585, 2.40271831464729, 2.40467132674051, 2.38328459466145, 2.42611434532794, 2.4584221412084, 2.45137660784115, 2.43821051696408, 2.40280152687628, 2.38278384879461, 2.49163869716708, 2.41946515131631, 2.43810471173564, 2.34739225709733, 2.47624159318175, 2.4964470472715, 2.3653257324884, 2.38401125012935, 2.51877756911575, 2.43939687760124, 2.41644348610733, 2.37421125876658, 2.45504202752404, 2.47168248828113, 2.39391106020836, 2.35462044113925, 2.40118372390009, 2.39866390359313, 2.41714657953075, 2.39163078279305, 2.33935472524114, 2.46512897105643, 2.47136769691569, 2.43960040643729, 2.4784200249298, 2.36425885781525, 2.4691915183115, 2.44685210941404, 2.37499435273031, 2.43621011360441, 2.40726862859139, 2.38361047510563, 2.45701656308728, 2.41187392796192, 2.37761184832487, 2.40378441797082, 2.42668727078288, 2.42983013376845, 2.4313042370319, 2.4338831793927, 2.44616982265902, 2.45824758346389, 2.34622080519487, 2.43557379923108, 2.33034120658501, 2.42600710393765, 2.47813864445272, 2.27210054562014, 2.37184122020502, 2.40058075286463, 2.41585791471561, 2.38326051053278, 2.44091582399657, 2.43708730530628, 2.37006703008916, 2.45793405902487, 2.40171921617121, 2.42471117261185, 2.41727595685722, 2.523797591338, 2.47581474331037, 2.50260824215688, 2.39570662222914, 2.47562003406199, 2.45598607126005, 2.43429990395648, 2.44862810145643, 2.47312647921581
        ]

        
        # parsing simulations
        n1_ci_overlap_count = 0
        n2_ci_overlap_count = 0
        n3_ci_overlap_count = 0
        root_age_ci_overlap_count = 0
        for i, batch in enumerate(sim_batches):
            n1s = [ann_tr.state_count_dict[0] for ann_tr in batch]
            n2s = [ann_tr.state_count_dict[1] for ann_tr in batch]
            n3s = [ann_tr.state_count_dict[2] for ann_tr in batch]
            root_ages = [ann_tr.root_age for ann_tr in batch]
            
            mean_n1 = statistics.mean(n1s)
            mean_n2 = statistics.mean(n2s)
            mean_n3 = statistics.mean(n3s)
            mean_root_ages = statistics.mean(root_ages)

            # debugging
            # print("batch " + str(i) + ", tree = ")
            # print(batch[0].tree.as_string(schema="newick", suppress_annotations=False))
            # print("PJ's vs. diversitree mean state 0 count: " + str(mean_n1) + " <-> " + str(n1_mean_maxt_divtree[i]))
            # print("PJ's vs. diversitree mean state 1 count: " + str(mean_n2) + " <-> " + str(n2_mean_maxt_divtree[i]))
            # print("PJ's vs. diversitree mean state 2 count: " + str(mean_n3) + " <-> " + str(n3_mean_maxt_divtree[i]))
            # print("PJ's vs. diversitree mean root age: " + str(mean_root_ages) + " <-> " + str(root_ages_mean_maxt_divtree[i]))
            
            stdevs_n1 = statistics.stdev(n1s)
            stdevs_n2 = statistics.stdev(n2s)
            stdevs_n3 = statistics.stdev(n3s)
            stdevs_root_ages = statistics.stdev(root_ages)
            
            sterr_n1 = stdevs_n1 / math.sqrt(n_sim)
            sterr_n2 = stdevs_n2 / math.sqrt(n_sim)
            sterr_n3 = stdevs_n3 / math.sqrt(n_sim)
            sterr_root_ages = stdevs_root_ages / math.sqrt(n_sim)
            
            n1_ci_width_maxt = 1.96 * sterr_n1
            n2_ci_width_maxt = 1.96 * sterr_n2
            n3_ci_width_maxt = 1.96 * sterr_n3
            root_ages_ci_width_maxtaxa = 1.96 * sterr_root_ages

            if abs(mean_n1 - n1_mean_maxt_divtree[i]) <= (n1_ci_width_maxt + n1_ci_width_maxt_divtree[i]):
                n1_ci_overlap_count += 1

            if abs(mean_n2 - n2_mean_maxt_divtree[i]) <= (n2_ci_width_maxt + n2_ci_width_maxt_divtree[i]):
                n2_ci_overlap_count += 1

            if abs(mean_n3 - n3_mean_maxt_divtree[i]) <= (n3_ci_width_maxt + n3_ci_width_maxt_divtree[i]):
                n3_ci_overlap_count += 1

            if abs(mean_root_ages - root_ages_mean_maxt_divtree[i]) <= (root_ages_ci_width_maxtaxa + root_ages_ci_width_maxt_divtree[i]):
                root_age_ci_overlap_count += 1

        # [==== * ====][.... + ....] if we take '+' to be the "truth" of the '*' interval, + cannot be more than '====' away from '*' 95% of the time
        # then abs('+' - '*') can be at most ('====' + '....'). '....' can be added because we still are guaranteed to see '+' falling within that range
        # 95% of the time
    
        print("\n95% CIs of simulations here and from diversitree overlapped " + str(n1_ci_overlap_count) + " times for state 1 count.")
        print("\n95% CIs of simulations here and from diversitree overlapped " + str(n2_ci_overlap_count) + " times for state 2 count.")
        print("\n95% CIs of simulations here and from diversitree overlapped " + str(n3_ci_overlap_count) + " times for state 3 count.")
        print("\n95% CIs of simulations here and from diversitree overlapped " + str(root_age_ci_overlap_count) + " times for root age.")
        exp_count = int(0.95 * n_batches)
        a_delta = math.ceil(0.07 * exp_count)
        self.assertAlmostEqual(n3_ci_overlap_count, exp_count,
                                msg="Mean absolute difference must be 1.96 * (stderr_python + stderr_divtree) apart " + str(exp_count) + " (+/- " + str(a_delta) + ") out of 100 times.", delta=a_delta)
        self.assertAlmostEqual(n2_ci_overlap_count, exp_count,
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
    # $ python3.9 tests/distribution/test_dn_discrete_sse_classe.py
    # 
    # or
    #
    # $ python3.9 -m tests.distribution.test_dn_discrete_sse_classe
    #
    # or 
    #
    # $ python3.9 -m unittest tests.distribution.test_dn_discrete_sse_classe.TestClaSSETrees.test_tree_size_state_count_max_t_classe

    unittest.main()