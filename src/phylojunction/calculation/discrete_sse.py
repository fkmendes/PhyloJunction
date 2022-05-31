import typing as ty
import numpy as np
import enum
import random
import math

# pj imports
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

MacroevolEvent = enum.Enum("Macroevol. event", ["W_SPECIATION", "BW_SPECIATION", "ASYM_SPECIATION", "EXTINCTION", "ANAGENETIC_TRANSITION"], start=0)

##############################################################################

class AtomicSSERateParameter():
    """Main class for discrete-state macroevolutionary parameters

    Supports vectorization of values only.
    """
    
    value: ty.List[float]
    
    def __init__(self, val: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]], event: MacroevolEvent, name: str="", states: ty.List[int]=[0,0,0]):
        self.name = name
        self.value = []
        
        # scalar
        # if not type(val) in (list, np.ndarray):
        # if type(val) in (int, float, str):
        if isinstance(val, (float, int, str)):
            self.value = [float(val)] # vectorizing if input wasn't in list or ndarray form
            # if it is an instance of an object like RandomVariablePGM

        # list
        # TODO: make sure that parametric distributions get output converted to list 
        # inside their own classes so the line below is not necessary (and mypy passes)
        # or type(val) == np.ndarray:
        elif isinstance(val, list):
            # if isinstance(val, np.ndarray):
            #     val = val.tolist()
            if isinstance(val[0], (int, float, str, np.float64)):
                self.value = [float(v) for v in val]
            else:
                raise ec.SSEAtomicRateMisspec(message="Could not recognize type of value parameter. Cannot initialize. Exiting...")

        else:
            raise ec.SSEAtomicRateMisspec(message="Argument to value parameter is not already in vectorized or scalar form (it is likely an object). Cannot initialize. Exiting...")

        self.state_tuple = tuple(int(s) for s in states)
        self.departing_state = self.state_tuple[0]
        self.arriving_states = self.state_tuple[1:] # does not include time_slice_index 1
        self.event = event

        self.str_representation = "AtomicRateParam\n"
        self.str_representation += "   Name:                " + self.name + "\n"
        self.str_representation += "   Value:               " + ", ".join(str(v) for v in self.value) + "\n"
        self.str_representation += "   Departing state:     " + str(self.departing_state) + "\n"
        self.str_representation += "   Arriving state(s):   " + ", ".join(str(v) for v in self.arriving_states) + "\n"
        self.str_representation += "   Associated event:    " + str(self.event) + "\n\n"

        # making sure inputs are ok
        if self.event == MacroevolEvent.W_SPECIATION and \
            len(set(self.arriving_states)) != 1 and \
                not self.departing_state in self.arriving_states:
            raise ec.SSEAtomicRateMisspec(message="Within-state speciation requires the departing and arriving states to be all the same. Exiting...")

        if self.event == MacroevolEvent.BW_SPECIATION and \
            len(set(self.arriving_states)) != 2 and \
                self.departing_state in self.arriving_states:
            raise ec.SSEAtomicRateMisspec(message="Between-state speciation requires the departing and the two arriving states to be all different. Exiting...")

        if self.event == MacroevolEvent.ASYM_SPECIATION and \
            len(set(self.arriving_states)) != 2 and \
                not self.departing_state in self.arriving_states:
            raise ec.SSEAtomicRateMisspec(message="Asymmetric-state speciation requires the departing and one of the arriving states to be the same, and " + \
                "the other arriving state to be different. Exiting...")

        if self.event == MacroevolEvent.ANAGENETIC_TRANSITION and \
            len(set(self.arriving_states)) != 1 and \
                self.departing_state in self.arriving_states:
            raise ec.SSEAtomicRateMisspec(message="State transition requires the departing and arriving states to be different. Exiting...")

    def __str__(self):
        return self.str_representation

    def __repr__(self):
        return self.str_representation

    def sample(self):
        super().sample()

    def get_length(self):
        return super().get_length()

    def get_gcf(self):
        pass

##############################################################################

class FIGRatesManager:
    """Stash for discrete-state macroevolutionary parameters and time slices

    At the moment, this class does not care about vectorization. It manipulates
    whole instances of AtomicSSERateParameters, which _then in turn_
    contain multiple values if vectors have been passed by user.

    Later, vectorization of seed ages for time slicing and time slices themselves
    might become vectorized.
    """

    slice_t_ends: ty.List[ty.Optional[float]]

    # NOTE: This class is flexible in that it allows different parameter
    # numbers per time slice, for whatever that is worth. However, the
    # user interface has a check inside make_MacroEvolEventHandler()
    # that forces the user to specify the same number of parameters in
    # all time slices

    def __init__(self, matrix_atomic_rate_params: ty.List[ty.List[AtomicSSERateParameter]], total_state_count: int, seed_age_for_time_slicing: ty.Optional[float]=None, list_time_slice_age_ends: ty.Optional[ty.List[float]]=None, epsilon: float=1e-12):
        # how do we need to access parameters and values during simulation?
        # we probably need two types of parameter containers.
        # type 1: records "high-level" governing parameters, rho, sigma, phi
        # type 2: records "low-level" event-specific rates, which are used
        #         to simulate next event times, classes, and outcomes
        
        self.atomic_rate_params_matrix = matrix_atomic_rate_params # 1D: time slices, 2D: list of atomic rate params
        self.state_count = total_state_count
        # TODO: make this vectorized next, so a prior can be put on seed age
        self.seed_age = seed_age_for_time_slicing # this is the origin or root age, and it is used to anchor the user-specified ages to convert it to time
        self.n_time_slices = 1 # default is one slice
        self.slice_age_ends = [ 0.0 ] # default slice ends at present (age = 0.0)
        self.slice_t_ends = [ self.seed_age ] # default slice t's at present (t = seed_age = stop_value with "age" condition)

        # age ends (larger in the past, 0.0 in the present)
        if list_time_slice_age_ends:
            self.n_time_slices += len(list_time_slice_age_ends)
            self.slice_age_ends = list_time_slice_age_ends
            self.slice_age_ends.append(0.0) # appends age end of time slice ending in present

            # convert into time ends (0.0 at origin, seed_age = stop_val for "age" stop condition in the present)
            if self.seed_age:
                self.slice_t_ends = [self.seed_age - age_end for age_end in self.slice_age_ends] # no need to append seed_age, because self.slice_age_ends already has 0.0 in it

        self.atomic_rate_params_dict: ty.Dict[int, ty.List[ty.List[AtomicSSERateParameter]]] = dict((s, [[] for j in range(self.n_time_slices)]) for s in range(self.state_count))
        self.init_atomic_rate_param_dict(matrix_atomic_rate_params) # side effect: initializes self.atomic_rate_params_dict
        # self.atomic_rate_params_dict =
        # { state0:
        #           [ [ #slice1#; atomic_param1, atomic_param2, ...] [ #slice2#; atomic_param1, atomic_param2 ] ], ... ]
        #   state1:
        #           [ [ #slice1#; atomic_param1, atomic_param2, ...] [ #slice2#; atomic_param1, atomic_param2 ] ], ... ]
        # }

        self.epsilon = epsilon


    def init_atomic_rate_param_dict(self, matrix_atomic_rate_params: ty.List[ty.List[AtomicSSERateParameter]]):
        # original implementation (2nd dimension were states, rather than all rates from all states together)
        # for k, s_state_list in enumerate(matrix_atomic_rate_params):

        # k-th time slice
        for k, list_atomic_rate_params in enumerate(matrix_atomic_rate_params):

            # original implementation had this extra loop
            # for s, list_atomic_rate_params in enumerate(s_state_list):

            for atomic_rate_param in list_atomic_rate_params:
                try:
                    self.atomic_rate_params_dict[atomic_rate_param.departing_state][k].append(atomic_rate_param)
                except:
                    exit("Parameter " + atomic_rate_param.name + "'s associated departing state is " + \
                            str(atomic_rate_param.departing_state) +
                            ", but this departing state is not represented between 0 and " +
                            str(self.state_count) +
                            ". Likely the total number of states was misspecified.")


    def atomic_rate_params_at_time(self, atomic_rate_params_matrix, a_time: float):
        # see where a_time falls within (get time_slice_index)
        # grab list of atomic rate params at that time_slice_index
        # print("self.n_time_slices = " + str(self.n_time_slices))

        time_slice_index = -1
        # try:
        while time_slice_index < self.n_time_slices:
            time_slice_index += 1

            try:
                if isinstance(self.slice_t_ends, list) and isinstance(self.slice_t_ends[time_slice_index], float):
                        time_slice_t_end = ty.cast(float, self.slice_t_ends[time_slice_index])
                        
                        if a_time > time_slice_t_end or (abs(a_time - time_slice_t_end) <= self.epsilon):
                            continue

            # self.slice_t_ends will be None if no seed age or time slice age ends are provided
            except:
                time_slice_index = 0 # in which case we want the index to just be 0

            break        

        # adjusting
        # if time_slice_index > (self.n_time_slices - 1):
        #     time_slice_index = self.n_time_slices - 1

        return atomic_rate_params_matrix[time_slice_index]

##############################################################################

class MacroEvolEventHandler():
    """Class for sampling of discrete-state macroevolutionary events

    Some of the methods in this class are vector-aware (through parameter
    value_idx), because instances of this class will be directly called by
    dn_sse.simulate() method, which in turn will try to access parameter values
    stored in a #-simulations-sized list.
    """

    # NOTE: This class depends on FIGRatesManager, which allows different
    # parameter numbers per time slice, for whatever that is worth. However,
    # the user interface has a check inside make_MacroEvolEventHandler()
    # that forces the user to specify the same number of parameters in
    # all time slices

    slice_t_ends: ty.List[ty.Optional[float]]

    def __init__(self, a_fig_rates_manager: FIGRatesManager) -> None:
        self.fig_rates_manager = a_fig_rates_manager
        self.state_count = self.fig_rates_manager.state_count
        self.n_time_slices = self.fig_rates_manager.n_time_slices
        self.seed_age = self.fig_rates_manager.seed_age
        self.slice_age_ends = self.fig_rates_manager.slice_age_ends
        self.slice_t_ends = self.fig_rates_manager.slice_t_ends

        self.str_representation = "MacroEvolEventHandler"
        # state s
        for s, atomic_rates_state_mat in self.fig_rates_manager.atomic_rate_params_dict.items():
            self.str_representation += "\n  State " + str(s) + ":\n"

            # time slice k
            for k, list_atomic_rates_slice in enumerate(atomic_rates_state_mat):
                if self.seed_age and isinstance(self.slice_t_ends, list):
                    time_slice_t_end = ty.cast(float, self.slice_t_ends[k])
                    self.str_representation += "    Time slice " + str(k + 1) + " (time = " + str(round(time_slice_t_end,4)) + ", age = " + str(round(self.slice_age_ends[k],4)) + ")\n"
                else:
                    self.str_representation += "    Time slice " + str(k + 1) + " (age = " + str(round(self.slice_age_ends[k],4)) + ")\n"

                # atomic rate 'ar'
                for ar in list_atomic_rates_slice:
                    self.str_representation += "      " + ar.name + " = "
                    self.str_representation += ", ".join(str(v) for v in ar.value) + "\n"

    # this function deals with vectorization
    def total_rate(self, a_time: float, state_representation_dict: ty.Dict[int, ty.Set[str]], value_idx: int=0, departing_state: ty.Optional[int]=None, debug: bool=False) -> ty.Union[float, ty.Tuple[float, ty.List[float]]]:
        """Calculate total rate for either any event, or for events conditioned on a specific state.

        Args:
            a_time (float): Forward time (not age!) we want to recover rates with.
            state_representation_dict (dict): Dictionary with (int) states as keys and list of (str) node labels at that state.
            value_idx (int, optional): Index specifying which parameter value (within a vector of #-of-simulations-size) we care about. Defaults to 0 (first element).
            departing_state (int, optional): State we can condition rates on. Defaults to None.
            debug (bool, optional): Flag for printing debugging messages. Defaults to False.

        Returns:
            (float): If no departing_state is provided, returns the global rate. Otherwise, returns a state-conditioned total rate and a list with each state total rate.
        """

        if debug:
            if departing_state is None:
                print("\nCalculating rate for exponential:")
            else:
                print("\nCalculating denominator for sampling event:")

        total_rate = 0.0
        state_rates = [0.0 for i in range(len(state_representation_dict))]

        for state_idx, nd_labels in state_representation_dict.items():
            # if departing state is provided, then we do not care about rates not departing from it
            if departing_state is not None and departing_state != state_idx:
                continue

            n_lineages_in_state = len(nd_labels)

            if debug:
                if departing_state is None:
                    print("  state " + str(state_idx) + " is represented by " + str(n_lineages_in_state) + " lineages (weight).")

            # scoped to total_rate
            atomic_rate_params_matrix = self.fig_rates_manager.atomic_rate_params_dict[state_idx] # conditioning
            atomic_rate_params_at_time = self.fig_rates_manager.atomic_rate_params_at_time(atomic_rate_params_matrix, a_time)

            for atomic_rate_param in atomic_rate_params_at_time:
                w = 1.0 # weight

                # if we already drew the time to the next event and a node
                # representation does not matter, the total is just the
                # sum of all rates of the specified state (this is guaranteed
                # by the if statement above)

                # otherwise, we are getting the total rate for the whole tree,
                # before picking a node; the weight is then the number of
                # targetable nodes at a given state
                if departing_state is None:
                    w = float(n_lineages_in_state)

                if debug:
                    print("    This rate's value is " + str(atomic_rate_param.value[value_idx]))

                # value_idx is what introduces vectorization here
                # if we have 100 lambda values (b/c we are doing 100 simulations)
                # then we will be only looking at the i-th (0 < i < 99) lambda at a time
                state_rates[state_idx] += atomic_rate_param.value[value_idx]
                total_rate += w * atomic_rate_param.value[value_idx]

        if debug:
            print("Total rate = " + str(total_rate))

        if departing_state is None:
            return total_rate, state_rates
        else:
            return total_rate


    # this function deals with vectorization
    def sample_event_atomic_parameter(self, denominator: float, a_time: float, state_indices: ty.List[int], value_idx: int=0, a_seed: ty.Optional[float]=None, debug: bool=False):
        """Return one-sized list with a sampled macroevolutionary event

        Args:
            denominator (float): Normalization term for computing weights of different events.
            a_time (float): When the event will take place.
            state_indices (int): Event is conditioned on departing from this state(s).
            value_idx (int, optional): Index specifying which parameter value (within a vector of #-of-simulations-size) we care about. Defaults to 0 (first element).
            a_seed (int, optional): Random seed for simulation (not working at the moment). Defaults to None.
            debug (bool, optional): Flag for printing debugging messages. Defaults to False.

        Returns:
            MacroevolEvent: Sampled macroevolutionary event that will take place
        """

        if debug:
            print("    Will sample event from state_indices = " + ", ".join(str(i) for i in state_indices))

        if a_seed:
            random.seed(a_seed)

        all_states_atomic_rate_params = list()
        ws = list() # weights for sampling proportional to rate value
        for state_idx in state_indices:
            atomic_rate_params_matrix = self.fig_rates_manager.atomic_rate_params_dict[state_idx]
            this_state_atomic_rate_params = self.fig_rates_manager.atomic_rate_params_at_time(atomic_rate_params_matrix, a_time)
            all_states_atomic_rate_params += this_state_atomic_rate_params

            # total rate of outcomes must depend on "adjacent" states across events
            for atomic_rate_param in this_state_atomic_rate_params:
                ws.append( (atomic_rate_param.value[value_idx] / denominator) )

        if debug:
            print("    The following events are allowed:")
            print("\n    ".join(str(ap) for ap in all_states_atomic_rate_params))

        return random.choices( [ atomic_param for atomic_param in all_states_atomic_rate_params ], weights=ws)

    def __str__(self):
        return self.str_representation

    def __repr__(self):
        return self.str_representation

    def get_gcf(self):
        pass

    def get_length(self):
        pass

##############################################################################

class StateIntoPatternConverter:
    """Stash and machinery for checking and converting character compound-states into bit patterns,
    and vice-versa
    
    For example, if we are thinking about regions as characters, then compound-state means different ranges (e.g., "A", "AB"),
    while "state" is the number of different values a single character can take.
    In the case of 1-character-1-region, then the number of states is 2, because a species is either present or not in a region.
    """

    # e.g., for 2 states
    # Region | compound-state  | Bit pattern 
    # --------------------------------------
    # Dead   |        0        | 000
    # A      |        1        | 100
    # B      |        2        | 010
    # C      |        3        | 001
    # AB     |        4        | 110
    # AC     |        5        | 101
    # BC     |        6        | 011
    # ABC    |        7        | 111

    n_char: int
    n_states_per_char: int
    n_states: int
    int2set_dict: ty.Dict[str, str] = dict()
    set2int_dict: ty.Dict[str, str] = dict()

    def __init__(self, n_characters: int, n_states_per_char: int) -> None:
        self.n_char = n_characters
        self.n_states_per_char = n_states_per_char
        self.n_states = int(n_states_per_char ** n_characters)

        # side-effect
        self.populate_dicts()

    def populate_dicts(self) -> None:
        """Non-recursively generate all bit patterns for a given number of characters and states per characters

        Note that if we are thinking about regions as characters, then "compound state" means the range (e.g., "AB"),
        and state is the number of options per character (e.g., binary if present/absent in a region)
        """
        
        # non-recursive method
        _list_for_sorting: ty.List[ty.List[str]] = [[] for i in range(self.n_char + 1)]
        for compound_state_idx in range(self.n_states):
            _bit_pattern = [0 for i in range(self.n_char)]
            _z = compound_state_idx

            for bit_idx in range(self.n_char-1, -1, -1):
                _v = self.n_states_per_char ** bit_idx # n_states_per_char here is present/absent for regions (i.e., binary state)
                
                if _z >= _v:
                    _bit_pattern[bit_idx] = 1
                    _z -= _v
                
                else:
                    _bit_pattern[bit_idx] = 0

            _bit_pattern_str = "".join(str(b) for b in _bit_pattern)
            _n_bits_on = sum(_bit_pattern)
            _list_for_sorting[_n_bits_on].append(_bit_pattern_str)

        # will store all bit patterns after sorting first by number of bits that are on, and then by the decimal value underlying the bit pattern
        sorted_list_patterns: ty.List[str] = [] 
        for list_of_patterns in _list_for_sorting:
            sorted_list_patterns.extend(list_of_patterns)

        # getting the other dict by reversing this one
        self.int2set_dict = dict((str(idx), bit_pattern) for idx, bit_pattern in enumerate(sorted_list_patterns))
        for k, v in self.int2set_dict.items():
            self.set2int_dict[v] = k


##############################################################################

if __name__ == "__main__":
    # can be called from calculation/
    # $ python3 discrete_sse.py
    # 
    # can also be called from phylojunction/
    # $ python3 calculation/discrete_sse.py
    # or
    # $ python3 -m calculation.discrete_sse
    #
    # can also be called from VS Code, if open folder is phylojuction/
    
    # yule, one state, one epoch
    # arp = AtomicSSERateParameter(1.0, MacroevolEvent.W_SPECIATION, name="lambda")
    # print(arp)
    # rates_t0_s0 = [ arp ]
    # matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice        
    # fig_rates_manager = FIGRatesManager(matrix_atomic_rate_params, 1)

    n_characters = 3 # regions A and B
    n_states_per_char = 2 # presence/absence
    svc = StateIntoPatternConverter(n_characters, n_states_per_char)

    print(svc.int2set_dict)
    print(svc.set2int_dict)