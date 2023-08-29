import typing as ty
import numpy as np
import enum
import random
import math
from abc import ABC, abstractmethod
from copy import deepcopy

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.utility.helper_functions as pjh

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

# MacroevolEvent = enum.Enum("Macroevol. event", ["W_SPECIATION", "BW_SPECIATION", "ASYM_SPECIATION", "EXTINCTION", "ANAGENETIC_TRANSITION", "ANCESTOR_SAMPLING"], start=0)

class MacroevolEvent(enum.Enum):
    W_SPECIATION = 0
    BW_SPECIATION = 1
    ASYM_SPECIATION = 2
    EXTINCTION = 3
    ANAGENETIC_TRANSITION = 4
    ANCESTOR_SAMPLING = 5 

##############################################################################


class DiscreteStateDependentParameter():
    """
    Main class for discrete state-dependent parameters

    Supports vectorization of values only.
    """

    value: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]]
    name: str
    state: int=0  # associated state

    def __init__(self,
        val: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]],
        name: str="", state: int=0):

        self.name = name
        self.state = state
        self.value = []

        # scalar
        if isinstance(val, (float, int, str)):
            # vectorizing if input wasn't in list or ndarray form
            self.value = [float(val)]

        # list
        elif isinstance(val, list):
            if isinstance(val[0], (int, float, str, np.float64)):
                self.value = [float(v) for v in val]
            
            else:
                raise ec.StateDependentParameterMisspec(
                    message="Could not recognize type of value argument, " \
                        " for parameter " + self.name + ". Cannot initialize " + name + ".")

        else:
            raise ec.StateDependentParameterMisspec(
                message="Argument to value parameter is not either in scalar or " \
                    + "vectorized form (it is likely an object). " \
                    + " Cannot initialize. " + name)

    @abstractmethod
    def _initialize_str_representation(self) -> None:
        pass

    def __len__(self):
        if isinstance(self.value, list):
            return len(self.value)
        else:
            return 1

##############################################################################


class DiscreteStateDependentRate_old_and_working():
    """
    Main class for discrete state-dependent rate parameters

    Supports vectorization of values only.
    """
    
    value: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]]
    event: MacroevolEvent
    name: str
    states: ty.List[int]

    def __init__(self,
        val: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]],
        event: MacroevolEvent, name: str="", states: ty.List[int]=[0,0,0]):

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
                raise ec.StateDependentParameterMisspec(
                    message="Could not recognize type of value argument, " \
                        " for parameter " + self.name + ". Cannot initialize " + name + ".")

        else:
            raise ec.StateDependentParameterMisspec(
                message="Argument to value parameter is not either in scalar or " \
                    + "vectorized form (it is likely an object). " \
                    + " Cannot initialize. " + name)

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
        if self.event == MacroevolEvent.EXTINCTION and self.arriving_states:
            raise ec.StateDependentParameterMisspec(message="Extinction only takes a " \
                + "departing state (and no arriving states).")

        if self.event == MacroevolEvent.W_SPECIATION and \
            len(set(self.arriving_states)) != 1 and \
                not self.departing_state in self.arriving_states:
            raise ec.StateDependentParameterMisspec(message="Within-state speciation " \
                + "requires the departing and arriving states to be all " \
                + "the same.")

        if self.event == MacroevolEvent.BW_SPECIATION and \
            len(set(self.arriving_states)) != 2 and \
                self.departing_state in self.arriving_states:
            raise ec.StateDependentParameterMisspec(message="Between-state speciation " \
                + "requires the departing and the two arriving states to be " \
                + "all different.")

        if self.event == MacroevolEvent.ASYM_SPECIATION and \
            len(set(self.arriving_states)) != 2 and \
                not self.departing_state in self.arriving_states:
            raise ec.StateDependentParameterMisspec(message="Asymmetric-state " \
                + "speciation requires the departing and one of the " \
                + "arriving states to be the same, and the other " \
                + "arriving state to be different. Exiting...")

        if self.event == MacroevolEvent.ANAGENETIC_TRANSITION and \
            len(set(self.arriving_states)) != 1 and \
                self.departing_state in self.arriving_states:
            raise ec.StateDependentParameterMisspec(message="State transition " \
                + "requires the departing and arriving states to be " \
                + "different.")

        if self.event == MacroevolEvent.ANCESTOR_SAMPLING and self.arriving_states:
            raise ec.StateDependentParameterMisspec(message="Ancestor sampling " \
                + "only takes a departing state (and no arriving " \
                + "states).")

    def __str__(self):
        return self.str_representation

    def __repr__(self):
        return self.str_representation

    def __len__(self) -> int:
        if isinstance(self.value, list):
            return len(self.value)
        else:
            return 1

##############################################################################


class DiscreteStateDependentRate(DiscreteStateDependentParameter):
    """
    Main class for discrete state-dependent rate parameters

    Supports vectorization of values only.
    """
    
    state_tuple: ty.Tuple[int]
    
    def __init__(self,
        val: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]],
        event: MacroevolEvent, name: str="",
        states: ty.List[int]=[]) -> None:

        # side-effect: initializes self.value, self.name, self.state
        self.states: ty.List[int] = []
        if len(states) == 0:
            self.states = [0, 0, 0]

        else:
            self.states = states

        state = int(self.states[0])
        super().__init__(val=val, name=name, state=state)

        self.state_tuple = tuple(int(s) for s in self.states)
        self.departing_state = self.state_tuple[0]
        self.arriving_states = self.state_tuple[1:] # does not include time_slice_index 1
        
        # making sure inputs are ok
        self.event = event
        if self.event == MacroevolEvent.EXTINCTION and self.arriving_states \
            and (len(states) > 0):
            raise ec.StateDependentParameterMisspec(message="Extinction only takes a " \
                + "departing state (and no arriving states).")

        if self.event == MacroevolEvent.W_SPECIATION and \
            len(set(self.arriving_states)) != 1 and \
                not self.departing_state in self.arriving_states:
            raise ec.StateDependentParameterMisspec(message="Within-state speciation " \
                + "requires the departing and arriving states to be all " \
                + "the same.")

        if self.event == MacroevolEvent.BW_SPECIATION and \
            len(set(self.arriving_states)) != 2 and \
                self.departing_state in self.arriving_states:
            raise ec.StateDependentParameterMisspec(message="Between-state speciation " \
                + "requires the departing and the two arriving states to be " \
                + "all different.")

        if self.event == MacroevolEvent.ASYM_SPECIATION and \
            len(set(self.arriving_states)) != 2 and \
                not self.departing_state in self.arriving_states:
            raise ec.StateDependentParameterMisspec(message="Asymmetric-state " \
                + "speciation requires the departing and one of the " \
                + "arriving states to be the same, and the other " \
                + "arriving state to be different. Exiting...")

        if self.event == MacroevolEvent.ANAGENETIC_TRANSITION and \
            len(set(self.arriving_states)) != 1 and \
                self.departing_state in self.arriving_states:
            raise ec.StateDependentParameterMisspec(message="State transition " \
                + "requires the departing and arriving states to be " \
                + "different.")

        if self.event == MacroevolEvent.ANCESTOR_SAMPLING and self.arriving_states:
            raise ec.StateDependentParameterMisspec(message="Ancestor sampling " \
                + "only takes a departing state (and no arriving " \
                + "states).")

        self._initialize_str_representation()


    def _initialize_str_representation(self):
        self.str_representation = "Discrete state-dependent rate\n" \
            + "   Name:                " + self.name + "\n" \
            + "   Value:               " \
            + ", ".join(str(v) for v in self.value) + "\n" \
            + "   Departing state:     " + str(self.departing_state) + "\n" \
            + "   Arriving state(s):   " \
            + ", ".join(str(v) for v in self.arriving_states) + "\n" \
            + "   Associated event:    " + str(self.event) + "\n\n"


    def __str__(self):
        return self.str_representation


    def __repr__(self):
        return self.str_representation


    def __len__(self) -> int:
        return super().__len__()

##############################################################################


class DiscreteStateDependentProbability(DiscreteStateDependentParameter):
    """
    Main class for discrete state-dependent probability parameters

    Supports vectorization of values only.
    """

    state_representation: str

    def __init__(self,
        val: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]],
        name: str="", state: int=0):
        
        # side-effect: populates self.values (list of floats)
        super().__init__(val=val, name=name, state=state)

        for v in self.value:
            if v < 0.0:
                raise ec.NotBetweenZeroAndOneError(name, "negative")

            elif v > 1.0:
                raise ec.NotBetweenZeroAndOneError(name, "positive")

        self._initialize_str_representation()


    def _initialize_str_representation(self):
        self.str_representation = "Discrete state-dependent probability\n" \
            + "   Name:   " + self.name + "\n" \
            + "   Value:  " \
            + ", ".join(str(v) for v in self.value) + "\n" \
            + "   State:  " + str(self.state) + "\n\n"


    def __str__(self):
        return self.str_representation


    def __repr__(self):
        return self.str_representation


    def __len__(self) -> int:
        return super().__len__()

##############################################################################


class DiscreteStateDependentParameterManager:
    """Stash for discrete state-dependent parameters and time slices

    At the moment, this class does not care about vectorization.
    It manipulates whole instances of DiscreteStateDependentParameter,
    which then in turn contain multiple values if vectors have been passed
    by user.

    Later, vectorization of seed ages for time slicing and time slices themselves
    might become vectorized.
    """

    matrix_state_dep_params: \
        ty.List[ty.List[DiscreteStateDependentParameter]]
    seed_age: ty.Optional[float]
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.List[ty.Optional[float]]
    epsilon: float

    # NOTE: This class is flexible in that it allows different parameter
    # numbers per time slice, for whatever that is worth. However, the
    # user interface has a check inside make_SSEStash()
    # that forces the user to specify the same number of parameters in
    # all time slices

    def __init__(self,
                 matrix_state_dep_params: \
                    ty.List[ty.List[DiscreteStateDependentParameter]],
                 total_state_count: int,
                 seed_age_for_time_slicing: ty.Optional[float] = None,
                 list_time_slice_age_ends: ty.Optional[ty.List[float]] = None,
                 epsilon: float=1e-12):

        # 1D: time slices
        # 2D: list of atomic rate params
        self.matrix_state_dep_params = matrix_state_dep_params 

        self._check_single_parameter_type()
        
        self.state_count = total_state_count

        # this is the origin or root age, and it is used to anchor
        # the user-specified ages to convert it to time
        self.seed_age = seed_age_for_time_slicing 
        
        # default is one slice, but this gets updated below
        self.n_time_slices = 1
        
        # default slice ends at present (age = 0.0)
        self.slice_age_ends = [ 0.0 ]

        # default slice t's at present
        # (t = seed_age = stop_value with "age" condition)
        self.slice_t_ends = [ self.seed_age ]

        # age ends (larger in the past, 0.0 in the present)
        if list_time_slice_age_ends:
            self.n_time_slices += len(list_time_slice_age_ends)

            # we do not want to change this object outside
            # of this class instance, because other state-dep
            # param managers may be created (and will be taking
            # 'list_time_slice_age_ends' as arguments)
            self.slice_age_ends = deepcopy(list_time_slice_age_ends)
            
            # appends age end of time slice ending in present
            self.slice_age_ends.append(0.0)

            # convert into time ends
            # (0.0 at origin, seed_age = stop_val for "age" stop
            # condition in the present)
            if self.seed_age:
                self.slice_t_ends = \
                    [self.seed_age - age_end for \
                        age_end in self.slice_age_ends]
                        # no need to append seed_age, because
                        # self.slice_age_ends already has 0.0 in it

        if self.n_time_slices > 1:
            self._check_all_states_in_all_time_slices()

        # TODO: at this point, we should check that we have
        # a # of parameters that is a multiple of the number of
        # states
        # 
        # if manager is holding rates, we need to build a dictionary
        # of state-triplets [(0,0,0), (2,0,1)] and make sure all triplets
        # appear in all time slices ()
        #
        # some of this testing is already being done at the interface level
        # within dn_fn_discrete_sse.py.make_SSEStash(), but
        # that testing is weaker, because it just checks that the number
        # of specified parameters is a multiple of the number of time slices
        # irrespective of what parameter states are

        self.state_dep_params_dict: \
            ty.Dict[int, ty.List[ty.List[DiscreteStateDependentParameter]]] = \
                dict((s, [[] for j in range(self.n_time_slices)]) \
                for s in range(self.state_count))
        
        # side effect: initializes self.atomic_rate_params_dict
        self._init_matrix_state_dep_params_dict()
        # self.state_dep_params_dict =
        # { state0:
        #           [ [ #slice1#; atomic_param1, atomic_param2, ...] [ #slice2#; atomic_param1, atomic_param2 ] ], ... ]
        #   state1:
        #           [ [ #slice1#; atomic_param1, atomic_param2, ...] [ #slice2#; atomic_param1, atomic_param2 ] ], ... ]
        # }

        self.epsilon = epsilon


    def _check_single_parameter_type(self):
        different_par_type_set: ty.Set(DiscreteStateDependentParameter) \
            = set()
        
        # k-th time slice
        for k, list_params in enumerate(self.matrix_state_dep_params):
            for param in list_params:
                this_param_type = type(param)
                
                if not this_param_type in different_par_type_set:
                    different_par_type_set.add(this_param_type)

        if len(different_par_type_set) >= 2:
            raise(ec.RequireSameParameterType(
                "DiscreteStateDependentParameterManager",
                len(different_par_type_set))
            )


    def _check_all_states_in_all_time_slices(self):
        states_per_epoch_dict = \
            dict((k, set()) for k in range(self.n_time_slices))
        
        # k-th time slice
        for k, list_params in enumerate(self.matrix_state_dep_params):
            for param in list_params:
                sts = tuple()

                # try states (if rate)
                try:
                    sts = param.state_tuple

                # if not, it's a probability, then state
                except:
                    sts = param.state

                if sts in states_per_epoch_dict[k]:
                    raise ec.RepeatedStateDependentParameterError(
                        k, sts, ""
                    )

                states_per_epoch_dict[k].add(sts)

        # debugging
        # print("states_per_epoch_dict after populating")
        # print(states_per_epoch_dict)

        oldest_epoch_state_set = states_per_epoch_dict[0]
        for k, state_tup_set in states_per_epoch_dict.items():
            if k > 0:
                if oldest_epoch_state_set != state_tup_set:
                    raise ec.MissingStateDependentParameterError(
                        k,
                        pjh.symmetric_difference(
                            oldest_epoch_state_set,
                            state_tup_set
                        ),
                        ""
                    )


    def _init_matrix_state_dep_params_dict(self):

        # k-th time slice
        for k, list_params in enumerate(self.matrix_state_dep_params):

            if (k - self.n_time_slices) > -1:
                raise ec.ObjInitIncorrectDimensionError(
                    "DiscreteStateDependentParameterManager",
                    self.n_time_slices,
                    len(self.matrix_state_dep_params),
                    message=str(self.n_time_slices) + " time slices was" \
                        + (" (were) specified, but from the matrix of"
                           " state-dependent parameters, it looks like"
                           " there are ") \
                        + str(len(self.matrix_state_dep_params)) + " slices"
                )

            for param in list_params:
                try:
                    self.state_dep_params_dict[param.state][k].append(param)
                
                except:
                    exit("Parameter " + param.name + "'s associated " \
                        + "(possibly departing) state is " + str(param.state) +
                        ", but this state is not represented between 0 and " +
                        str(self.state_count) +
                        ". Likely the total number of states was misspecified.")


    def state_dep_params_at_time(self,
        a_time: float,
        params_matrix: ty.Optional[
            ty.List[ty.List[DiscreteStateDependentParameter]]
            ]=None) -> ty.List[DiscreteStateDependentParameter]:
        # see where a_time falls within (get time_slice_index)
        # grab list of atomic rate params at that time_slice_index
        # print("self.n_time_slices = " + str(self.n_time_slices))

        time_slice_index = -1
        # try:
        while time_slice_index < self.n_time_slices:
            time_slice_index += 1

            try:
                if isinstance(self.slice_t_ends, list) and \
                    isinstance(self.slice_t_ends[time_slice_index], float):
                        time_slice_t_end = \
                            ty.cast(float, self.slice_t_ends[time_slice_index])
                        
                        if a_time > time_slice_t_end or \
                            (abs(a_time - time_slice_t_end) <= self.epsilon):
                            continue

            # self.slice_t_ends will be None if no seed age
            # or time slice age ends are provided
            except:
                # in which case we want the index to just be 0
                time_slice_index = 0

            break

        # adjusting
        # if time_slice_index > (self.n_time_slices - 1):
        #     time_slice_index = self.n_time_slices - 1

        # params_matrix will have been provided as a
        # matrix conditioned on a departing state
        if params_matrix != None:
            return params_matrix[time_slice_index]
        
        # when no params_matrix is specified, just grab
        # all parameters at time t, but warning:
        # no control for the order (in different time slices)
        # the parameters have been passed by user!!!!
        else:
            return self.matrix_state_dep_params[time_slice_index]


    def __len__(self) -> int:
        if isinstance(self.matrix_state_dep_params, list):
            return sum([ len(list_atomic_rates_in_time_slice) for list_atomic_rates_in_time_slice in self.matrix_state_dep_params])
        else:
            return 0

##############################################################################


class DiscreteStateDependentParameterManager_old_working:
    """Stash for discrete state-dependent parameters and time slices

    At the moment, this class does not care about vectorization. It manipulates
    whole instances of DiscreteStateDependentParameter, which then in turn
    contain multiple values if vectors have been passed by user.

    Later, vectorization of seed ages for time slicing and time slices themselves
    might become vectorized.
    """

    seed_age: ty.Optional[float]
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.List[ty.Optional[float]]

    # NOTE: This class is flexible in that it allows different parameter
    # numbers per time slice, for whatever that is worth. However, the
    # user interface has a check inside make_SSEStash()
    # that forces the user to specify the same number of parameters in
    # all time slices

    def __init__(self,
                 matrix_atomic_rate_params: ty.List[ty.List[DiscreteStateDependentRate]],
                 total_state_count: int,
                 seed_age_for_time_slicing: ty.Optional[float]=None,
                 list_time_slice_age_ends: ty.Optional[ty.List[float]]=None,
                 epsilon: float=1e-12):
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

        self.atomic_rate_params_dict: ty.Dict[int, ty.List[ty.List[DiscreteStateDependentRate]]] = dict((s, [[] for j in range(self.n_time_slices)]) for s in range(self.state_count))
        
        self.init_atomic_rate_param_dict(matrix_atomic_rate_params) # side effect: initializes self.atomic_rate_params_dict
        # self.atomic_rate_params_dict =
        # { state0:
        #           [ [ #slice1#; atomic_param1, atomic_param2, ...] [ #slice2#; atomic_param1, atomic_param2 ] ], ... ]
        #   state1:
        #           [ [ #slice1#; atomic_param1, atomic_param2, ...] [ #slice2#; atomic_param1, atomic_param2 ] ], ... ]
        # }

        self.epsilon = epsilon


    def __len__(self) -> int:
        if isinstance(self.atomic_rate_params_matrix, list):
            return sum([ len(list_atomic_rates_in_time_slice) for list_atomic_rates_in_time_slice in self.atomic_rate_params_matrix])
        else:
            return 0


    def init_atomic_rate_param_dict(self, matrix_atomic_rate_params: ty.List[ty.List[DiscreteStateDependentRate]]):
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


class MacroevolEventHandler():
    """Class for sampling of discrete-state macroevolutionary events

    Some of the methods in this class are vector-aware (through parameter
    value_idx), because instances of this class will be directly called by
    dn_sse.simulate() method, which in turn will try to access parameter values
    stored in a #-simulations-sized list.
    """

    # NOTE: This class depends on DiscreteStateDependentParameterManager, which allows different
    # parameter numbers per time slice, for whatever that is worth. However,
    # the user interface has a check inside make_SSEStash()
    # that forces the user to specify the same number of parameters in
    # all time slices

    state_dep_rate_manager: DiscreteStateDependentParameterManager
    state_count: int
    n_time_slices: int
    seed_age: ty.Optional[float]
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.List[ty.Optional[float]]
    str_representation: str

    def __init__(self,
        state_dep_rate_manager: DiscreteStateDependentParameterManager) \
            -> None:

        self.state_dep_rate_manager = state_dep_rate_manager
        self.state_count = self.state_dep_rate_manager.state_count
        self.n_time_slices = self.state_dep_rate_manager.n_time_slices
        self.seed_age = self.state_dep_rate_manager.seed_age
        self.slice_t_ends = self.state_dep_rate_manager.slice_t_ends
        self.slice_age_ends = self.state_dep_rate_manager.slice_age_ends

        self.str_representation = "State-dependent rates"
        # state s
        for s, atomic_rates_state_mat in self.state_dep_rate_manager.state_dep_params_dict.items():
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
    def total_rate(self,
        a_time: float,
        state_representation_dict: ty.Dict[int, ty.Set[str]],
        value_idx: int=0,
        departing_state: ty.Optional[int]=None,
        debug: bool=False) -> \
            ty.Union[float, ty.Tuple[float, ty.List[float]]]:
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
            atomic_rate_params_matrix = self.state_dep_rate_manager.state_dep_params_dict[state_idx] # conditioning
            atomic_rate_params_at_time = self.state_dep_rate_manager.state_dep_params_at_time(a_time, params_matrix=atomic_rate_params_matrix)

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
    def sample_event_atomic_parameter(self,
        denominator: float,
        a_time: float,
        state_indices: ty.List[int],
        value_idx: int=0,
        a_seed: ty.Optional[float]=None,
        debug: bool=False):
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
            atomic_rate_params_matrix = self.state_dep_rate_manager.state_dep_params_dict[state_idx]
            this_state_atomic_rate_params = self.state_dep_rate_manager.state_dep_params_at_time(a_time, params_matrix=atomic_rate_params_matrix)
            all_states_atomic_rate_params += this_state_atomic_rate_params

            # total rate of outcomes must depend on "adjacent" states across events
            for atomic_rate_param in this_state_atomic_rate_params:
                ws.append( (atomic_rate_param.value[value_idx] / denominator) )

        if debug:
            print("    \nThe following events are allowed:")
            print("   " + "   ".join(str(ap) for ap in all_states_atomic_rate_params))

        return random.choices( [ atomic_param for atomic_param in all_states_atomic_rate_params ], weights=ws)


    def __len__(self) -> int:
        if self.state_dep_rate_manager:
            return len(self.state_dep_rate_manager)

        else:
            return 0

    def __str__(self) -> str:
        return self.str_representation

    def __repr__(self) -> str:
        return self.str_representation

    def get_gcf(self):
        pass

    def get_length(self):
        pass

##############################################################################


class DiscreteStateDependentProbabilityHandler():
    """Class for managing time-heterogeneous state-dependent probabilities
    
    Time slices do not have to be the same as those for state-dependent
    rates.
    """

    # NOTE: This class depends on DiscreteStateDependentParameterManager,
    # which allows different parameter numbers per time slice, for whatever
    # that is worth.
    # 
    # However, the user interface has a check inside
    # make_SSEStash() that forces the user to specify the same
    # number of parameters in all time slices

    state_dep_prob_manager: DiscreteStateDependentParameterManager
    state_count: int
    n_time_slices: int
    seed_age: ty.Optional[float]
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.List[ty.Optional[float]]
    str_representation: str

    def __init__(self,
        state_dep_prob_manager: DiscreteStateDependentParameterManager) \
            -> None:

        self.state_dep_prob_manager = state_dep_prob_manager
        self.state_count = self.state_dep_prob_manager.state_count
        self.n_time_slices = self.state_dep_prob_manager.n_time_slices
        self.seed_age = self.state_dep_prob_manager.seed_age
        self.slice_t_ends = self.state_dep_prob_manager.slice_t_ends
        self.slice_age_ends = self.state_dep_prob_manager.slice_age_ends

        # side-effect: initializes self.str_representation
        self._initialize_str_representation()


    def _initialize_str_representation(self) -> None:
        self.str_representation = "State-dependent probabilities"
        # state s
        for s, atomic_rates_state_mat in \
            self.state_dep_prob_manager.state_dep_params_dict.items():
            
            self.str_representation += "\n  State " + str(s) + ":\n"

            # time slice k
            for k, list_state_dep_params_slice in enumerate(atomic_rates_state_mat):
                if self.seed_age and isinstance(self.slice_t_ends, list):
                    time_slice_t_end = ty.cast(float, self.slice_t_ends[k])
                    self.str_representation += "    Time slice " \
                        + str(k + 1) + " (time = " + str(round(time_slice_t_end,4)) \
                        + ", age = " + str(round(self.slice_age_ends[k],4)) + ")\n"
                
                else:
                    self.str_representation += "    Time slice " \
                        + str(k + 1) + " (age = " \
                        + str(round(self.slice_age_ends[k],4)) + ")\n"

                for state_dep_param in list_state_dep_params_slice:
                    self.str_representation += "      " + state_dep_param.name + " = "
                    self.str_representation += ", ".join(
                        str(v) for v in state_dep_param.value) + "\n"


    def _state_dep_prob_at_time(self, a_time, state_idx) \
        -> ty.List[DiscreteStateDependentParameter]:

        # scoped to total_rate
        state_cond_prob_matrix = self.state_dep_prob_manager.state_dep_params_dict[state_idx] # conditioning

        # in the context of rates (multiple rates per departing state)
        # state_dep_params_at_time returns a list of rates,
        # but here, we have a single probability inside a list
        state_cond_prob_at_time = self.state_dep_prob_manager.state_dep_params_at_time(a_time, params_matrix=state_cond_prob_matrix)
        
        return state_cond_prob_at_time


    def randomly_decide_taxon_sampling_at_time_at_state(
        self, a_time, state_idx, sample_idx) \
        -> bool:
        
        prob_at_state_in_slice_list = \
            self._state_dep_prob_at_time(a_time, state_idx)

        if len(prob_at_state_in_slice_list) > 1: 
            exit("Should only have one probability deterministic node " \
                + " per state per slice. Exiting...")

        prob_at_state_in_slice_ith_sample = \
            prob_at_state_in_slice_list[0].value[sample_idx]

        if prob_at_state_in_slice_ith_sample == 1.0:
            return True

        elif prob_at_state_in_slice_ith_sample == 0.0:
            return False
        
        elif np.random.binomial(1, prob_at_state_in_slice_ith_sample):
            return True
        
        else:
            return False


    def __len__(self) -> int:
        if self.state_dep_prob_manager:
            return len(self.state_dep_prob_manager)

        else:
            return 0


    def __str__(self) -> str:
        return self.str_representation

##############################################################################


class SSEStash():
    """
    Class for stashing state-dependent rate and probability handlers
    """

    meh: MacroevolEventHandler
    prob_handler: DiscreteStateDependentProbabilityHandler
    str_representation: str

    def __init__(self,
        macroevol_event_handler: MacroevolEventHandler,
        state_dep_prob_handler: ty.Optional[DiscreteStateDependentProbabilityHandler]=None) \
            -> None:

        self.meh = macroevol_event_handler
        self.prob_handler = state_dep_prob_handler

        if state_dep_prob_handler == None:
            self._initialize_missing_prob_handler()

        self._init_str_representation()


    # side-effect: creates and populates self.prob_handler if user did
    # not provide it
    def _initialize_missing_prob_handler(self) -> None:
        expected_n_prob \
            = self.meh.state_count * self.meh.state_dep_rate_manager.n_time_slices
            # have to use the rate manager number of slices, as that
            # does not reflect that last 0.0
        n_prob_per_slice = int(expected_n_prob / self.meh.n_time_slices)
        
        matrix_state_dep_probs: \
            ty.List[ty.List[DiscreteStateDependentProbability]] = []
        
        # TODO: this code is wrong, fix it!
        for i in range(0, expected_n_prob, n_prob_per_slice):

            state_dep_probs: \
                ty.List[DiscreteStateDependentProbability] = []
            for j in range(self.meh.state_count):
                st = j % n_prob_per_slice
                state_dep_prob = DiscreteStateDependentProbability(
                    name="rho" + str(st) + "_t" + str(j // n_prob_per_slice),
                    val=1.0, state=st)
                    
                # so mypy won't complain
                if isinstance(state_dep_prob,
                    DiscreteStateDependentProbability):
                    state_dep_probs.append(state_dep_prob)

            # appending a list of DiscreteStateDependentProbability
            # 1D: time slices, 2D: state-dep prob list
            matrix_state_dep_probs.append(state_dep_probs)

        # we have to throw away the last element in slice_age_ends
        # because it has already been added a 0.0 for the youngest
        # time slice when 'meh' was initialized; if we don't do this
        # the number of time slices inside 'state_dep_probs_manager'
        # will be 1 too large
        state_dep_probs_manager = \
            DiscreteStateDependentParameterManager(
                matrix_state_dep_probs,
                self.meh.state_count,
                seed_age_for_time_slicing=self.meh.seed_age,
                list_time_slice_age_ends=self.meh.slice_age_ends[:-1]
            )

        self.prob_handler = DiscreteStateDependentProbabilityHandler(
            state_dep_probs_manager
        )


    def _init_str_representation(self) -> None:
        self.str_representation = str(self.meh)

        if self.prob_handler:
            self.str_representation += "\n" + str(self.prob_handler)


    def get_meh(self) -> MacroevolEventHandler:
        return self.meh


    def get_prob_handler(self) -> DiscreteStateDependentProbabilityHandler:
        return self.prob_handler


    def __str__(self):
        return self.str_representation

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
    # arp = DiscreteStateDependentRate(1.0, MacroevolEvent.W_SPECIATION, name="lambda")
    # print(arp)
    # rates_t0_s0 = [ arp ]
    # matrix_atomic_rate_params = [ rates_t0_s0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice        
    # state_dep_par_manager = DiscreteStateDependentParameterManager(matrix_atomic_rate_params, 1)

    # n_characters = 3 # regions A and B
    # n_states_per_char = 2 # presence/absence
    # svc = StateIntoPatternConverter(n_characters, n_states_per_char)

    # print(svc.int2set_dict)
    # print(svc.set2int_dict)

    n_characters = 8 # regions A and B
    n_states_per_char = 2 # presence/absence
    svc = StateIntoPatternConverter(n_characters, n_states_per_char)

    # print(svc.int2set_dict)
    for k, v in svc.int2set_dict.items():
        print(k + "\t" + v)
    # print(svc.set2int_dict)