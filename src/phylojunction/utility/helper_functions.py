import sys
sys.path.extend(["../", "../phylojunction"]) # necessary to run it as standalone on command line (from phylojunction/ or phylojunction/utility/)
import time
import typing as ty
import numpy as np # type: ignore
from collections import defaultdict

# pj imports
import utility.exception_classes as ec

def print_progress(idx: int, iterator_len: int) -> None:
    """Print progress bar on iteration for longer tests

    Args:
        idx (int): Index to keep track of where in iteration one is
        iterator_len (int): Size of the iterator whose progress is keeping track of
    """

    j = (idx + 1) / iterator_len
    sys.stdout.write('\r')
    sys.stdout.write("[%-20s] %d%%" % ('='*int(20*j), 100*j))
    sys.stdout.flush()
    time.sleep(0.25)


def autovivify(levels=1, final=dict) -> ty.DefaultDict:
    return (defaultdict(final) if levels < 2 else
            defaultdict(lambda: autovivify(levels - 1, final)))


def verify_or_convert2_vector(param_list, dn_name, size_to_grow=1) -> ty.List[ty.List[ty.Union[int, float, str]]]:
    """

    Args:
        size_to_grow (int): This will be the number of simulations (n_draws inside a DistributionPGM object)
    """
    param_list_array = np.array(param_list, dtype=object) # using numpy array as a hack to get nested-ness of param_list
    n_params = len(param_list_array.shape) # each entry in shape is a dimension

    vectorizable_param_list = param_list
    vectorized_param_list: ty.List[ty.List[ty.Union[int, float, str]]] = list()

    one_scalar_provided: bool = False
    one_scalar_in_list_provided: bool = False
    one_par_right_value_count: bool = False

    # scalar provided by itself, e.g., rate = 1.0
    if isinstance(param_list, (int, float, str)):
        vectorizable_param_list = [vectorizable_param_list]
        one_scalar_provided = True

    # dealing with weird cases -- must set some flags
    # and/or initialize variables appropriately
    elif type(param_list) == list:
        # single scalar provided by itself inside list, e.g, param_list = [1.0]
        if len(param_list) == 1:
            one_scalar_in_list_provided = True
        # we have a single parameter inside a list of size_to_grow already,
        # but no nestedness; this happens when we do
        # rv ~ someDn(n=5, param=[1, 2, 3, 4, 5])
        #
        # here, we want to just paste those values as-is, inside vectorized_param_list,
        # but we do need to initialize this second dimension already so the
        # nestedness is correct
        elif len(param_list) == size_to_grow and n_params == 1:
            one_par_right_value_count = True
            vectorized_param_list.append([])

    # each v here will be a different parameter (e.g., mean and sd) in case there
    # is more than one parameter
    for v in vectorizable_param_list:
        # a single value was provided as scalar, either by itself, or inside a list by itself
        if isinstance(v, (int, float, str)):
            # param list is a single scalar, e.g., mean=1.0
            if one_scalar_provided:
                vectorized_param_list.append([v for i in range(size_to_grow)])
            # param list is a single scalar, by itself, inside a list, e.g., mean=[[1.0]]
            elif one_scalar_in_list_provided:
                vectorized_param_list.append([v] * size_to_grow)
            # param list contains a single parameter, and we already have the right number of values
            elif one_par_right_value_count:
                vectorized_param_list[0].append(v)

        elif type(v) == list:
            n_val = len(v)
            # more values than specified number of samples
            if n_val > size_to_grow:
                raise ec.DimensionalityError(dn_name)

            # don't know how to multiply if more than one element
            elif n_val > 1 and n_val < size_to_grow:
                raise ec.DimensionalityError(dn_name)

            # a single value was provided as list
            elif n_val == 1:
                vectorized_param_list.append([v[0] for i in range(size_to_grow)])

            # more than 1 value, but same as specified number of samples
            else:
                vectorized_param_list.append(v)

    return vectorized_param_list


def get_ellapsed_time_in_minutes(start: float, end: float) -> int:
        """Calculate ellapsed time

        Args:
            start (float): Start of time window
            end (float): End of time window

        Returns:
            int: Ellapsed time in minutes in time window
        """
        ellapsed_minutes, ellapsed_secs = divmod(end-start, 60) # returns (min, sec)

        # return int(ellapsed_minutes * 60 + ellapsed_secs) # in seconds
        return int(ellapsed_minutes)

##############################################################################

if __name__ == "__main__":
    # can be called from utility/
    # $ python3 helper_functions.py
    # 
    # can also be called from phylojunction/
    # $ python3 utility/helper_functions.py
    # or
    # $ python3 -m utility.helper_functions
    #
    # can also be called from VS Code, if open folder is phylojuction/
    
   d = autovivify(levels=3, final=dict)
   d[0][1][2] = 3
   print(d)

   print_progress(49, 100)

   print(get_ellapsed_time_in_minutes(0.0, 125)) # 2 minutes (and 5 seconds)