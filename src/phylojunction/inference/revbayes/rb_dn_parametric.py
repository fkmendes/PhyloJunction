import sys
# sys.path.extend(["../../", "../phylojunction"]) # necessary to run it as standalone on command line (from its own path all the way up to phylojunction/)
import typing as ty
import math

# pj imports
import pgm.pgm as pgm

def get_rev_str_from_dn_parametric_obj(dn_obj: pgm.DistributionPGM) -> ty.Tuple[int, int, ty.List[str]]:
    
    rev_str_list: ty.List[str]
    n_sim: int = dn_obj.n_draws
    n_repl: int = dn_obj.n_repl
    rev_str_list = dn_obj.get_rev_inference_spec_info()
    
    return n_sim, n_repl, rev_str_list

def get_unif_rev_inference_spec_info(n_draws: int) -> ty.List[str]:
    pass

def get_exponential_rev_inference_spec_info(n_draws: int) -> ty.List[str]:
    pass

def get_gamma_rev_inference_spec_info(n_draws: int) -> ty.List[str]:
    pass

def get_normal_rev_inference_spec_info(n_draws: int) -> ty.List[str]:
    pass

def get_ln_rev_inference_spec_info(n_draws: int, param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]], parent_node_tracker: ty.Optional[ty.Dict[str, str]]) -> ty.List[str]:
    rev_str_list: ty.List[str] = [] # return
    real_mean_list = param_dict["mean_param"] # one per simulation
    real_sd_list = param_dict["sd_param"] # one per simulation

    # so mypy won't complain
    if isinstance(real_mean_list, list) and isinstance(real_sd_list, list) and \
        isinstance(parent_node_tracker, dict):
        real_mean_list = [ math.exp(float(i)) for i in real_mean_list ]
        real_sd_list = [ math.exp(float(i)) for i in real_sd_list ]
        
        # real_mean_list and real_sd_list will have n_sim values inside, even if they are all the same
        for ith_sim in range(n_draws):
            ith_sim_str = "dnLognormal(mean="

            # if we can find a node that holds the value of the mean, we use it
            try:
                ith_sim_str += parent_node_tracker["mean"] # returns NodePGM, and we grab its name
            except:
                ith_sim_str += str(real_mean_list[ith_sim])

            ith_sim_str += ", sd="

            # if we can find a node that holds the value of the sd, we use it
            try:
                ith_sim_str += parent_node_tracker["sd"] # returns NodePGM, and we grab its name
            except:
                ith_sim_str += str(real_sd_list[ith_sim])
            ith_sim_str += ")"

            rev_str_list.append(ith_sim_str)

    return rev_str_list