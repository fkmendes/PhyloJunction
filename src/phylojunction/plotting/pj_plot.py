import matplotlib.pyplot as plt # type: ignore
import seaborn as sns # type: ignore
import pandas as pd # type: ignore
import typing as ty
from tabulate import tabulate

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.cmd.cmd_parse as cmd
import phylojunction.readwrite.pj_write as pjwrite
import phylojunction.readwrite.pj_read as pjread
import phylojunction.plotting.pj_draw as pjdraw

def join_dataframes(pj_df: pd.DataFrame, compare_to_df: pd.DataFrame, value_to_compare: str, summaries_avg_over_repl: bool=False) -> pd.DataFrame:
    """_summary_

    Args:
        pj_df (pd.DataFrame): _description_
        compare_to_df (pd.DataFrame): _description_
        value_to_compare (str): _description_
        summaries_avg_over_repl (bool, optional): _description_. Defaults to False.

    Returns:
        pd.DataFrame: _description_
    """

    # debugging
    # print(tabulate(compare_to_df, compare_to_df.head(), tablefmt="plain", showindex=False).lstrip())
    # print(tabulate(pj_df, pj_df.head(), tablefmt="plain", showindex=False).lstrip())

    sub_pj_df = pd.DataFrame()
    sub_compare_to_df = pd.DataFrame()

    # TODO:
    if summaries_avg_over_repl:
        pass
    #     pj_df will be program, sample, value_to_compare, make sure value to compare is float
    #     compare_to_df will be program, sample, value_to_compare, make sure value to compare is float
    #
    else:
        type_dict = { "program": str, "sample": int, "replicate": int, locals()["value_to_compare"]: float }
        
        try:
            sub_pj_df = pj_df.loc[:,["program", "sample", "replicate", value_to_compare]]
            sub_pj_df = sub_pj_df.astype(type_dict)
            # print(tabulate(sub_pj_df, sub_pj_df.head(), tablefmt="pretty", showindex=False).lstrip())
        
        except:
            print("Could not build joint dataframe from PJ in order to compare " + value_to_compare + ". For now, ignoring. Later implement Exception.")
            return pd.DataFrame()

        try:
            sub_compare_to_df = compare_to_df.loc[:,["program", "sample", "replicate", value_to_compare]]
            sub_compare_to_df = sub_compare_to_df.astype(type_dict)
            # print(tabulate(sub_compare_to_df, sub_compare_to_df.head(), tablefmt="pretty", showindex=False).lstrip())
        
        except:
            print("Could not build sub-dataframe from other simulator data in order to compare " + value_to_compare + ". For now, ignoring. Later implement Exception.")
            return pd.DataFrame()
    
    joint_df = pd.concat([sub_pj_df, sub_compare_to_df], axis=0)

    # debugging
    # print(tabulate(joint_df, joint_df.head(), tablefmt="pretty", showindex=False).lstrip())

    return joint_df

def draw_compare_fig(scalar_output_stash: ty.List[ty.Union[pd.DataFrame, ty.Dict[int, pd.DataFrame]]], tree_output_stash: ty.List[ty.Dict[str, pd.DataFrame]], pgm_obj: pgm.ProbabilisticGraphicalModel, ax: plt.Axes, csv_fp: str, summarize_replicates: bool=False) -> None:

    # first load the dataframe containing the data we want to compare against
    compare_to_df = pjread.read_csv_into_dataframe(csv_fp)
        
    # now we prepare PJ's dataframes, if they have not already been prepared
    output_df_list: ty.List[pd.DataFrame]
    if not scalar_output_stash or not tree_output_stash:
        scalar_output_stash, tree_output_stash = pjwrite.prep_data_df(pgm_obj) # still prefixless
        _, output_df_list = pjwrite.prep_data_filepaths_dfs(scalar_output_stash, tree_output_stash)

    


if __name__ == "__main__":
    # initializing figure
    comparison_fig = plt.figure() # note that pjgui uses matplotlib.figure.Figure (which is part of Matplotlib's OOP class library)
                              # here, we instead use pyplot's figure, which is the Matlab-like state-machine API
    comparison_fig_axes = comparison_fig.add_axes([0.25, 0.2, 0.5, 0.6])
    comparison_fig_axes.patch.set_alpha(0.0)
    comparison_fig_axes.xaxis.set_ticks([])
    # comparison_fig_axes.yaxis.set_ticks([])
    comparison_fig_axes.spines['left'].set_visible(True)
    comparison_fig_axes.spines['bottom'].set_visible(True)
    comparison_fig_axes.spines['right'].set_visible(False)
    comparison_fig_axes.spines['top'].set_visible(False)

    # initializing model
    model_fp = "examples/multiple_scalar_tree_plated.pj"
    pgm_obj = cmd.script2pgm(model_fp, in_pj_file=True)

    # getting pj dataframes
    scalar_output_stash: ty.List[ty.Union[pd.DataFrame, ty.Dict[int, pd.DataFrame]]]
    tree_output_stash: ty.List[ty.Dict[str, pd.DataFrame]]
    output_df_list: ty.List[pd.DataFrame]
    scalar_output_stash, tree_output_stash = pjwrite.prep_data_df(pgm_obj) # still prefixless
    _, output_df_list = pjwrite.prep_data_filepaths_dfs(scalar_output_stash, tree_output_stash)
    
    for df in output_df_list:
        df["program"] = "PJ"

    # getting dataframe to compare to
    # csv_fp = "examples/compare_files/comparison_multiple_scalar_tree_plated_rv5.csv"
    csv_fp = "examples/compare_files/comparison_multiple_scalar_tree_plated_trs2.csv"
    compare_to_df = pjread.read_csv_into_dataframe(csv_fp)
    
    # try to join dataframes
    # output_df_list[2] has data for "rv5", output_df_list[8] has data for "trs2"
    # joint_df = join_dataframes(output_df_list[2], compare_to_df, value_to_compare="rv5", summaries_avg_over_repl=False)
    joint_df = join_dataframes(compare_to_df, output_df_list[8], value_to_compare="root_age", summaries_avg_over_repl=False)

    # trying some plotting
    # sns.violinplot(x="program", y="rv5", data=joint_df, ax=comparison_fig_axes)
    # sns.violinplot(x="program", y="root_age", data=joint_df, ax=comparison_fig_axes)
    # pjdraw.plot_violins(comparison_fig, comparison_fig_axes, joint_df, "program", "rv5", xlab="Program", ylab="rv5", color1="blue", color2="red")
    pjdraw.plot_violins(comparison_fig, comparison_fig_axes, joint_df, "program", "root_age", xlab="Program", ylab="Root age")
    plt.show()
