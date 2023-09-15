import argparse
import os
import typing as ty

# pj imports
import phylojunction.interface.cmdbox.cmd_parse as cmd
import phylojunction.pgm.pgm as pgm
import phylojunction.readwrite.pj_read as pjr
import phylojunction.readwrite.pj_write as pjw
import phylojunction.interface.pjcli.cli_plotting as cliplt

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


def execute_pj_script(
        model: str,
        prefix: str = "",
        out_dir: str = "./",
        write_data: bool = False,
        write_figures: str = "",
        write_inference: bool = False,
        write_nex_states: bool = False) -> None:
    """
    Execute .pj script

    This is called by application 'pjcli' from the terminal,
    but can also be called from a .py script importing phylojunction
    """

    # Reading model #
    pgm_obj: pgm.ProbabilisticGraphicalModel = \
        pgm.ProbabilisticGraphicalModel()

    print("Reading script " + model)
    pgm_obj = cmd.script2pgm(model, in_pj_file=True)
    print("    ... done!")

    n_samples = pgm_obj.sample_size

    fig_obj, fig_axes = cliplt.start_fig_and_axes()

    # making sure the PGM has at least a note in it #
    if pgm_obj.n_nodes > 0:

        ########################################
        # Sorting out and creating directories #
        ########################################
        output_dir: str = out_dir
        if not out_dir.endswith("/"):
            output_dir += "/"

        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        # debugging (looking at model)
        # for node_name, node_pgm in pgm_obj.node_name_val_dict.items():
        #     print("\nnode name = " + node_name)
        #     print(node_pgm.value)

        # Writing data #
        if write_data:
            data_dir: str = output_dir

            if not os.path.isdir(data_dir):
                os.mkdir(data_dir)

            pjw.dump_pgm_data(
                data_dir,
                pgm_obj,
                prefix,
                write_nex_states)

        # Writing figures #
        if write_figures:
            node_range_dict: ty.Dict[str, ty.Tuple[int]] = \
                pjr.parse_cli_str_write_fig(write_figures)
            fig_dir: str = output_dir + "figures/"

            if not os.path.isdir(fig_dir):
                os.mkdir(fig_dir)

            cliplt.call_node_plot_cli(
                fig_dir,
                pgm_obj,
                n_samples,
                node_range_dict,
                fig_obj,
                fig_axes)
            
        # Writing inference files #
        if write_inference:
            inference_dir: str = out_dir + "inference_files/"

            if not os.path.isdir(inference_dir):
                os.mkdir(inference_dir)
            # TODO


# the function that pjcli application on terminal calls
def call_cli() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "script",
        action="store",
        type=str,
        help="Path to phylojunction script specifying a model")
    parser.add_argument(
        "-d",
        "--data-output",
        dest="write_data",
        action="store_true",
        default=False,
        help="Toggle data output")
    parser.add_argument(
        "-f",
        "--figure-output",
        dest="write_figures",
        type=str,
        action="store",
        help=("String specifying which stochastic nodes to draw figures "
              "for, and in what range of the sample (e.g., "
              "\'birth_rate,tr;0-1,0-10\')"))
    parser.add_argument(
        "-i",
        "--inference-output",
        dest="write_inference",
        action="store_true",
        default=False,
        help="Toggle inference script output")
    parser.add_argument(
        "-o",
        "--output-dir",
        dest="out_dir",
        type=str,
        default="./",
        help=("Path to project root directory, where automatic "
              "subdirectories will be created"))
    parser.add_argument(
        "-p",
        "--prefix",
        dest="prefix",
        type=str,
        default="",
        help="Prefix to be used when naming output files")
    parser.add_argument(
        "-nex", "--nexus-state-output",
        dest="write_nex_states",
        action="store_true",
        default=False,
        help="Toggle states nexus file output")

    args = parser.parse_args()

    execute_pj_script(
        args.script,
        prefix=args.prefix,
        out_dir=args.out_dir,
        write_data=args.write_data,
        write_figures=args.write_figures,
        write_inference=args.write_inference,
        write_nex_states=args.write_nex_states
    )

# if one wants to run pj_cli.py for some reason
if __name__ == "__main__":
    call_cli()