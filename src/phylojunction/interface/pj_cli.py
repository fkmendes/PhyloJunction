import argparse
import os

# pj imports
import phylojunction.interface.cmdbox.cmd_parse as cmd
import phylojunction.readwrite.pj_write as pjw


def execute_pj_script(
        model: str,
        prefix: str="",
        root_dir: str="./",
        write_data: bool=False,
        write_inference: bool=False,
        write_nex_states: bool=False) -> None:
    """
    Execute .pj script

    This is called by application 'pjcli' from the terminal,
    but can also be called from a .py script importing phylojunction
    """

    ########################################
    # Sorting out and creating directories #
    ########################################
    root_dir = root_dir
    if not root_dir.endswith("/"):
        root_dir += "/"

    if not os.path.isdir(root_dir):
        os.mkdir(root_dir)

    #################
    # Reading model #
    #################
    pgm_obj = cmd.script2pgm(model, in_pj_file=True)

    # debugging (looking at model)
    # for node_name, node_pgm in pgm_obj.node_name_val_dict.items():
    #     print("\nnode name = " + node_name)
    #     print(node_pgm.value)

    ################
    # Writing data #
    ################
    if write_data:
        data_dir = root_dir + "simulated_tables/"

        if not os.path.isdir(data_dir):
            os.mkdir(data_dir)

        pjw.dump_pgm_data(data_dir, pgm_obj, prefix, write_nex_states)

    ###################
    # Writing figures #
    ###################

    # TODO

    ###########################
    # Writing inference files #
    ###########################
    if write_inference:
        inference_dir = root_dir + "inference_files/"

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
        help="String specifying which stochastic nodes to draw figures for, " +
             "and in what range of the sample " +
             "(e.g., \'birth_rate,tr;0-1,0-10\')")
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
        help="Path to project root directory, where automatic subdirectories " +
             "will be created")
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
        root_dir=args.out_dir,
        write_data=args.write_data,
        write_inference=args.write_inference,
        write_nex_states=args.write_nex_states
    )

# if one wants to run pj_cli.py for some reason
if __name__ == "__main__":
    call_cli()