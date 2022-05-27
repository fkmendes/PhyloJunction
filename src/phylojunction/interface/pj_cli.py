import argparse

# pj imports
import phylojunction.interface.cmd.cmd_parse as cmd
import phylojunction.readwrite.pj_write as pjw

def call_cli():

    parser = argparse.ArgumentParser()
    parser.add_argument("model", action="store", type=cmd.script2pgm, help="Path to phylojunction script specifying a model")
    parser.add_argument("-d", "--data-output", dest="write_data", action="store_true", default=False, help="Toggle data output")
    parser.add_argument("-i", "--inference-output", dest="write_inference", action="store_true", default=False, help="Toggle inference script output")
    parser.add_argument("-o", "--output-dir", dest="out_dir", type=str, default="./", help="Path to project root directory, where automatic subdirectories will be created")
    parser.add_argument("-p", "--prefix", dest="prefix", type=str, default="", help="Prefix to be used when naming output files")
    
    args = parser.parse_args()

    ###########################
    # Sorting out directories #
    ###########################
    prefix = args.prefix
    root_dir = args.out_dir
    data_dir = "simulated_tables/"
    inference_dir = "inference_files/"
    if not root_dir.endswith("/"): root_dir + "/"

    #################
    # Reading model #
    #################
    pgm_obj = cmd.script2pgm(args.model)

    # debugging (looking at model)
    for node_pgm_name, node_pgm in pgm_obj.node_name_val_dict.items():
        print("\nnode name = " + node_pgm_name)
        print(node_pgm.value)

    ################
    # Writing data #
    ################
    if args.write_data:
        pjw.dump_pgm_data(data_dir, pgm_obj, prefix)

# call GUI
if __name__ == "__main__":
    call_cli()