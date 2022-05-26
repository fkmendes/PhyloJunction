import argparse

# pj imports
import phylojunction.interface.cmd.cmd_parse as cmd

def call_cli():

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model', dest="pgm_obj", type=cmd.script2pgm, help="Path to phylojunction script specifying a model")
    parser.add_argument('-o', '--output-dir', dest="out_dir", type=str, default="./", help="Path to output simulated data")
    
    args = parser.parse_args()

    pgm_obj = args.pgm_obj

    for node_pgm_name, node_pgm in pgm_obj.node_name_val_dict.items():
        print("\nnode name = " + node_pgm_name)
        print(node_pgm.value)

# call GUI
if __name__ == "__main__":
    call_cli()