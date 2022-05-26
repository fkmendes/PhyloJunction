import argparse

# pj imports
import phylojunction.interface.cmd.cmd_parse as cmd

def call_cli():

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--model', dest="pgm_obj", type=cmd.script2pgm, help="Path to phylojunction script specifying a model")
    
    args = parser.parse_args()

    print(args.pgm_obj)

# call GUI
if __name__ == "__main__":
    call_cli()