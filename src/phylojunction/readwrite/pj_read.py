import sys
sys.path.extend(["../", "../phylojunction"])
import os
import typing as ty
import numpy as np
import pandas as pd # type: ignore

def read_text_file(fp_string: str) -> ty.List[str]:
    cmd_line_list = list()
    
    with open(fp_string, "r") as infile:
        for line in infile:
            line = line.rstrip()

            if line: cmd_line_list.append(line)
    
    return cmd_line_list