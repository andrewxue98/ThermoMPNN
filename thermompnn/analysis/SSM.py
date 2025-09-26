import pandas as pd
from tqdm import tqdm
import torch
from omegaconf import OmegaConf

from .thermompnn_benchmarking import (
    ALPHABET,
)


def get_ssm_mutations(pdb):
    # make mutation list for SSM run
    mutation_list = []
    for seq_pos in range(len(pdb["seq"])):
        wtAA = pdb["seq"][seq_pos]
        # check for missing residues
        if wtAA != "-":
            # add each mutation option
            for mutAA in ALPHABET[:-1]:
                mutation_list.append(wtAA + str(seq_pos) + mutAA)
        else:
            mutation_list.append(None)

    return mutation_list
