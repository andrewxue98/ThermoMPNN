import sys

sys.path.append("../")

import numpy as np
import os
import pandas as pd
from tqdm import tqdm
from Bio.Align import PairwiseAligner

from protein_mpnn_utils import parse_PDB


def get_global_alignment(seq1, seq2):
    """Perform global sequence alignment using Bio.Align"""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    return next(aligner.align(seq1, seq2))


def load_csv_data(csv_loc):
    """Load CSV data and key colun unique IDs"""
    df = pd.read_csv(csv_loc)
    protein_names = df["protein_name"].unique()
    pdb_ids = df["pdb_id_corrected"].unique()
    sequences = df["sequence"]
    return df, protein_names, pdb_ids, sequences


def get_df_slice(pdb_id, df):
    pdb_id = pdb_id.split("_")[0]
    data_idx = df["pdb_id_corrected"].str.contains(pdb_id)
    df_slice = {
        "wt_AA": df["wild_type"][data_idx],
        "mut_AA": df["mutation"][data_idx],
        "pos": df["position"][data_idx],
        "seq": df["sequence"][data_idx].unique()[0],
    }
    return df_slice, data_idx


def seq1_index_to_seq2_index(alignment, index):
    """Convert index from seq1 to seq2 after alignment"""
    cur_seq1_index = 0

    # first find the aligned index
    for aln_idx, (char1, char2) in enumerate(zip(alignment.target, alignment.query)):
        if char1 != "-":
            cur_seq1_index += 1
        if cur_seq1_index > index:
            break

    # now get the index in seq2 corresponding to aligned index
    if alignment.query[aln_idx] == "-":
        return None

    seq2_to_idx = alignment.query[: aln_idx + 1]
    seq2_idx = aln_idx
    for char in seq2_to_idx:
        if char == "-":
            seq2_idx -= 1

    if seq2_idx < 0:
        return None

    return seq2_idx


class CustomAlignment:

    def __init__(self, csv_loc, out_loc, pdb_loc):

        # load CSV data
        self.df, self.protein_names, self.pdb_ids, self.sequences = load_csv_data(
            csv_loc
        )

        # set up new columns
        pdb_position = [[] for _ in range(self.df.shape[0])]
        pdb_sequence = [[] for _ in range(self.df.shape[0])]

        # tracker variables
        self.successful = [False for _ in range(self.df.shape[0])]
        mutations = 0  # recording # of checks (can be more than # muts)

        used_pdbs = sorted(os.listdir(pdb_loc))
        for pdb in tqdm(self.pdb_ids):
            pdb_fname = pdb + ".pdb"
            if pdb_fname not in used_pdbs:
                continue
            # print(pdb)
            # retrieve PDB data
            pdb_data = parse_PDB(os.path.join(pdb_loc, pdb_fname))
            pdb_seq = pdb_data[0]["seq"].replace("-", "X")
            # retrieve CSV data
            df_slice, data_idx = get_df_slice(pdb, self.df)
            alignment = get_global_alignment(df_slice["seq"], pdb_seq)
            iterrows = data_idx[
                data_idx == True
            ].index  # iterate through rows for this seq and check mutations

            # iterate over mutations and search through each one
            for idx in iterrows:
                csv_pos = df_slice["pos"][idx] - 1
                if (
                    df_slice["seq"][csv_pos] != df_slice["wt_AA"][idx]
                ):  # CSV should be internally consistent
                    continue

                new_position = seq1_index_to_seq2_index(align, csv_pos)
                if new_position is not None:  # happens on hetero-oligomers, mostly
                    mutations += 1
                    # check if offset is correct and log if so
                    if pdb_seq[new_position] == df_slice["seq"][csv_pos]:
                        self.successful[idx] = True
                        # add successful chain/position to list
                        pdb_position[idx] = new_position
                        # add position in concatenated sequence and increment by chain length

            # collect all results for a given protein and update the CSV to include the PDB-derived position and sequence
            for itt in iterrows:
                if self.successful[itt]:
                    pdb_sequence[itt] = pdb_seq.replace("-", "X")
                    pdb_position[itt] = pdb_position[itt]

            print("Total Mutations Processed:", mutations)
            print("Successful Mutation Alignments:", sum(self.successful))

        # update DF and save changes
        self.df["pdb_position"] = pdb_position
        self.df["pdb_sequence"] = pdb_sequence
        # drop those cases that failed to align (if any)
        self.df = self.df.loc[self.successful, :].reset_index(drop=True)
        print("Dataset size after alignment and quality check:", self.df.shape)
        self.df.to_csv(out_loc)


csv_loc = "../data/fireprot/1_fireprotDB_cleaned.csv"
out_loc = "../data/fireprot/2_fireprotDB_aligned.csv"
pdb_loc = "../data/fireprot/pdbs/monomers"

a = CustomAlignment(csv_loc, out_loc, pdb_loc)
