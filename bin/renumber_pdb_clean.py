#!/usr/bin/env python2

import swalign
from . import generic_functions

# global paramaters
SW_SCORE = swalign.NucleotideScoringMatrix(2, -1)
SMITH_WATERMAN = swalign.LocalAlignment(SW_SCORE)


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def grab_anarci_numbers(anarci_file):
    start = False

    alpha_num = []
    beta_num = []

    alpha_seq = ""
    beta_seq = ""
    previous_number = ""

    with open(anarci_file) as f:
        for line in f:
            if line.startswith("D") or line.startswith("A"):
                split_line = line.split()

                number = split_line[1]
                aa     = split_line[2]

                if aa is not "-":
                    if len(split_line) == 4:
                        number += split_line[2]
                    alpha_num.append(number)
                    alpha_seq += aa
                
                previous_number = number
 
            if line.startswith("B"):
                split_line = line.split()

                number = split_line[1]
                aa     = split_line[2]

                if aa is not "-":
                    if len(split_line) == 4:
                        number += split_line[2]
                    beta_num.append(number)
                    beta_seq += aa
                    
                previous_number = number
    return alpha_num, beta_num, alpha_seq, beta_seq


def filter_alignment(anarci, pdb):
    swout = SMITH_WATERMAN.align(anarci, pdb)
    offset = swout.q_pos - swout.r_pos
    
    attrs = vars(swout)
    print(', '.join("%s: %s" % item for item in list(attrs.items())))

    return offset


def get_pdb_list(pdb_file):
    pdb = []

    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                pdb.append(line)
    return pdb


def get_end_of_constant(pdb_nums, anarci_nums):
    #convert to int
    check_range = list(map(int, pdb_nums))

    #convert to int, for ANARCI files we need to remove AB conformations in CDR loops
    repls = ("A", ""), ("B", "")
    anarci_ints = [x.replace("A", "") for x in anarci_nums]
    anarci_ints = [x.replace("B", "") for x in anarci_ints]
    anarci_ints = list(map(int, anarci_ints))

    #get last variable ANARCI annotation
    final_variable = max(anarci_ints)

    #get number of AA past Anarci
    num_aa = len(list([x for x in check_range if x > final_variable]))


    #make a range from end of constant -> end of sequence
    start_constant = final_variable + 1
    end_constant = start_constant + num_aa

    constant_list = list(range(start_constant, end_constant))
    constant_list = list(map(str, constant_list))

    return constant_list



def renumber_chain(anarci_num, pdb_num, aa, current_chain, pdb_list):
    """
    make a dictionary that the key is the original residue number and amino acid, and the value is the new number,
    if a tuple of resnum and aa in dictionary, replace the number and add to a new list, return list twice in
    two functions and then make a second function getting mhca, mhcb and peptide
    """

    reference = {}
    out = []

    constant_list = get_end_of_constant(pdb_num, anarci_num)
    anarci_num = anarci_num + constant_list

    for x, y, z in zip(anarci_num, pdb_num, aa):
        reference[(y, z)] = x

    for line in pdb_list:
        chain = line[21]
        
        if chain == current_chain:
            amino_acid = generic_functions.three2one[line[17:20].strip()]
            res_num = line[22:26].strip()

            if (res_num, amino_acid) in reference:
                new_number = reference[(res_num, amino_acid)]

                line = list(line)

                for i in range(22, 26):
                    line[i] = " "
                
                start_idx = 23

                for i in range(len(new_number)):
                    line[start_idx] = new_number[i]
                    start_idx += 1
                
                out.append("".join(line))
    return out
            

def get_remaining_chains(pdb_list, mhca_chain, mhcb_chain, pep_chain):
    out = []

    chains = [mhca_chain, mhcb_chain, pep_chain]

    for line in pdb_list:
        chain = line[21]
        if chain in chains:
            out.append(line)
    
    return out

def renumber(anarci_file, pdb_file, alpha_chain, beta_chain, mhca_chain, mhcb_chain, pep_chain, outfile):
    alpha_anarci_num, beta_anarci_num, alpha_anarci_seq, beta_anarci_seq = grab_anarci_numbers(anarci_file)
    alpha_pdb_num, beta_pdb_num, alpha_pdb_seq, beta_pdb_seq = generic_functions.grab_pdb_numbers(pdb_file, alpha_chain, beta_chain)

    alpha_offset = filter_alignment(alpha_anarci_seq, alpha_pdb_seq)
    beta_offset = filter_alignment(beta_anarci_seq, beta_pdb_seq)

    alpha_pdb_num = alpha_pdb_num[alpha_offset:]
    alpha_pdb_seq = alpha_pdb_seq[alpha_offset:]

    beta_pdb_num = beta_pdb_num[beta_offset:]
    beta_pdb_seq = beta_pdb_seq[beta_offset:]

    pdb_list = get_pdb_list(pdb_file)

    alpha = renumber_chain(alpha_anarci_num, alpha_pdb_num, alpha_pdb_seq, alpha_chain, pdb_list)
    beta = renumber_chain(beta_anarci_num, beta_pdb_num, beta_pdb_seq, beta_chain, pdb_list)

    remaining = get_remaining_chains(pdb_list, mhca_chain, mhcb_chain, pep_chain)

    out_item = alpha + beta + remaining

    with open(outfile, "w") as f:
        for i in out_item:
            f.write(i)

