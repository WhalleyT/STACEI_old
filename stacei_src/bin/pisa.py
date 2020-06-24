#!/usr/bin/env python

from Bio.PDB import *
from Bio.SeqUtils import *
import subprocess
import os
import re


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


def get_peptide_lines(pdb_name, pdb_file, peptide_res):
    # get the peptide sequences
    parser = PDBParser()
    struct = parser.get_structure(pdb_name, pdb_file)
    pep_chain = struct[0][peptide_res]
    pep_residues = list(Selection.unfold_entities(pep_chain, 'R'))
    return pep_residues


def get_peptide_residues(pep_residues):
    pep_sequence = []
    for i in range(0, len(pep_residues)):
        if is_aa(pep_residues[i]):
            line = str(pep_residues[i])
            line = line[8:]
            line = line.replace(' het=  resseq=', ' ')
            line = line.replace(' icode= >', ' ')
            line = ''.join([i for i in line if not i.isdigit()])
            line = line.strip()
            pep_sequence.append(line)
    return pep_sequence


def pep_list_to_string(pep_sequence):
    pepstring = ""
    for i in pep_sequence:
        i = seq1(i)
        pepstring += i

    peplen = len(pepstring)
    return pepstring, peplen


def get_ATOMIC_data_PDB(pdb_file, peptide_res, MHCa_res, MHCb_res):
    file = open(pdb_file, "r")
    ATOM_lines = []

    for line in file:
        if line.startswith('ATOM'):
            line_p = re.sub(r'\s+', ' ', line)
            letter = line_p.split()[4]
            if letter == peptide_res or letter == MHCa_res or letter == MHCb_res:
                ATOM_lines.append(line)
    file.close()
    return ATOM_lines


def get_ATOMIC_TCR(pdb_file, tcra_res, tcrb_res):
    file = open(pdb_file, "r")
    atom_lines = []

    for line in file:
        if line.startswith('ATOM'):
            line_p = re.sub(r'\s+', ' ', line)
            letter = line_p.split()[4]
            if letter == tcra_res or letter == tcrb_res:
                atom_lines.append(line)
    file.close()
    return atom_lines


def get_header_PDB(pdb_file):
    file2 = open(pdb_file, "r")
    header_lines = []

    for line in file2:
        if not line.startswith('ATOM'):
            header_lines.append(line)
        else:
            break
    file2.close()
    return header_lines


def get_tail_PDB(pdb_file):
    last_atom = 0
    tail_lines = []

    # get the stuff at the end of the file, PISA doesn't seem to like it missing
    with open(pdb_file) as myFile:
        for num, line in enumerate(myFile, 1):
            if 'ATOM' or 'ANISOU' in line:
                last_atom = num
    myFile.close()

    with open(pdb_file) as myFile:
        for num, line in enumerate(myFile, 1):
            if num > last_atom:
                tail_lines.append(line)

    return last_atom, tail_lines


def write_parsed_pdb_file(header_lines, tail_lines, atom_lines, pisa_path):
    f = open(pisa_path + '/parsed_pdb.pdb', "w")
    for line in header_lines:
        f.write(line)
    for line in atom_lines:
        f.write(line)
    for line in tail_lines:
        f.write(line)
    f.close()
    return


def run_PISA(pisa_path):
    subprocess.call('pisa structural_tool -analyse %s/parsed_pdb.pdb' % str(pisa_path), shell=True)

    # get the monomers names <- don't think this may be needed
    subprocess.call('pisa  structural_tool  -list monomers > %s/monomers.txt' % str(pisa_path), shell=True)

    # extract monomer 3 details
    subprocess.call('pisa structural_tool -detail monomers 3 > %s/peptide_BSA.txt' % str(pisa_path), shell=True)


def run_TCR_pisa(pisa_path, pdb_file):
    # Call PISA
    print("Calling PISA for TCR")
    subprocess.call('pisa structural_tool2 -analyse %s' % str(pdb_file), shell=True)

    # get the monomers names <- don't think this may be needed
    subprocess.call('pisa  structural_tool2  -list monomers > %s/monomers.txt' % str(pisa_path), shell=True)

    # extract monomer 3 details
    subprocess.call('pisa structural_tool2 -detail monomers 2 > %s/alpha_peptide_BSA.txt' % str(pisa_path), shell=True)
    subprocess.call('pisa structural_tool2 -detail monomers 1 > %s/beta_peptide_BSA.txt' % str(pisa_path), shell=True)


def clean_BSA(pisa_path):
    peptide_BSA = []

    with open(pisa_path + '/peptide_BSA.txt') as peptide_file:
        for num, line in enumerate(peptide_file):
            peptide_BSA.append(line)
            if '  #' in line:
                start = num
            if '       Total ' in line:
                BSA = line
                BSA = BSA.split()[3]
                # print "Buried surface area of peptide is %s" %str(BSA)
                with open(pisa_path + "/BSA.txt", "w") as f:
                    f.write(BSA)
                f.close()
                end = num

    peptide_BSA_filtered = []

    for i in range(start + 2, end - 1):
        regexed_line = peptide_BSA[i].replace('|', '')
        peptide_BSA_filtered.append(regexed_line)

    f = open(pisa_path + '/peptide_BSA.txt', 'w')
    for line in peptide_BSA_filtered:
        f.write(line)
    f.close()
    return


def clean_ab_BSA(pisa_path):
    files = [pisa_path + "/alpha_peptide_BSA.txt", pisa_path + "/beta_peptide_BSA.txt"]
    out = open(pisa_path + "/TCR_BSA.txt", "w")
    for file in files:
        lines = []
        with open(file, "r") as f:
            start = False
            for line in f:
                if "   #" in line:
                    start = True

                x = line.strip()
                if len(x) > 0:
                    if start and x[0].isdigit():
                        line = line.replace('|', '')
                        line = re.sub(r"\s+", " ", line)
                        line = re.sub(":", " ", line)
                        out.write(line + "\n")


# print "alpha and beta TCR PISA output cleaned"


def plot(pisa_path):
    # print "Plotting peptide BSA in R"
    subprocess.call('Rscript %s/bin/R/peptide_BSA.R %s' % (str(os.getcwd()), pisa_path), shell=True)


# print "Peptide BSA plot written to %s" %pisa_path

def plot_TCR(TCR_path, alpha, beta, fasta):
    # print "Plotting TCR BSA in R"
    subprocess.call('Rscript %s/bin/R/TCR_BSA.R %s %s %s %s' % (
    str(os.getcwd()), TCR_path, alpha, beta, str(os.getcwd() + "/" + fasta)), shell=True)


def _find_between(s, first, last):
    try:
        start = s.index(first) + len(first)
        end = s.index(last, start)
        return s[start:end]
    except ValueError:
        return ""


def _get_CDR_loop(chain, chain_str):
    chain = chain.split("|")
    ranges = []

    loops = ["CDR1", "CDR2", "CDR3", "FW"]
    loops = [s + chain_str for s in loops]

    cdrs = []
    for loop in loops:
        for element in chain:
            if loop in element:
                cdrs.append(element)

    for cdr in cdrs:
        # print cdr
        loop = _find_between(cdr, "[", "]")
        loop = loop.split(",")
        loop = loop[::len(loop) - 1]
        ranges.append(loop)

    return ranges


def _get_lines(ranges, chain, file):
    sums = []
    for rang in ranges:
        BSA_sum = 0
        with open(file) as f:
            for line in f:
                line = line.split()
                if line[1] == chain:
                    if line[3].isdigit():
                        if int(rang[0]) <= int(line[3]) <= int(rang[1]):
                            BSA_sum += float(line[5])
        sums.append(BSA_sum)
    return sums


def sum_BSA(alpha, beta, BSA, fasta):
    # get lines we care about
    alpha_line = ""
    beta_line = ""

    with open(fasta) as fasta_file:
        for name, sequence in read_fasta(fasta_file):
            if name.split("|")[1] == alpha:
                alpha_line = name
            elif name.split("|")[1] == beta:
                beta_line = name

    alpha_ranges = _get_CDR_loop(alpha_line, "a")
    beta_ranges = _get_CDR_loop(beta_line, "b")

    alpha_sums = _get_lines(alpha_ranges, alpha, BSA + "/TCR_BSA.txt")
    beta_sums = _get_lines(beta_ranges, beta, BSA + "/TCR_BSA.txt")

    out_b = ["CDR1a: ", "CDR2a: ", "CDR3a: ", "CDRFWa: "]
    out_a = ["CDR1b: ", "CDR2b: ", "CDR3b: ", "CDRFWb: "]

    out = open(BSA + "/CDR_loop_BSA.txt", "w")

    for string, summ in zip(out_b, alpha_sums):
        out.write(string + str(summ) + "\n")
    for string, summ in zip(out_a, beta_sums):
        out.write(string + str(summ) + "\n")

# print "BSA sums for CDR loops written to 'CDR_loop_BSA.txt'"

def plot_sum(pdb):
    subprocess.call('Rscript %s/bin/R/CDR_loop_BSA.R %s' % (str(os.getcwd()), str(os.getcwd()) + "/" + pdb), shell=True)
