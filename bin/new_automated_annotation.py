import operator
import subprocess
import math
import os
import swalign
import sys
import re
import shutil

from Bio.SeqIO import convert
from itertools import product
from data.mhc_sequences import mhc_fastas


def _three_dimensional_distance(a, b):
    #distance = root(x^2 + y^2 + z2) where x = x1 -x2 etc
    #https://math.stackexchange.com/questions/42640/calculate-distance-in-3d-space
    if len(a) != len(b):
        raise Exception("These two vectors appear to be of different lengths")

    if len(a) != 3 or len(b) != 3:
        raise Exception("One of your vectors are not of length 3, are you sure these are 3d coordinates")

    difference = map(operator.sub, a, b)
    diffsquare = [x*x for x in difference]
    return math.sqrt(sum(diffsquare))

def _read_fasta(fp):
        name, seq = None, []
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name:
            yield (name, ''.join(seq))


def _anarci_installed():
    try:
        subprocess.call(["ANARCI"])
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            print("ANARCI is not installed")
            sys.exit()
        else:
            print("Something else is wrong")
            raise


def _get_chains(pdb):
    chains = []

    with open(pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21]
                chains.append(chain)

    return list(set(chains))

def _clean_pdb_new(pdb, outfile):

    out = open(outfile, "w")

    solvents_to_remove = ["NAG", "NA", "EMC", "GOL", "SO4", "PGE", "7PE", "TAM",
                          "EPE", "EDO", "PG4", "BMA", "MAN", "MLI"]
    second_conformation = ["CD2", "CE2", "OE2", "NE2", "OG2", "OD2", "NH2",
                           "NE2", "ND2", "CG2 ", "CZ2", "CZ3", "CH2"]

    amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX',
                   'CYS', 'GLU', 'GLN', 'GLX', 'GLY',
                   'HIS', 'ILE', 'LEU', 'LYS', 'MET',
                   'PHE', 'PRO', 'SER', 'THR', 'TRP',
                   'TYR', 'VAL']

    with open(pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                aa = line[16:20]
                if aa.startswith(" ") is False:
                    if aa.startswith("A"):
                        aa = " " + aa[1:]
                        line = line[:16] + aa + line[20:]

                if aa.startswith("B") is False:
                    out.write(line)

def untangle(i):
    # if our coordinates are joined togther, seperate them
    if len(i) < 12:
        new = []
        for item in i:
            if item.count(".") > 1:
                splitter = item.find(".") + 3
                first = item[0:splitter]
                second = item[splitter:]
                new.append(first)
                new.append(second)
            else:
                new.append(item)
        return new
    else:
        return i


def _list_to_pdb_spec(pdb, out):

    outfile = open(out, "w")

    for i in pdb:
        i = i.split()
        if i[-1].isalpha():
            count = 12
        else:
            count = 11

        while len(i) < count:
            i = untangle(i)

        if len(i) == 12:
            outstr = i[0] + (7 - len(str(i[1]))) * " " + str(i[1]) + 2 * " " + str(i[2]) + \
                     (4 - len(str(i[2]))) * " " + str(i[3]) + " " + str(i[4]) + " " * (4- len(str(i[5]))) + \
                     str(i[5]) + (11 - len(str(i[6]))) * " " + str(i[6]) + (8 - len(str(i[7]))) * " " + str(i[7]) + (8 - len(str(i[8]))) * \
                     " " + str(i[8]) + (6 - (len(str(i[9])))) * " " + str(i[9]) + (6 - (len(str(i[10])))) * " " \
                     + str(i[10]) + (12 - len(str(i[11]))) * " " + str(i[11])
        else:
            outstr = i[0] + (7 - len(str(i[1]))) * " " + str(i[1]) + 2 * " " + str(i[2]) + \
                     (4 - len(str(i[2]))) * " " + str(i[3]) + " " + str(i[4]) + " " * (4- len(str(i[5]))) + \
                     str(i[5]) + (11 - len(str(i[6]))) * " " + str(i[6]) + (8 - len(str(i[7]))) * " " + str(i[7]) + (8 - len(str(i[8]))) * \
                     " " + str(i[8]) + (6 - (len(str(i[9])))) * " " + str(i[9]) + (6 - (len(str(i[10])))) * " " \
                     + str(i[10])
        outfile.write(outstr + "\n")


def _read_chain(pdb, chain):
    pdb_list = []
    with open(pdb) as infile:
        for line in infile:
            if line.startswith("ATOM") and line.split()[4] == chain:
                pdb_list.append(line)
    return pdb_list


def _get_2dim_list(pdb_list):

    first = True

    amino_acid = []
    chains = []

    previous_residue = ""
    for line in pdb_list:

        residue = line.split()[5]
        if first:
            amino_acid.append(line)
        else:
            if residue == previous_residue:
                amino_acid.append(line)
            else:
                chains.append(amino_acid)
                amino_acid = [line]

        first = False
        previous_residue = residue
    chains.append(amino_acid)
    return chains


def _renumber_pdb(pdb, chains, outname):
    atoms = []

    for annot in chains:
        pdb_list = _read_chain(pdb, annot)
        chain = _get_2dim_list(pdb_list)

        if chain[0][0].split()[3] == "MET":
            index = 0
        else:
            index = 1

        for AA in chain:
            for line in AA:
                split = line.split()
                old = " " + str(split[5]) + " "
                new = " " + str(index) + " "
                len_difference = len(new) - len(old)
                if len_difference == 1:
                    new = new[1:]

                line = line.replace(old, new)
                atoms.append(line)
            index += 1

    _list_to_pdb_spec(atoms, outname)



def _call_ANARCI(fasta, chains):
    subprocess.call(["ANARCI", "--sequence", fasta, "--outfile", "ANARCI.txt"])

    alphas = []
    betas = []

    tcr = None
    current_AA = None

    with open("ANARCI.txt") as f:
        for line in f:
            if line.startswith("#"):
                if ":" in line:
                    current_AA = line.split(":")[1]
                    tcr = None
                if "|" in line:
                    split = line.split("|")[2]
                    if split in ["A", "B", "D"]:
                        tcr = split

            if tcr is not None and current_AA is not None:
                if tcr is "A" or tcr is "D": #some delta TCRs get miss-annotated due to TRAV29-DV5 similarity
                    alphas.append(current_AA.strip())
                if tcr is "B":
                    betas.append(current_AA.strip())

    return list(set(alphas)), list(set(betas))


def _non_tcr_chains(chains, alphas, betas):
    return list(set(chains) - set(alphas + betas))


def is_pair_joined(pair):
    print "dunno"

def _pair_tcrs(alphas, betas, pdb):
    if len(alphas) is 1 and len(betas) is 1:
        print "There is only one TCR alpha and beta chain pairing possible"
        return alphas[0], betas[0]

    #if not we have to start pairing them.
    pairs = [list(product([ai], betas)) for ai in alphas]
    tcrs  = []
    for pair in pairs:
        tcrs.append(pair[0])
        tcrs.append(pair[1])

    contacts = []

    for pair in tcrs:
        command = "ncont XYZIN %s <<eof > ab_contact.txt \n" \
                      "source /*/%s\n" \
                      "target /*/%s\n" \
                      "mindist 0.0\n" \
                      "maxdist 4.0" %(pdb, pair[0], pair[1])
        os.system(command)

        with open("ab_contact.txt") as f:
            for line in f:
                if "Total" in line:
                    contacts.append(line.split()[1])

    pair = tcrs[contacts.index(max(contacts))]

    return [x for y, x in sorted(zip(contacts, tcrs), reverse=True)]


def  _get_mhc_pep(none_tcr, fasta):
    MATCH = 2
    MISMATCH = -1
    SW_SCORE = swalign.NucleotideScoringMatrix(MATCH, MISMATCH)
    SMITH_WATERMAN = swalign.LocalAlignment(SW_SCORE)
    b2m_found = False

    proteins = {}

    peptides = []
    mhcas = []
    mhcbs = []

    mhc_class = None

    with open(fasta) as f:
        for name, seq in _read_fasta(f):
            letter = name.split(":")[1]
            if letter in none_tcr:
                proteins[letter] = seq

    for test in proteins:
        scores = []
        refs = []
        for reference in mhc_fastas:
            refs.append(reference)
            ref_protein = mhc_fastas[reference]
            hit_protein = proteins[test]

            swout = SMITH_WATERMAN.align(ref_protein, hit_protein)
            swout = swout.score
            scores.append(swout)

        max_score = max(scores)
        max_index = scores.index(max(scores))

        #print "Maximum score for protein %s is %i, which corresponds to %s" %(test, max_score, refs[max_index])

        if max_score < 50:
            peptides.append(test)
        else:
            ref_name = refs[max_index]

            if ref_name == "B2M":
                mhcbs.append(test)
                b2m_found = True
            elif "D" in ref_name.split("*")[0] and "B" in ref_name.split("*")[0]: #contains a d for class II and b for beta chain
                mhcbs.append(test)
            else:
                mhcas.append(test)

    if b2m_found:
        mhc_class = 1
    else:
        mhc_class = 2

    return peptides, mhcas, mhcbs, mhc_class


def _pair_mhcs(alphas, betas, pdb, alpha_range, beta_range, distance):
    if len(alphas) is 1 and len(betas) is 1:
        print "There is only one MHC alpha and beta chain pairing possible"
        return [tuple([alphas[0], betas[0]])]

    potential_pairs = []

    #if not we have to start pairing them.
    pairs = [list(product([ai], betas)) for ai in alphas]
    mhcs  = []
    for pair in pairs:
        mhcs.append(pair[0])
        mhcs.append(pair[1])

    contacts = []

    for pair in mhcs:
        command = "ncont XYZIN %s <<eof > ab_contact.txt \n" \
                      "source /*/%s/%s\n" \
                      "target /*/%s/%s\n" \
                      "mindist 0.0\n" \
                      "maxdist %s" %(pdb, pair[0], alpha_range, pair[1], beta_range, distance)
        os.system(command)

        with open("ab_contact.txt") as f:
            for line in f:
                if "Total" in line:
                    contacts.append(line.split()[1])

    for a in alphas:
        subpair = []
        subcont = []

        for p, c in zip(mhcs, contacts):
            if a == p[0]:
                subpair.append(p)
                subcont.append(c)

        idx = subcont.index(max(subcont))
        potential_pairs.append(subpair[idx])

    return potential_pairs


def _map_tcr_to_mhc(mhc_pairs, alpha, beta, pdb):
    contacts = []
    if len(mhc_pairs) == 1:
       return mhc_pairs[0][0], mhc_pairs[0][1], True

    for pairing in mhc_pairs:
        command = "ncont XYZIN %s <<eof > ab_contact.txt \n" \
                  "source /*/%s/40-60\n" \
                  "source /*/%s/40-60\n" \
                  "target /*/%s\n" \
                  "target /*/%s\n" \
                  "mindist 0.0\n" \
                  "maxdist 8.0" %(pdb, alpha, beta, pairing[0], pairing[1])
        os.system(command)
        with open("ab_contact.txt") as f:
            for line in f:
                if "Total" in line:
                    contacts.append(int(line.split()[1]))

    pair = mhc_pairs[contacts.index(max(contacts))]
    if all(v == 0 for v in contacts):
        success = False
    else:
        success = True

    
    return pair[0], pair[1], success


def _map_peptide(tcra, tcrb, peptides, pdb):
    if len(peptides) == 1:
       return peptides[0]

    contacts = []

    for peptide in peptides:
        command = "ncont XYZIN %s <<eof > ab_contact.txt \n" \
                  "source /*/%s/95-120\n" \
                  "source /*/%s/95-120\n" \
                  "target /*/%s\n" \
                  "mindist 0.0\n" \
                  "maxdist 4.0" %(pdb, tcra, tcrb, peptide)

        os.system(command)

        with open("ab_contact.txt") as f:
            for line in f:
                if "Total" in line:
                    contacts.append(line.split()[1])
                    print peptide, tcra, tcrb, line
    pep = peptides[contacts.index(max(contacts))]
    return pep


def remove_conformations(pdb_list):
    """
    Remove uncertain conformations (messes up some tools)
    and replaces them just with the A set
    """
    lines = []

    for line in pdb_list:
        if line.startswith("ATOM"):
            if line[16] == " ":
                lines.append(line)
            if line[16] == "A":
                # python doesn't allow direct conversion; convert to list
                line = list(line)
                line[16] = " "
                line = "".join(line)
                lines.append(line)
    return lines


def get_default_residues(pdb_list, chain):
    """
    Gets length of PDB chain
    """
    started = False
    numbers = set()

    for line in pdb_list:
        if chain == line[21]:
            if started is False:
                start = line[22:27].strip()
                started = True

            if line[22:27].isalpha() == False:
                number = line[22:27].strip()
                numbers.add(number)
    return len(numbers), start


def renumber(pdb_list, chain, residue_nums, start):
    """
    Takes a pdb list file, a chain and a
    corresponding list of residues to replace
    in the PDB file
    """

    previous = start

    res_iter = iter(residue_nums)
    new_residue = next(res_iter)

    out = []

    for line in pdb_list:
        if chain == line[21]:
            current = line[22:27].strip()
            if current == previous:
                line = list(line)
                size = len(str(new_residue))
                spaces = 4 - size
                spaces = " " * spaces
                string = spaces + str(new_residue) + "  "

                start_char = 22
                for i in string:
                    line[start_char] = i
                    start_char += 1
                line = "".join(line)
            else:
                new_residue = next(res_iter)
                line = list(line)
                size = len(str(new_residue))
                spaces = 4 - size
                spaces = " " * spaces
                string = spaces + str(new_residue) + "  "

                start_char = 22
                for i in string:
                    line[start_char] = i
                    start_char += 1
                line = "".join(line)

            previous = current
            out.append(line)
    return out


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


def get_fasta_annot(fasta_file, chain):
    indexes = []
    with open(fasta_file) as f:
        for name, seq in read_fasta(f):
            if name.split("|")[:-2][1] == chain:
                #chains found
                if name.split("|")[2] in ["TCRA", "TCRB"]:
                    for i,j in enumerate(seq):
                        if j != ".":
                            indexes.append(i)
                else:
                    indexes = range(1, len(seq) + 1)
    return indexes


def clean_pdb(pdb_file, chains, fasta, linear, outname):
    pdb = []
    with open(pdb_file) as f:
        for line in f:
            pdb.append(line)

    pdb = remove_conformations(pdb)
    pdb_2D = []

    if linear:
        for chain in chains:
            size, start = get_default_residues(pdb, chain)
            resi = list(range(1, size + 1))
            subpdb = renumber(pdb, chain, resi, start)
            pdb_2D.append(subpdb)

            outfile = open(outname, "w")

            for i in pdb_2D:
                for j in i:
                    outfile.write(j)
    else:
        for chain in chains:
            size, start = get_default_residues(pdb, chain)
            resi = get_fasta_annot(fasta, chain)
            subpdb = renumber(pdb, chain, resi, start)
            pdb_2D.append(subpdb)

            outfile = open(outname, "w")

            for i in pdb_2D:
                for j in i:
                    outfile.write(j)


def clean_pdb_standard(infile, outfile):
    
    chains = set()
    #grab all possible chains
    with open(infile) as f:
        for line in f:
            if line.startswith("ATOM"):
                chains.add(line[21])

    out = open(outfile, "w")

    for chain in chains:
        index = 0
        previous_string = None

        with open(infile) as f:
            for line in f:
                if line.startswith("ATOM"):
                    if chain == line[21]:
                        amino_acid = list(line[16:20])
                        line       = list(line)
                        if amino_acid[0] == " " or amino_acid[0] == "A":
                            line[16] = " "
                            #now we have fixed the amino acid, now onto the number
                            current_string = "".join(line[23:27])

                            if current_string != previous_string:
                                index += 1

                            insert_index = list(str(index))
                            spaces = 3 - len(insert_index)
                            spaces = " " * spaces

                            insert_index = list(spaces) + insert_index + list(" ")

                            for i, j in zip(range(0, 4), range(23, 27)):
                                line[j] = insert_index[i]
                            
                            line = "".join(line)
                            out.write(line)

                            previous_string = current_string

def annotate_complex(pdb_file, filtered_name, numbered_name):
    chains = _get_chains(pdb_file)
    print "PDB contains %i chains" %len(chains)

    if len(chains) < 5:
        outdir = pdb_file.split("/")[-1].replace(".pdb", "")
        shutil.rmtree(outdir)
        sys.exit("Not enough chains")


    clean_pdb_standard(pdb_file, numbered_name)

    convert(numbered_name, "pdb-atom", "clean.fasta", "fasta")
    alphas, betas = _call_ANARCI("clean.fasta", chains)

    removable = alphas + betas

    not_tcrs = _non_tcr_chains(chains, alphas, betas)
    tcr_pairs = _pair_tcrs(alphas, betas, numbered_name)


    if len(tcr_pairs) > 2:
        alphas, betas = zip(*tcr_pairs)
    else:
        alphas = list(tcr_pairs[0])
        betas = list(tcr_pairs[1])
        tcr_pairs = [tcr_pairs]

    none_tcr = list(set(chains) - set(removable))
    peps, mhcas, mhcbs, mhc_class =_get_mhc_pep(none_tcr, "clean.fasta")

    #now let's pair our mhca and bs then see if they map to the tcr

    if not mhcas or not mhcbs:
        os.remove(numbered_name)
        os.remove("ANARCI.txt")
        os.remove("ab_contact.txt")
        os.remove("clean.fasta")

        outdir = pdb_file.split("/")[-1].replace(".pdb", "")
        shutil.rmtree(outdir)
        sys.exit("One of MHCa or MHCb could not be identified")
        
    if mhc_class == 1:
        print "Pairing MHC class I"
        alpha_range = "14-16"
        beta_range = "22-24"
        distance = "32.0"
        mhc_pairs = _pair_mhcs(mhcas, mhcbs, numbered_name, alpha_range, beta_range, distance)
    else:
        alpha_range = "28-30"
        beta_range = "38-40"
        distance = "22.0"
        mhc_pairs = _pair_mhcs(mhcas, mhcbs, numbered_name, alpha_range, beta_range, distance)
        print "Pairing MHC class II"

    for tcr in tcr_pairs:
        mhc_alpha, mhc_beta, success =_map_tcr_to_mhc(mhc_pairs, tcr[0], tcr[1], numbered_name)

        if success:
            print "Successfully paired a TCR to MHC"
            peptide = _map_peptide(tcr[0], tcr[1], peps, numbered_name)
            tcr_alpha = tcr[0]
            tcr_beta = tcr[1]
            return tcr_alpha, tcr_beta, peptide, mhc_alpha, mhc_beta, mhc_class
        else:
            print "Unsuccessfuly tried to pair a TCR and MHC, trying again"
            continue
    
    sys.exit("Could not successfully group a TCR-pMHC complex")
