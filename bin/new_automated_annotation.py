import subprocess
import os
import swalign
import sys
import shutil
import pymol

import numpy as np

from Bio.SeqIO import convert
from itertools import product
from .data.mhc_sequences import mhc_fastas

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

def _call_ANARCI(fasta):
    subprocess.call(["ANARCI", "-r", "tr", "--sequence", fasta, "--outfile", "ANARCI.txt"])

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


def _pair_tcrs(alphas, betas, pdb):
    if len(alphas) is 1 and len(betas) is 1:
        print("There is only one TCR alpha and beta chain pairing possible")
        return tuple([alphas[0], betas[0]])

    if len(alphas) is 0 or len(betas) is 0:
        os.remove("ANARCI.txt")
        os.remove("clean.fasta")
        os.remove(pdb)

        sys.exit("No pair of TCRs found")

    #if not we have to start pairing them.
    pairs = [list(product([ai], betas)) for ai in alphas]
    contacts = []

    tcrs = [item for sublist in pairs for item in sublist]

    for pair in tcrs:
        command = "ncont XYZIN %s <<eof > ab_contact.txt \n" \
                      "source /*/%s/85-110/\n" \
                      "target /*/%s/85-110/\n" \
                      "mindist 0.0\n" \
                      "maxdist 22.0" %(pdb, pair[0], pair[1])
        os.system(command)

        total = 0
        with open("ab_contact.txt") as f:
            for line in f:
                if "CYS" in line:
                    total += 1

        contacts.append(total)

    pair = tcrs[contacts.index(max(contacts))]

    return pair


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


def init_pymol():
    print("\nInitialising pymol...\n")
    pymol.finish_launching(['pymol', '-nqc'])
    pymol.cmd.reinitialize()

    return np.array([float(x), float(y), float(z)])


def dist(residue_one, residue_two):
    diff = residue_one - residue_two
    return np.sqrt(np.sum(diff * diff))


def COM(pdb, pairs_of_chains):
    pymol.cmd.delete("all")
    print("Loading PDB")
    pymol.cmd.load(pdb)

    distances = []
    for pair in pairs_of_chains:

        pymol.cmd.select("to", "chain %s" % pair[0])
        pymol.cmd.select("from", "chain %s" % pair[1])

        to_COM = np.asarray(pymol.cmd.centerofmass("to"), dtype=np.float32)
        from_COM = np.asarray(pymol.cmd.centerofmass("from"), dtype=np.float32)

        distance = dist(to_COM, from_COM)

        distances.append(distance)

    return distances



def _pair_mhcs(alphas, betas, pdb, mhc_class):
    if len(alphas) is 1 and len(betas) is 1:
        print("There is only one MHC alpha and beta chain pairing possible")
        return [tuple([alphas[0], betas[0]])]

    potential_pairs = []

    #if not we have to start pairing them.
    pairs = [list(product([ai], betas)) for ai in alphas]
    mhcs = [item for sublist in pairs for item in sublist]

    distances = COM(pdb, mhcs)


    for a in alphas:
        subpair = []
        subcont = []

        for p, c in zip(mhcs, distances):
            if a == p[0]:
                subpair.append(p)
                subcont.append(c)

        idx = subcont.index(min(subcont))
        potential_pairs.append(subpair[idx])

    return potential_pairs


def _map_tcr_to_mhc(mhc_pairs, alpha, beta, pdb, mhc_class):
    contacts = []

    """
    lazy to not account for single pairing but it shouldn't fail
    if len(mhc_pairs) == 1:
       return mhc_pairs[0][0], mhc_pairs[0][1]
    """
    #todo skip contacts for 5 chain complex

    if mhc_class == 1:
        a1 = "50-86"
        a2 = "140-176"
    elif mhc_class == 2:
        a1 = "46-78"
        a2 = "54-91"
    else:
        sys.exit("MHC class unknown")

    for pairing in mhc_pairs:

        if mhc_class == 1:
            command = "ncont XYZIN %s <<eof > ab_contact.txt \n" \
                      "source /*/%s/40-60\n" \
                      "source /*/%s/40-60\n" \
                    "target /*/%s/%s/\n" \
                    "target /*/%s/%s/\n" \
                    "mindist 0.0\n" \
                    "maxdist 4.0" %(pdb, alpha, beta, pairing[0], a1, pairing[0], a2)
        else:
            command = "ncont XYZIN %s <<eof > ab_contact.txt \n" \
                      "source /*/%s/40-60\n" \
                      "source /*/%s/40-60\n" \
                      "target /*/%s/%s/\n" \
                      "target /*/%s/%s/\n" \
                      "mindist 0.0\n" \
                      "maxdist 4.0" % (pdb, alpha, beta, pairing[0], a1, pairing[1], a2)

        os.system(command)
        with open("ab_contact.txt") as f:
            for line in f:
                if "Total" in line:
                    contacts.append(int(line.split()[1]))

    return mhc_pairs, contacts


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

    pep = peptides[contacts.index(max(contacts))]
    return pep


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

                            for i, j in zip(list(range(0, 4)), list(range(23, 27))):
                                line[j] = insert_index[i]
                            
                            line = "".join(line)
                            out.write(line)

                            previous_string = current_string


def annotate_complex(pdb_file, filtered_name, numbered_name):
    chains = _get_chains(pdb_file)
    print("PDB contains %i chains" %len(chains))

    outdir = pdb_file.split("/")[-1].replace(".pdb", "")

    if len(chains) < 5:
        shutil.rmtree(outdir)
        sys.exit("Not enough chains")


    clean_pdb_standard(pdb_file, numbered_name)

    convert(numbered_name, "pdb-atom", "clean.fasta", "fasta")
    alphas, betas = _call_ANARCI("clean.fasta")

    removable = alphas + betas

    not_tcrs = _non_tcr_chains(chains, alphas, betas)

    tcr_pair = _pair_tcrs(alphas, betas, numbered_name)

    print("The following is the TCR pair being used: %s and %s" %(tcr_pair[0], tcr_pair[1]))

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

    print("The following chains are peptides: %s" % (",".join(peps)))
    print("The following chains are MHCas: %s" % (",".join(mhcas)))
    print("The following chains are MHCbs/B2Ms: %s" % (",".join(mhcbs)))


    mhc_pairs = _pair_mhcs(mhcas, mhcbs, numbered_name, mhc_class)


    mhc_pairings, contact_no =_map_tcr_to_mhc(mhc_pairs, tcr_pair[0], tcr_pair[1], numbered_name, mhc_class)

    index = contact_no.index(max(contact_no))
    chosen_mhc = mhc_pairings[index]

    peptide = _map_peptide(tcr_pair[0], tcr_pair[1], peps, numbered_name)

    return tcr_pair[0], tcr_pair[1],  peptide, chosen_mhc[0], chosen_mhc[1], mhc_class
