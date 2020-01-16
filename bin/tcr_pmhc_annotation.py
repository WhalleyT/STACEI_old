import subprocess
import re
import os
import sys

from Bio.SeqIO import convert as bio_convert
from Bio.PDB import PDBParser

"""
##################################################
                  cleaning of PDB
##################################################
"""


def _split_joined_conformation_residue(line, amino_acids):
    for i in range(len(amino_acids)):
        if amino_acids[i] in line and " " + amino_acids[i] not in line:
            line = re.sub(amino_acids[i], " " + amino_acids[i], line)
    return line


def _parse_first_conformation(line, first_conformation):
    for i in range(len(first_conformation)):
        if first_conformation[i] in line:
            index = line.find(first_conformation[i])
            end = index + len(first_conformation[i]) - 1
            substitute = line[index:end]
            line = re.sub(substitute, substitute[:-1] + " ", line)
    return line


def clean_pdb(pdb_file, out_name):
    solvents_to_remove = ["NAG", "NA", "EMC", "GOL", "SO4", "PGE", "7PE", "TAM",
                          "EPE", "EDO", "PG4", "BMA", "MAN", "MLI"]
    second_conformation = ["CD2", "CE2", "OE2", "NE2", "OG2", "OD2", "NH2",
                           "NE2", "ND2", "CG2 ", "CZ2", "CZ3", "CH2"]

    amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX',
                   'CYS', 'GLU', 'GLN', 'GLX', 'GLY',
                   'HIS', 'ILE', 'LEU', 'LYS', 'MET',
                   'PHE', 'PRO', 'SER', 'THR', 'TRP',
                   'TYR', 'VAL']

    a_conformations = []
    b_conformations = []

    res = False
    for i in amino_acids:
        a_conf = "A"+i
        b_conf = "B"+i
        a_conformations.append(a_conf)
        b_conformations.append(b_conf)

    infile = open(pdb_file, "r")
    pdb = []

    for line in infile:
        pdb.append(line)
    infile.close()

    with open(out_name, 'w') as output:
        for line in pdb:
                if line.startswith("ATOM"):
                    line = _split_joined_conformation_residue(line, amino_acids)
                    res = line.split()[5].isdigit() is False
                    if res:
                        x = line


                if any(solvent in line for solvent in solvents_to_remove) is False and \
                   any(confor in line for confor in second_conformation) is False and res is False:
                    output.write(line)
                res = False
    return



def clean_pdb_new(pdb, outfile):

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




def remove_gap(pdb):
    """removes gap between MOLECULE ID and CHAIN ID"""
    out = open(pdb + "_test.pdb", "w")
    with open(pdb + ".pdb", "r") as pdb_file:
        hit_molecule = False
        previous_line = ""
        for line in pdb_file:
            if "MOLECULE: " in previous_line:
                hit_molecule = True
            if "CHAIN: " in line and hit_molecule:
                hit_molecule = False


            if hit_molecule is False:
                out.write(line)

            previous_line = line

def remove_gap_2(pdb):
    """removes gap between MOLECULE ID and CHAIN ID"""
    out = open(pdb + "_test.pdb", "w")
    with open(pdb + ".pdb", "r") as pdb_file:
        hit_molecule = False
        previous_line = ""
        for line in pdb_file:
            if "MOLECULE: " in previous_line:
                hit_molecule = True
                concat = previous_line
            if "CHAIN: " in line and hit_molecule:
                hit_molecule = False
                concat += previous_line
                concat = " ".join(concat.split())
                out.write(concat + "\n")
                concat = ""

            if hit_molecule is False:
                out.write(line)

            previous_line = line
"""
##################################################
                 Annotation of PDB
##################################################
"""


def can_be_fully_annotated(pdb):
    with open(pdb, "r") as structure:
        for line in structure:
            if "MOLECULE" in line:
                return True
    return False


def read_fasta(fp):
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


def anarci_installed():
    try:
        subprocess.call(["ANARCI"])
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            print("ANARCI is not installed")
            sys.exit()
        else:
            print("Something else is wrong")
            raise


def find_tcr(pdb):
    anarci_file = pdb + "_ANARCI.txt"

    alpha = ""
    beta = ""

    with open(anarci_file) as ANARCI:
        for line in ANARCI:
            if pdb in line:
                letter = line.split(pdb + ":")[-1].strip()
            if "|A|" in line:
                alpha = letter
            if "|B|" in line:
                beta = letter
            if alpha != "" and beta != "":
                return alpha, beta

    if alpha == "" or beta == "":
        print "One or more TCR chain could not be identified"
        sys.exit()


def find_extra_tcr(pdb, alpha, beta):
    anarci_file = pdb + "_ANARCI.txt"

    ignore = []

    with open(anarci_file) as ANARCI:
        for line in ANARCI:
            if pdb in line:
                letter = line.split(pdb + ":")[-1].strip()
            if "|A|" in line:
                if letter != alpha:
                    ignore.append(letter)
            if "|B|" in line:
                if letter != beta:
                    ignore.append(letter)
    return ignore


def find_peptide(pdb):
    pdb_fa = pdb + ".fasta"
    
    names = []
    seqs = []
    
    with open(pdb_fa, "r") as fasta:
        for name, sequence in read_fasta(fasta):
            names.append(name)
            seqs.append(len(sequence))

    protein = names[seqs.index(min(seqs))]
    return protein[-1]


def get_remaining_chains(pdb_name, tcr_and_peptide):
    # get letters in fastas
    fasta_residues = []

    with open(pdb_name + ".fasta", "r") as fasta:
        for line in fasta:
            if line.startswith(">") and "<unknown description" not in line:
                x = line.split(":")
                x = x[1]
                fasta_residues.append(x[0])
            if line.startswith(">") and "<unknown description" in line:
                fasta_residues.append(line[1])

    # run ANARCI again; removing tcrs we have chosen
    surplus_tcrs = find_extra_tcr(pdb_name, tcr_and_peptide[0], tcr_and_peptide[1])

    return list(set(fasta_residues) - set(tcr_and_peptide + surplus_tcrs))


def find_mhc(pdb, mhc_list):

    class1 = ["A", "B", "C", "E", "F", "G"]
    class1 = ["HLA-" + s for s in class1]
    class1 = class1 + ["CLASS ONE", "CLASS 1 ", "CLASS I ", "CLASS I-", "CLASS 1-"]
    b2m = ["B2M", "B-2-M", "MICROGLOBULIN"]

    class2 = ["DM", "DO", "DP", "DQ", "DR", "CLASS 2", "CLASS TWO",
              "CLASS II", "CLASS II-", "CLASS 2-", "A CHAIN ", "B CHAIN"]

    c2a = ["ALPHA"]
    c2b = ["BETA"]

    is_class1 = False
    is_b2m = False

    mhca = ""
    mhcb = ""

    potential_names = [">" + pdb + ":" + s for s in mhc_list]
    lens = []

    # let's check the length of our sequence is long enough to be MHC
    with open(pdb + ".fasta", "r") as fasta:
        for name, sequence in read_fasta(fasta):
            if name in potential_names:
                lens.append((name[-1], len(sequence)))

    mhc_list = []
    for double in lens:
        if double[1] > 50:
            mhc_list.append(double[0])

    a_found = False
    b_found = False

    with open(pdb + "_test.pdb", "r") as pdb:
        previous_line = ""
        for line in pdb:
            if "CHAIN: " in line:
                letter = line.split("CHAIN: ")[1].strip()[0]
                if letter in mhc_list:
                    if any(ID in previous_line for ID in class1) and not a_found:
                        is_class1 = True
                        mhca = letter
                        a_found = True
                    if any(ID in previous_line for ID in b2m) and not b_found:
                        is_b2m = True
                        mhcb = letter
                        b_found = True
                    if any(ID in previous_line for ID in class2):
                        if any(ID in previous_line for ID in c2a) and not a_found:
                            mhca = letter
                            a_found = True
                        if any(ID in previous_line for ID in c2b) and not b_found:
                            mhcb = letter
                            b_found = True

            previous_line = line

    if is_class1 and is_b2m:
        print "MHC is class I and MHC alpha is chain %s and beta-2-microglobulin is chain %s" %(mhca, mhcb)
        mhc_class = 1
        return mhca, mhcb, mhc_class
    if not is_class1 and not is_b2m:
        print "MHC is class II and MHC alpha is chain %s and MHC beta is chain %s" %(mhca, mhcb)
        mhc_class = 2
        return mhca, mhcb, mhc_class


def convert_pdb_to_fasta(pdb):
    pdb_name = pdb + "_filtered.pdb"
    pdb_fa = pdb + ".fasta"

    fail = False
    fasta = None
    
    try:
        fasta = bio_convert(pdb_name, "pdb-atom", pdb_fa, "fasta")
    except:
        # Change this exception when we find a case where the error happens
        print "Bio.convert was unsuccessful, trying by SeqRes"
        fail = True

    if fail:
        try:
            fasta = bio_convert(pdb_name, "pdb-seqres", pdb_fa, "fasta")
        except:
            print "Converting by SeqRes was also unsuccessful"
            sys.exit()
        return fasta
    else:
        return fasta


def clean_fasta(pdb):
    fasta = pdb + ".fasta"
    n = []
    s = []

    with open(fasta, "r") as f:
        for name, sequence in read_fasta(f):
            n.append(name)
            s.append(sequence)

    fasta_out = open(fasta, "w")

    for name, sequence in zip(n, s):
        if "<unknown description>" in name:
            name = ">" + pdb + ":" + name[1]
        if len(name) > len(pdb) + 3:
            before = name.split(":")[0]
            after = name.split(":")[1]
            name = before + ":" + after[0]
        fasta_out.write(name + "\n" + sequence + "\n")


def run_anarci(pdb):
    pdb_anarci = pdb + "_ANARCI.txt"    
    pdb_fa = pdb + ".fasta"
    subprocess.call(["ANARCI", "--sequence", pdb_fa, "--outfile", pdb_anarci])


def clean_up(pdb):
    os.remove(pdb + "_ANARCI.txt")
    os.remove(pdb + ".fasta")
    os.remove(pdb + "_filtered.pdb")
    os.remove(pdb + "_test.pdb")


def count_chains(pdb):
    chain_count = 0

    parser = PDBParser()
    structure = parser.get_structure('foobar', pdb)
    for model in structure:
        for chain in model:
            chain_count += 1
    
    if chain_count >= 5:
        print "There are %i chains in the structure" % chain_count
    else:
        print "There are < 5 chains in the structure %s, meaning it cannot be a TCR-pMHC" % pdb
        sys.exit()


def count_entities(pdb):
    entity_counter = 0

    with open(pdb, "r") as pdb_file:
        for line in pdb_file:
            if "CHAIN:" in line:
                entity_counter += 1

    if entity_counter < 5:
        print "There are %i entities." % entity_counter
        sys.exit()


def write(pdb, tcra, tcrb, peptide, mhca, mhcb, mhc_class):
    outfile = open(pdb, "w")

    outfile.write("TCR alpha chain: %s\n" % tcra)
    outfile.write("TCR beta chain: %s\n" % tcrb)
    outfile.write("peptide chain: %s\n" % peptide)

    if mhc_class == 1:
        outfile.write("MHC alpha chain: %s\n" % mhca)
        outfile.write("beta-2-microglobulin chain: %s\n" % mhcb)
        outfile.write("\n")
        outfile.write("MHC complex is class I\n")
        outfile.write("\n")
    elif mhc_class == 2:
        outfile.write("MHC alpha chain: %s\n" % mhca)
        outfile.write("MHC beta chain: %s\n" % mhcb)
        outfile.write("\n")
        outfile.write("MHC complex is class II\n")
        outfile.write("\n")


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


def _get_before_atoms(pdb):
    before = []
    with open(pdb) as infile:
        for line in infile:
            if line.startswith("ATOM"):
                return before
            before.append(line)


def _get_after_atoms(pdb):
    lines = []
    atom_start = False

    with open(pdb) as infile:
        for line in infile:
            if line.startswith("ATOM"):
                atom_start = True
            if atom_start:
                lines.append(line)
    # remove lines with atom in
    out = []
    for line in lines:
        if line.startswith("ATOM") is False:
            out.append(line)

    return out


def list_to_pdb_spec(pdb, out):

    outfile = open(out, "w")
    """
    ATOM    788  CD2 LEU A  99      88.036  -2.528  38.631  1.00 45.72           C
    ATOM    789  N   ARG A 100      89.799  -2.593  34.112  1.00 47.59           N
    ATOM    790  CA  ARG A 100      90.624  -2.462  32.905  1.00 50.00           C
    ATOM    791  C   ARG A 100      89.922  -2.918  31.614  1.00 47.10           C
    ATOM    792  O   ARG A 100      90.487  -2.806  30.524  1.00 46.66           O
    ATOM    793  CB  ARG A 100      91.923  -3.275  33.042  1.00 55.91           C
    """


    for i in pdb:
        i = i.split()

        #if our coordinates are joined togther, seperate them
        if len(i) < 12:
            new = []
            for item in i:
                if item.count(".") > 1:
                    splitter = item.find(".") + 3
                    first  = item[0:splitter]
                    second = item[splitter:]
                    new.append(first)
                    new.append(second)
                else:
                    new.append(item)
            i = new

        outstr = i[0] + (7 - len(str(i[1]))) * " " + str(i[1]) + 2 * " " + str(i[2]) + \
                 (4 - len(str(i[2]))) * " " + str(i[3]) + " " + str(i[4]) + " " * (4- len(str(i[5]))) + \
                 str(i[5]) + (11 - len(str(i[6]))) * " " + str(i[6]) + (8 - len(str(i[7]))) * " " + str(i[7]) + (8 - len(str(i[8]))) * \
                 " " + str(i[8]) + (6 - (len(str(i[9])))) * " " + str(i[9]) + (6 - (len(str(i[10])))) * " " \
                 + str(i[10]) + (12 - len(str(i[11]))) * " " + str(i[11])

        outfile.write(outstr + "\n")

def renumber_pdb(pdb, mhca, mhcb, peptide, tcra, tcrb, pdb_name):

    tcr_pmhc = [mhca, mhcb, peptide, tcra, tcrb]



    atoms = []

    for annot in tcr_pmhc:

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

    before = _get_before_atoms(pdb)
    after = _get_after_atoms(pdb)

    full_set = before + atoms + after

    list_to_pdb_spec(atoms, pdb_name + "_numbered.pdb")

