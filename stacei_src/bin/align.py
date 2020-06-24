import swalign


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


def get_before_after_sequence(before_file, after_file, chain):
    with open(before_file) as f:
        for name, seq in read_fasta(f):
            if name.split("|")[1] == chain:
                before = seq

    with open(after_file) as f:
        for name, seq in read_fasta(f):
            if name.split("|")[1] == chain:
                after = seq

    return before, after

def calculate_missing_start(old, new):
    """
    This function finds the start of the alignment for
    the new IMGT numbered sequence, in essence it jumps
    past the sequences lost in the HMM
    """

    MATCH = 2
    MISMATCH = -1
    SCORE = swalign.NucleotideScoringMatrix(MATCH, MISMATCH)

    sw = swalign.LocalAlignment(SCORE)

    new = new.replace(".", "")

    alignment = sw.align(old, new)
    offset = alignment.r_pos

    return offset

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
    start = None

    for line in pdb_list:
        if chain == line[21]:
            if started is False:
                start = int(line[22:27].strip())
                started = True

            number = int(line[22:27].strip())
            numbers.add(number)
    return len(numbers), start


def renumber(pdb_list, chain, residue_nums, start, offset):
    """
    Takes a pdb list file, a chain and a
    corresponding list of residues to replace
    in the PDB file
    """
    occurences = 0

    previous = start

    res_iter = iter(residue_nums)
    new_residue = next(res_iter)

    out = []
    residues = set()

    for line in pdb_list:
        if chain == line[21]:
            current = int(line[22:27].strip())

            residues.add(current)
            if len(residues) > offset:
                occurences += 1
                if current == previous or occurences == 1: #such that if we are skipping due to the offset we override the difference in residues
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
                    indexes = list(range(1, len(seq) + 1))
    return indexes


def imgt_clean_pdb(pdb_file, chains, new_fasta, old_fasta, outname):
    pdb = []
    with open(pdb_file) as f:
        for line in f:
            pdb.append(line)

    pdb = remove_conformations(pdb)
    pdb_2D = []


    for chain in chains:
        after, before = get_before_after_sequence(new_fasta, old_fasta, chain)
        offset = calculate_missing_start(before, after)

        size, start = get_default_residues(pdb, chain)
        resi = get_fasta_annot(new_fasta, chain)

        subpdb = renumber(pdb, chain, resi, start, offset)
        pdb_2D.append(subpdb)

        outfile = open(outname, "w")

    for i in pdb_2D:
        for j in i:
            outfile.write(j)

imgt_clean_pdb("3PL6_WT_nan_Y48A_B_WT_numbered.pdb", ["A", "B", "C", "D", "E"],
               "3PL6_WT_nan_Y48A_B_WT_ANARCI_IMGT_annotated.fasta",
               "3PL6_WT_nan_Y48A_B_WT/FASTAs/3PL6_WT_nan_Y48A_B_WT.fasta", "test.pdb")
