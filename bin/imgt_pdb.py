import re

def _extract_cdrs(tcr, names, chain):
    cdrs = ["CDR1", "CDR2", "CDR3", "FW"]
    cdrs = [e + chain for e in cdrs]
    nums = []
    lets = []
    for name in names:
        if tcr in name:
            split_name = name.split("|")
            for loop in cdrs:
                for element in split_name:
                    if loop in element:
                        letters = element.split('"')[1]
                        numbers = element.split("[")[1][:-1]
                        numbers = numbers.split(",")
                        nums.append(numbers)
                        lets.append(letters)
                        break
    return nums, lets


def _starts_with_met(pdb, alpha_chain, beta_chain):
    alpha_met = False
    beta_met = False

    found_a = False
    found_b = False

    for line in pdb:
        if line.startswith("ATOM"):
            if line.split()[4] is alpha_chain:
                if found_a is False:
                    if line.split()[3] is "MET":
                        alpha_met = True
                found_a = True

            if line.split()[4] is beta_chain:
                if found_b is False:
                    if line.split()[3] is "MET":
                        beta_met = True
                found_b = True

            if found_a and found_b:
                break
    return alpha_met, beta_met


def _get_start_offset(seq, start):
    if seq.startswith("."):
        for element in seq:
            if element is ".":
                start += 1
            else:
                break
    return start


def _write_indexes(nums, letters, start, chain_seq):
    # Reorder nums
    nums = [nums[0], nums[1], nums[3], nums[2]]

    # print len(nums)
    end = len(letters)
    indexes = []
    between = []


    for index, vector in enumerate(nums):
        # print vector
        if index > 0:
            between.append(range(int(previous_vector[-1]) + 1, int(vector[0])))
        previous_vector = vector

    start_list = range(start, int(nums[0][0]))
    #print start_list
    list_minus_end = start_list + nums[0] + between[0] + nums[1] + between[1] + nums[2] + between[2] + nums[3]
    list_minus_end = [int(e) for e in list_minus_end]

    indexes = list_minus_end + range(list_minus_end[-1] + 1, len(chain_seq))
    # print indexes
    return indexes

def _renumber(pdb, alpha_nums, alpha_letters, beta_nums, beta_letters, alpha_seq, beta_seq, achain, bchain):
    alpha_met, beta_met = _starts_with_met(pdb, achain, bchain)

    if alpha_met:
        alpha_start = 0
    else:
        alpha_start = 0

    if beta_met:
        beta_start = 0
    else:
        beta_start = 0

    alpha_start = _get_start_offset(alpha_seq, alpha_start)
    beta_start = _get_start_offset(beta_seq, beta_start)

    alpha_seq = alpha_seq.replace(".", "")
    beta_seq = beta_seq.replace(".", "")

    alpha_range = _write_indexes(alpha_nums, alpha_letters, alpha_start, alpha_seq)
    beta_range = _write_indexes(beta_nums, beta_letters, beta_start, beta_seq)
    #print alpha_range
    # print len(alpha_range)
    # print len(alpha_seq)
    return alpha_range, beta_range


def _get_2dim_list(pdb_list):

    first = True

    amino_acid = []
    chains = []
    for line in pdb_list:

        residue = line.split()[5]
        if first:
            counter = residue
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


def _rewrite_pdb(pdb_obj, a_range, b_range, achain, bchain):
    pdb_atom = [x for x in pdb_obj if x.startswith("ATOM")]
    pdb_alpha = [x for x in pdb_atom if achain == x.split()[4]]
    pdb_beta = [x for x in pdb_atom if bchain == x.split()[4]]

    pdb_2d_a = _get_2dim_list(pdb_alpha)
    pdb_2d_b = _get_2dim_list(pdb_beta)

    a_diff = len(pdb_2d_a) - len(a_range)
    b_diff = len(pdb_2d_b) - len(b_range)

    a_add = range(a_range[-1] + 1, a_range[-1] + 1 + a_diff)
    b_add = range(b_range[-1] + 1, b_range[-1] + 1 + b_diff)

    a_range = a_range + a_add
    b_range = b_range + b_add

    pdb_tcr = pdb_2d_a + pdb_2d_b
    ranges = a_range + b_range

    new = []
    length = 81

    for amino, replacement in zip(pdb_tcr, ranges):
        for atom in amino:
            split = re.split(r'(\s+)', atom)
            difference = len(str(replacement)) - len(split[10])
            offset = " " * difference
            split[10] = offset + str(replacement)
            line = ''.join(split)
            new.append(line)

    return new


def imgt_number_pdb(pdb_fasta, pdb_file, achain, bchain):
    names = []
    seqs = []
    pdb = []

    with open(pdb_fasta) as f:
        for line in f:
            if line.startswith(">"):
                names.append(line)
            else:
                seqs.append(line)

    with open(pdb_file) as p:
        for line in p:
            pdb.append(line)

    a_nums, a_lets = _extract_cdrs("|TCRA|", names, "a")
    b_nums, b_lets = _extract_cdrs("|TCRB|", names, "b")
    """
    for i,j in zip(a_nums, a_lets):
        print i
        print j
    """

    a_range, b_range = _renumber(pdb, a_nums, a_lets, b_nums, b_lets, seqs[3], seqs[4], achain, bchain)
    new_tcr = _rewrite_pdb(pdb, a_range, b_range, achain, bchain)

    # Get before and after our new TCR chains
    idx = 0
    pdb1 = []
    out = []
    for item in pdb:
        if item.startswith("ATOM") and item.split()[4] == achain:
            item = new_tcr[idx]
            idx += 1
        pdb1.append(item)

    for item in pdb1:
        if item.startswith("ATOM") and item.split()[4] == bchain:
            item = new_tcr[idx]
            idx += 1
        out.append(item)
    """
    outfile = open(pdb_file.rsplit(".", 1)[0] + "_imgt.pdb", "w")
    for o in out:
        outfile.write(o)
    """

def realign_pdb(pdb):
    before = 4
    total = 10
    after = 52

    arr = []
    with open(pdb) as file:
        for line in file:
            if line.startswith("ATOM"):
                split_line = re.split(r'(\s+)', line)
                split_line[10] = split_line[10].strip()
                before_len = len(split_line[9] + split_line[10])
                after_len = len("".join(split_line[11:]))
                if before_len != before:
                    before_diff = before - before_len
                    if before_diff > 0:
                        spaces = " " * before_diff
                        split_line[9] =  spaces + split_line[9]
                    if before_diff < 0:
                        split_line[9] = split_line[9][abs(before_diff):]
                if after_len != after:
                    after_diff = after - after_len
                    if after_diff > 0:
                        split_line[11] = split_line[11][after_diff:]
                    if after_diff < 0:
                        spaces = " " * after_diff
                        split_line[11] = spaces + split_line[11]
                line = ''.join(split_line)
            arr.append(line)

    out = open(pdb, "w")
    for line in arr:
        out.write(line)

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


def iterate_tcr(tcr, nums):
    for i, j in enumerate(tcr):
        print i,j

def find_start(tcr):
    if tcr[0] == "M":
        return -1
    else:
        return 0


def pop_list(numbers, letters):

    new_nums = []
    new_lets = []

    for i,j in zip(numbers, letters):
        if j != ".":
            new_nums.append(i)
            new_lets.append(j)

    return new_nums, new_lets


def get_pdb_list(pdb):
    l = []
    with open(pdb) as inf:
        for line in inf:
            l.append(line)

    return l


def renumber_chain(pdb, chain, chain_pairs):
    chain_iter = iter(chain_pairs)

    previous = 9999
    start = True

    out = []

    for line in pdb:
        if line.startswith("ATOM"):
            split = line.split()

            start_not_met = start and split[3] == "MET"
            #print start_not_met
            if start_not_met == False:
                if split[4] == chain:
                    #print "right chain"
                    current = int(split[5])

                    if current != previous:
                        #print "not previous"
                        start = False
                        chain_pair = chain_iter.next()
                        #print chain_pair

                    #print split
                    split[5] = str(chain_pair[0])
                    #print split

                    previous = current

                out.append("   ".join(split) + "\n")
    return out



def write_atoms(pdb_list, name):
    out = open(name, "w")

    for i in pdb_list:
        if i.startswith("ATOM"):
            out.write(i)


def get_pMHC(letters, pdb):
    out = []
    with open(pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                if line.split()[4] in letters:
                    out.append(line)
    return out


def renumber_ANARCI(fasta, pdb, tcra, tcrb, peptide, mhca, mhcb):
    outfile = pdb.rsplit(".", 1)[0] + "_imgt.pdb"

    tcra_seq = ""
    tcrb_seq = ""
    with open(fasta) as f:
        for name, seq in read_fasta(f):
            if name.split("|")[1] == tcra:
                tcra_seq = seq
            if name.split("|")[1] == tcrb:
                tcrb_seq = seq


    a_start = find_start(tcra_seq)
    b_start = find_start(tcrb_seq)

    alpha_range = range(a_start, len(tcra_seq))
    beta_range = range(b_start, len(tcrb_seq))

    a_nums, a_lets = pop_list(alpha_range, tcra_seq)
    b_nums, b_lets = pop_list(beta_range, tcrb_seq)

    pdb_list = get_pdb_list(pdb)

    a_pair = map(lambda x,y:(x,y), a_nums, a_lets)
    a_pair = a_pair + [(a_nums[-1] + 1, "X")]
    b_pair = map(lambda x,y:(x,y), b_nums, b_lets)
    b_pair = b_pair +  [(b_nums[-1] + 1, "X")]

    pdb_list = renumber_chain(pdb_list, tcra, a_pair)
    pdb_list = renumber_chain(pdb_list, tcrb, b_pair)

    #pmhc = get_pMHC([peptide, mhca, mhcb], pdb)
    #pdb_list = pmhc + pdb_list

    list_to_pdb_spec(pdb_list, outfile)


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