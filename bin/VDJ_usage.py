import swalign
import subprocess
import sys

from Bio.SeqIO import convert as bio_convert

# Global Parameters
MATCH = 2
MISMATCH = -1
SW_SCORE = swalign.NucleotideScoringMatrix(MATCH, MISMATCH)
SMITH_WATERMAN = swalign.LocalAlignment(SW_SCORE)


def _read_fasta(fp):
    name = None
    seq = []

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


def _missing_elements(nums):
    nums = map(int, nums)
    start, end = nums[0], nums[-1]
    return sorted(set(xrange(start, end + 1)).difference(nums))


def _x_to_dot(sequence):
    return sequence.replace("X", ".")


def first_annotation_fasta(chains_class, fasta_name, pdb_name, full_file):
    chains = chains_class.mhca + chains_class.mhcb + chains_class.peptide + chains_class.tcra + chains_class.tcrb
    names = []
    seqs = []

    _convert_pdb_to_fasta(full_file, fasta_name)


    with open(fasta_name) as fp:
        for name, seq in _read_fasta(fp):
            names.append(name)
            seqs.append(_x_to_dot(seq))

    out = open(pdb_name + "/FASTAs/" + pdb_name + ".fasta", 'w')

    # convert letter index to their numeric equivalent e.g. A == 1
    nums = []
    for letter in chains:
        for index, name in enumerate(names):
            if name[-1] == letter:
                nums.append(index)

    new_names = ["MHCA", "MHCB", "peptide", "TCRA", "TCRB"]

    index = 0
    for num in nums:
        sequence = seqs[num]
        name = ">" + pdb_name + "|" + chains[index] + "|" + new_names[index] + "|\n"
        out.write(name + sequence + "\n")
        index += 1


def _convert_pdb_to_fasta(pdb_name, pdb_fa):
    convert_fail = False

    try:
        bio_convert(pdb_name, "pdb-atom", pdb_fa, "fasta")
    except:
        print "Bio.convert was unsuccessful, trying by SeqRes"
        convert_fail = True

    if convert_fail:
        try:
            bio_convert(pdb_name, "pdb-seqres", pdb_fa, "fasta")
        except:
            print "Converting by SeqRes was also unsuccessful"
            sys.exit()


def run_anarci(anarci_in, anarci_out):
    subprocess.call(["ANARCI", "--sequence", anarci_in, "--scheme", "imgt",
                     "--assign_germline", "--outfile", anarci_out])
