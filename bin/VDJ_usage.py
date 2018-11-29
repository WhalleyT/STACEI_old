import swalign
import subprocess
import sys
import re
import argparse

from Bio.SeqIO import convert as bio_convert
from Bio.SeqIO import parse as bio_parse

import bin.data.TCRgeneDictionary as gene_dict


# Global Parameters
MATCH = 2
MISMATCH = -1
SW_SCORE = swalign.NucleotideScoringMatrix(MATCH, MISMATCH)
SMITH_WATERMAN = swalign.LocalAlignment(SW_SCORE)

def chain_class(value):
    """acceptable 5 letter argument for chain class"""
    if len(value) != 5 or value.isalpha() == False:
        raise argparse.ArgumentTypeError("%s is an invalid argument. It must be a string of 5 characters" % value)
    return value


# left over from doing it by header, but may be useful
def check_aa_header(pdb_file):
    x = False
    with open(pdb_file) as pdb:
        for line in pdb:
            if "SEQRES" in line:
                x = True
                break
    return x


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


def missing_elements(L):
    L = map(int, L)
    start, end = L[0], L[-1]
    return sorted(set(xrange(start, end + 1)).difference(L))


def find_missing(chains, pdb):
    starts, ends, missings = ([] for i in range(3))

    for chain in chains:
        numbers = []
        with open(pdb, "r") as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM"):
                    splitline = line.split()
                    if chain == splitline[4].strip():
                        resnum = splitline[5]
                        resnum = re.sub("\D", "", resnum)
                        numbers.append(resnum)

        # get unique residues and re-sort them into ascending order
        unique_residues = sorted(list(set(numbers)), key=int)

        start = unique_residues[0]
        end = unique_residues[-1]

       # print "Start: %i" % int(start), "End: %i" % int(end)

        missing = missing_elements(unique_residues)

        starts.append(start)
        ends.append(end)
        missings.append(missing)

    return starts, ends, missings


def xToDot(sequence):
    return sequence.replace("X", ".")


def generate_fasta(file_name, full_file, chains):
    fasta_name = file_name + "/FASTAs/" + file_name + ".fasta"

    bio_convert(full_file, "pdb-atom", fasta_name, "fasta")

    names = []
    seqs = []
   # print names, seqs
    with open(fasta_name) as fp:
        for name, seq in read_fasta(fp):
            names.append(name)
            seqs.append(seq)

    out = open(file_name + "/FASTAs/" + file_name + ".fasta", 'w')

    counter = 0
    chain_annotations = ["MHCA", "MHCB", "peptide", "TCRA", "TCRB"]
    for chain in chains:
        for i in range(0, len(names)):
            if names[i].split(":")[1] == chain:
                out.write(">" + file_name + "|" + chain + "|" + chain_annotations[counter] + "|\n")
                out.write(seqs[i] + "\n")
                counter += 1

    fasta_path = file_name + "/FASTAs/" + file_name + ".fasta"
   # print "fasta written to file " + fasta_path
    return fasta_path


def first_annotation_fasta(chains_class, fasta_name, pdb_name, full_file):
    chains = chains_class.mhca + chains_class.mhcb + chains_class.peptide + chains_class.tcra + chains_class.tcrb
    names = []
    seqs = []

    _convert_pdb_to_fasta(full_file, fasta_name)

    with open(fasta_name) as fp:
        for name, seq in read_fasta(fp):
            names.append(name)
            seqs.append(xToDot(seq))
   # print len(seqs)

    out = open(pdb_name + "/FASTAs/" + pdb_name + ".fasta", 'w')

    # convert letter index to their numeric equivalent e.g. A == 1
    nums = []
    for letter in chains:
        for index, name in enumerate(names):
            if name[-1] == letter:
                nums.append(index)

    new_names = ["MHCA", "MHCB", "peptide", "TCRA", "TCRB"]
    idx = 0

    # get our start, end and missing residues
    starts, ends, missing = find_missing(chains, full_file)

    index = 0
    for num in nums:
        sequence = seqs[num]
        name = ">" + pdb_name + "|" + chains[index] + "|" + new_names[index] + "|\n"
        out.write(name + sequence + "\n")
        index += 1

   # print "fasta written to file " + fasta_name


def _convert_pdb_to_fasta(pdb_name, pdb_fa):
    convert_fail = False

    try:
        bio_convert(pdb_name, "pdb-atom", pdb_fa, "fasta")
    except:
        # Change this exception when we find a case where the error happens
        print "Bio.convert was unsuccessful, trying by SeqRes"
        convert_fail = True

    if convert_fail:
        try:
            bio_convert(pdb_name, "pdb-seqres", pdb_fa, "fasta")
        except:
            print "Converting by SeqRes was also unsuccessful"
            sys.exit()


def run_anarci(anarci_in, anarci_out):
    print anarci_in, anarci_out
    subprocess.call(["ANARCI", "--sequence", anarci_in, "--scheme", "imgt",
                     "--assign_germline", "--outfile", anarci_out])


def read_file(infile, file_type):
    if infile is None:
        sys.exit('No file to read')
    if infile.split('.')[-1].lower() != str(file_type):
        sys.exit("\nFile extension must be of type ." + str(file_type) + "\n")
    else:
       # print 'Reading file: ' + str(infile)
        return open(infile, "r")


def anarci_block_parser(infile):
    everything = []
    block = []
    for line in infile.readlines():
        if "//" not in line:
            block.append(line[:-1])
        else:
            everything.append(block)
            block = []

   # print '\n' + "Input file contains " + str(len(everything)) + ' blocks'
    return everything


def fasta_parser(file_name):
    entries = []
    for seq_record in bio_parse(file_name + ".fasta", "fasta"):
        entry = []
        ident = "".join(seq_record.id)
        sequence = "".join(seq_record.seq)
        entry.append(ident)
        entry.append(sequence)
        entries.append(entry)
    return entries


def unpack(entries):
    new_entries = []
    for entry in entries:
        indent = entry[0]
        sequence = entry[1]
        new_entry = []
        id_terms = indent.split("|")
        # print ID.split("|")
        for term in id_terms:
            if len(term) != 0:
                new_entry.append(term)
        new_entry.append(sequence)

        new_entries.append(new_entry)

    return new_entries


def strip_annotation_from_block(block):
    annot = []
    seq = []
    for line in block:
        if line[0] == "#":
            annot.append(line)
        else:
            seq.append(line)
    return annot, seq


def cheatsblastp(query, imgt_seq):
   # print "\nWhere does " + imgt_seq + " fit?"
   # print "\nQuerying against " + query + "\n"
    alignment = SMITH_WATERMAN.align(query, imgt_seq)
    cigar = alignment.cigar
    #alignment.dump()

   # print "\n----------------------------------------\n"

    return alignment.r_pos, alignment.r_end, cigar


def remove_dots(input_sequence):
    input_no_dots = re.split('[^a-zA-Z]', input_sequence)
    output = []
    for part in input_no_dots:
        if len(part) != 0:
            output.append(part)
   # print output
    return output


def concat_locations(locations):
    output = '['
    for loc in locations:
        output += str(loc)
        output += ","
    output = output[:-1]
    output += ']'
    return output


def concat_seq(query, locations):
    output = ''
    for loc in locations:
        output += query[loc]
    return output


def add_gene_usage(block_annot, fasta_entry):
    output = []
    seq = fasta_entry[-1]
    annots = fasta_entry[:-1]
    usage = block_annot[8][1:]
    usage_terms = usage.split("|")

    output += annots
    output += usage_terms[1:-2]
    output.append(seq)

    return output


def add_germline_annotations(tcr_header, tcr_imgt, annotation):
    output = ''
    start = gene_dict.TCRdict[annotation][0]
    end = gene_dict.TCRdict[annotation][-1] + 1
    tracker = start
    locations = []
    sequence = []

    for aa in tcr_imgt[start:end]:
        if aa.isalpha():
            locations.append(tracker)
            sequence.append(tcr_imgt[tracker])
        tracker += 1

    output += annotation + "="
    output += "\"" + concat_seq(tcr_imgt, locations) + "\""
    output += concat_locations(locations)

   # print output
    return output


def add_cdr3_annot(tcr_header, tcr_imgt, annotation, where_idx):

   # print "F at ", where_idx
    output = ''
    start = gene_dict.TCRdict[annotation][0]
    end = where_idx + 1
    tracker = start
    locations = []
    sequence = []

    for aa in tcr_imgt[start:end]:
        if aa.isalpha():
            locations.append(tracker)
            sequence.append(tcr_imgt[tracker])
        tracker += 1

    output += annotation + "="
    output += "\"" + concat_seq(tcr_imgt, locations) + "\""
    output += concat_locations(locations)

   # print output
    return output


def extract_gene(tcr_entry):
   # print tcr_entry
    trv = tcr_entry[4].rsplit("*")[0].replace("-", "_")
    trj = tcr_entry[6].rsplit("*")[0].replace("-", "_")
   # print trv, trj
    return trv, trj


def anarci2fasta(block_seq):
    output = '.'
    for line in block_seq:
        output += line[-1]
    output = output.replace("-", ".")
   # print "\n", output
    return output


def parse_anarci(chains, fasta, anarci):

    fasta_file = read_file(fasta, "fasta")
    anarci_file = read_file(anarci, "txt")

    fasta_file_name = fasta.rsplit('.', 1)[0]
    file_name = fasta_file_name.rsplit('/')[0]
   # print file_name
    fasta_entries = fasta_parser(fasta_file_name)
    fasta_entries = unpack(fasta_entries)

    anarci_data = anarci_block_parser(anarci_file)

    for block in anarci_data:
        if block[0].split("|")[1] == chains.tcra:
            tcra_block = block
        if block[0].split("|")[1] == chains.tcrb:
            tcrb_block = block
    for entry in fasta_entries:
        if entry[1] == chains.tcra:
            tcra_entry = entry
        if entry[1] == chains.tcrb:
            tcrb_entry = entry

    tcra_entry = add_gene_usage(tcra_block, tcra_entry)
    tcrb_entry = add_gene_usage(tcrb_block, tcrb_entry)

    tcr_a_annot, tcra_anarci_seq = strip_annotation_from_block(tcra_block)
    tcr_b_annot, tcrb_anarci_seq = strip_annotation_from_block(tcrb_block)

    tcra_imgt = anarci2fasta(tcra_anarci_seq)
    tcrb_imgt = anarci2fasta(tcrb_anarci_seq)

    tcra_header = []

    tcra_annots = ["CDR1a", "CDR2a", "FWa", "Cys1a", "Cys2a"]
    for annot in tcra_annots:
        tcra_header.append(add_germline_annotations(tcra_entry, tcra_imgt, annot))

    tcrb_header = []
    tcrb_annots = ["CDR1b", "CDR2b", "FWb", "Cys1b", "Cys2b"]
    for annot in tcrb_annots:
        tcrb_header.append(add_germline_annotations(tcrb_entry, tcrb_imgt, annot))

    try:
        trav = gene_dict.TRAVgenes[extract_gene(tcra_entry)[0]]
    except KeyError as e:
        print e
        print "Trying truncated version of gene"
        index = extract_gene(tcra_entry)[0]
        index = index.split("_")[0]
        trav = gene_dict.TRAVgenes[index]

    try:
        traj = gene_dict.TRAJgenes[extract_gene(tcra_entry)[1]]
    except KeyError as e:
        print e
        print "Trying truncated version of gene"
        index = extract_gene(tcra_entry)[0]
        index = index.split("_")[0]
        traj = gene_dict.TRAJgenes[index]

    try:
        trbv = gene_dict.TRBVgenes[extract_gene(tcrb_entry)[0]]
    except KeyError, e:
        print e
        print "Trying truncated version of gene"
        index = extract_gene(tcrb_entry)[0]
        index = index.split("_")[0]
        trbv = gene_dict.TRBVgenes[index]

    try:
        trbj = gene_dict.TRBJgenes[extract_gene(tcrb_entry)[1]]
    except KeyError, e:
        print e
        print "Trying truncated version of gene"
        index = extract_gene(tcrb_entry)[0]
        index = index.split("_")[0]
        trbj = gene_dict.TRBJgenes[index]

    # Find the end of the CDR3 #

    traj_find = traj[-11:]
    traj_f, traj_end, traj_find_cigar = cheatsblastp(tcra_imgt, traj_find)
    cdr3_a_header = add_cdr3_annot(tcra_header, tcra_imgt, "CDR3a", traj_f)

    trbj_find = trbj[-10:]
    trbj_f, trbj_end, trbj_find_cigar = cheatsblastp(tcrb_imgt, trbj_find)
    cdr3_b_header = add_cdr3_annot(tcrb_header, tcrb_imgt, "CDR3b", trbj_f)

    # Find TRC #
    traj_f, traj_end, traj_find_cigar = cheatsblastp(tcra_entry[-1], traj_find)
   # print traj_end
    trac = tcra_entry[-1][traj_end:]

    trbj_f, trbj_end, trbj_find_cigar = cheatsblastp(tcrb_entry[-1], trbj_find)
   # print trbj_end
    trbc = tcrb_entry[-1][trbj_end:]

    tcra_imgt_final = tcra_imgt + trac
    tcrb_imgt_final = tcrb_imgt + trbc

    tcra_final = tcra_entry[:6] + tcra_header[:3] + [cdr3_a_header] + tcra_header[3:] + [tcra_imgt_final]
    tcrb_final = tcrb_entry[:6] + tcrb_header[:3] + [cdr3_b_header] + tcrb_header[3:] + [tcrb_imgt_final]

    fasta_entries_new = []


    for entry in fasta_entries:
        if entry[1] == chains.tcra:
            fasta_entries_new.append(tcra_final)

        elif entry[1] == chains.tcrb:
            fasta_entries_new.append(tcrb_final)

        else:
            fasta_entries_new.append(entry)

   # print "Repacking FASTA headers!"
    fasta_entries_out = []
    for entries in fasta_entries_new:
        new_id = ">"
        for term in entries[:-1]:
            new_id += term+"|"
        new_entry = []
        new_entry.append(new_id)
        new_entry.append(entries[-1])
        fasta_entries_out.append(new_entry)

    # Repacakge FASTA file
    fasta_out_t = ""
    for entry in fasta_entries_out:
        for terms in entry:
            fasta_out_t += terms + "\n"

   # print "Repacking FASTA File Done!"

    # Write FASTA file

    fasta_out = open(file_name + "_ANARCI_IMGT_annotated.fasta", "w")
    fasta_out.write(fasta_out_t)
   # print "\nOutputting FASTA file as: \n"
   # print fasta_out_t
    fasta_out.close()

   # print('     ~  End of parseARANCI.py v0.1 BETA  ~')
