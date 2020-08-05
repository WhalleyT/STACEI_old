import swalign
import sys
import warnings

from . import generic_functions

#GLOBAL PARAMETERS
#PLEASE BE CAREFUL OF THEM
MATCH = 2
MISMATCH = -1
SW_SCORE = swalign.NucleotideScoringMatrix(MATCH, MISMATCH)
SMITH_WATERMAN = swalign.LocalAlignment(SW_SCORE)



def grab_VDJ_from_ANARCI(anarci):
    gene_keys = {}
    

    with open(anarci) as f:
        tra_found, trb_found, trd_found = False, False, False
        for line in f:
            split_line = line.split("|")
            
            if "TRB" in line:
                gene_keys["TRBV"] = split_line[2]
                gene_keys["TRBJ"] = split_line[4]
                gene_keys["B_species"] = split_line[1]
                gene_keys["B_score"] = split_line[3]
                trb_found = True
            elif "TRA" in line:
                gene_keys["TRAV"] = split_line[2]
                gene_keys["TRAJ"] = split_line[4]
                gene_keys["A_species"] = split_line[1]
                gene_keys["A_score"] = split_line[3]
                tra_found = True
            elif "TRD" in line:
                gene_keys["TRAV"] = split_line[2]
                gene_keys["TRAJ"] = split_line[4]
                gene_keys["A_species"] = split_line[1]
                gene_keys["A_score"] = split_line[3]
                trd_found = True
                tra_found = True

        if trd_found:
            warnings.warn("Warning. The ANARCI file annotated the TCRa as TCRd, this seems to happen on some "
                          "occasions, the numbering is still correct, but pay specific attention to it.")

        if not tra_found or not trb_found:
            sys.exit("Error, both the TCRa and TCRb were not found.")
    return gene_keys


def cheatsblastp(query, imgt_seq):
    alignment = SMITH_WATERMAN.align(query, imgt_seq)

    return alignment.r_pos, alignment.r_end, alignment.cigar

def grab_cdr_loops(pair, start, end, loop_string):
    nums = [i[0] for i in pair]
    
    #get start and end
    start = nums.index(str(start))
    end = nums.index(str(end)) + 1

    chunk = pair[start:end]
    
    #get the seqs and nums
    seq = "".join([i[1] for i in chunk])
    num = ",".join([i[0] for i in chunk])
    return loop_string + '="' + seq + '"[' + num + "]" 


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


def anarci_to_imgt_fasta(pdb, tcra, tcrb, peptide, mhca, mhcb, anarci, 
                        name, original_fasta, out_fasta):

    alpha_pdb_num, beta_pdb_num, alpha_pdb_seq, beta_pdb_seq = generic_functions.grab_pdb_numbers(pdb, tcra, tcrb)

    #now get VDJ usage
    gene_keys = grab_VDJ_from_ANARCI(anarci)

    #get start and end
    alpha = list(zip(alpha_pdb_num, alpha_pdb_seq))
    beta = list(zip(beta_pdb_num, beta_pdb_seq))

    for i in alpha:
        num = i[0]
        let = i[1]

        if num == "104":
            cys_alpha = let

    for i in beta:
        num = i[0]
        let = i[1]

        if num == "104":
            cys_beta = let
    #CDR1
    cdr1a = grab_cdr_loops(alpha, 27, 39, "CDR1a")
    cdr1b = grab_cdr_loops(beta, 27, 39, "CDR1b")   
    
    #CDR2
    cdr2a = grab_cdr_loops(alpha, 56, 65, "CDR2a")
    cdr2b = grab_cdr_loops(beta, 56, 65, "CDR2b")    
    
    #CDR3
    cdr3a = grab_cdr_loops(alpha, 104, 118, "CDR3a")
    cdr3b = grab_cdr_loops(beta, 104, 118, "CDR3b")

    #framework
    fwa = grab_cdr_loops(alpha, 81, 86, "FWa")
    fwb = grab_cdr_loops(beta, 79, 88, "FWb")

    mhca_name = ">" + name + "|" + mhca + "|MHCA|"
    mhcb_name = ">" + name + "|" + mhcb + "|MHCB|"
    pep_name = ">" + name + "|" + peptide + "|peptide|"
    
    tcra_name = ">" + name + "|"  + gene_keys["A_species"] + "|" + tcra + "|TCRA|" + gene_keys["TRAV"] + "|" + gene_keys["A_score"] + "|" + cdr1a + \
        "|" + cdr2a + "|" + fwa + "|" + cdr3a + '|Cys1b="C"[23]|Cys2b="' + cys_alpha + '"[104]|'

    tcrb_name = ">" + name + "|" + gene_keys["B_species"] + "|" + tcrb + "|TCRB|" + gene_keys["TRBV"] + "|" + gene_keys["B_score"] + "|" + cdr1b + \
        "|" + cdr2b + "|" + fwb + "|" + cdr3b + '|Cys1b="C"[23]|Cys2b="' + cys_beta + '"[104]|'
    
    names = [mhca_name, mhcb_name, pep_name, tcra_name, tcrb_name]
    seqs = []

    with open(original_fasta) as f:
        for name, seq in read_fasta(f):
            seqs.append(seq)
    
    outfile = open(out_fasta, "w")

    for name, seq in zip(names, seqs):
        outfile.write(name + "\n" + seq + "\n")
    
    print("Gene usage is as follows:\n\t-%s\n\t-%s\n\t-%s\n\t-%s" %(gene_keys["TRAV"], gene_keys["TRAJ"], gene_keys["TRBV"], gene_keys["TRBJ"]))
    return gene_keys