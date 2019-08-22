import swalign
import data.TCRgeneDictionary as gene_seq_dict

#GLOBAL PARAMETERS
#PLEASE BE CAREFUL OF THEM
three2one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
MATCH = 2
MISMATCH = -1
SW_SCORE = swalign.NucleotideScoringMatrix(MATCH, MISMATCH)
SMITH_WATERMAN = swalign.LocalAlignment(SW_SCORE)

def grab_pdb_numbers(pdb_file, alpha_chain, beta_chain):
    alpha_num = []
    beta_num = []

    alpha_seq = ""
    beta_seq = ""

    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21]

                if chain == alpha_chain:
                    amino_acid = three2one[line[17:20].strip()]
                    res_num = line[22:27].strip()
                    
                    if res_num not in alpha_num:
                        alpha_num.append(res_num)
                        alpha_seq += amino_acid

                if chain == beta_chain:
                    amino_acid = three2one[line[17:20].strip()]
                    res_num = line[22:27].strip()
                    
                    if res_num not in beta_num:
                        beta_num.append(res_num)
                        beta_seq += amino_acid
    
    return alpha_num, beta_num, alpha_seq, beta_seq


def grab_VDJ_from_ANARCI(anarci):
    gene_keys = {}
    

    with open(anarci) as f:
        for line in f:
            split_line = line.split("|")
            
            if "TRB" in line:
                gene_keys["TRBV"] = split_line[2]
                gene_keys["TRBJ"] = split_line[4]
                gene_keys["B_species"] = split_line[1]
                gene_keys["B_score"] = split_line[3]   

            if "TRA" in line:
                gene_keys["TRAV"] = split_line[2]
                gene_keys["TRAJ"] = split_line[4]
                gene_keys["A_species"] = split_line[1]
                gene_keys["A_score"] = split_line[3]  

    return gene_keys


def cheatsblastp(query, imgt_seq):
    alignment = SMITH_WATERMAN.align(query, imgt_seq)
    cigar = alignment.cigar
    alignment.dump()

    return alignment.r_pos, alignment.r_end, cigar

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

    alpha_pdb_num, beta_pdb_num, alpha_pdb_seq, beta_pdb_seq = grab_pdb_numbers(pdb, tcra, tcrb)

    #now get VDJ usage
    gene_keys = grab_VDJ_from_ANARCI(anarci)

    #get start and end
    alpha = zip(alpha_pdb_num, alpha_pdb_seq)
    beta = zip(beta_pdb_num, beta_pdb_seq)

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
    
    print "Gene usage is as follows:\n\t-%s\n\t-%s\n\t-%s\n\t-%s" %(gene_keys["TRAV"], gene_keys["TRAJ"], gene_keys["TRBV"], gene_keys["TRBJ"])
    return gene_keys["TRAV"], gene_keys["TRAJ"], gene_keys["TRBV"], gene_keys["TRBJ"]