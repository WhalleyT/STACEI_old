three2one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

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