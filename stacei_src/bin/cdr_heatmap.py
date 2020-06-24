import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re

class FastaInformation:
    def __init__(self, residues, names, sequences):
        self.residues = residues
        self.names = names
        self.sequences = sequences


class TCRpMHC:
    def __init__(self, tcra, tcrb, peptide, mhca, mhcb):
        self.tcra = tcra
        self.tcrb = tcrb
        self.peptide = peptide
        self.mhca = mhca
        self.mhcb = mhcb


def get_ranges(fasta):
    names = []
    seqs = []
    lens = []
    with open(fasta) as infile:
        for line in infile:
            if line.startswith(">"):
                names.append(line.split("|")[2])
            else:
                seqs.append(line)

    for name, seq in zip(names, seqs):
        if seq.startswith("M"):
            start = 0
            end = len(seq) - 1
        else:
            start = 1
            end = len(seq)

        lens.append(list(range(start, end)))

    return FastaInformation(lens, names, seqs)


def find_hi_low(text):

    text = text.split("=", 1)[-1]
    text = text.split(",")
    new = []
    for elem in text:
        new.append(int((re.sub("[^0-9]", "", elem))))

    return max(new), min(new)



def get_CDR_range(fasta):
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                if line.split("|")[2] == "TCRA":
                    alpha = line.split("|")
                if line.split("|")[2] == "TCRB":
                    beta = line.split("|")

    cdrs = ["CDR1", "CDR2", "CDR3", "FW"]
    loops = []
    sequences = []
    starts = []
    ends = []
    for cdr in cdrs:
        cdra = cdr + "a"
        for chunk in alpha:
            if chunk.startswith(cdra):
                loops.append(chunk)

        cdrb = cdr + "b"
        for chunk in beta:
            if chunk.startswith(cdrb):
                loops.append(chunk)

    for item in loops:
        sequences.append(re.findall(r'"([^"]*)"', item)[0])
        #print item
        h, l = find_hi_low(item)
        starts.append(l)
        ends.append(h)

    newseqs = []

    return sequences, starts, ends



def read_tsv(file):
    with open(file, 'r') as ins:
        return [[n for n in line.split("\t")] for line in ins]


def generate_map(fasta_file, tcr_file, mhc_file):
    fasta = get_ranges(fasta_file)
    tcr = read_tsv(tcr_file)
    mhc = read_tsv(mhc_file)
    #print mhc
    seqs, starts, ends = get_CDR_range(fasta_file)

    loops = ["CDR1a", "CDR1b", "CDR2a", "CDR2b", "CDR3a", "CDR3b", "FWa", "FWb"]

    resnum_list = []
    CDR_list = []
    pepnum_list = []
    contacts_list = []

    residues = []
    for line in mhc[1:]:
        #print line
        residues.append(line[9])

    pep_start = int(min(residues))
    pep_end = int(max(residues))
    peplen = pep_end - pep_start
    #todo get peptide length based on PDB

    pep_start = 1
    pep_end = 11
    peplen = pep_end - pep_start + 1
    peprange = list(range(pep_start, pep_end + 1))

    counter = 0
    for start, end in zip(starts, ends):
        cdrange = list(range(int(start), int(end) + 1))
        for num in cdrange:
            res_chunk = [num] * peplen
            pep_chunk = peprange
            con_chunk = [0] * peplen
            cdr_chunk = [loops[counter]] * peplen

            resnum_list += res_chunk
            CDR_list += cdr_chunk
            pepnum_list += pep_chunk
            contacts_list += con_chunk

        counter += 1

    for i in range(0, len(contacts_list)):
        cdr_residue = resnum_list[i]
        cdr_loop = CDR_list[i]
        pep_residue = pepnum_list[i]

        for line in tcr:
            if line[3] == cdr_loop:
                if int(line[2]) == cdr_residue:
                    if line[8] == "peptide":
                        if int(line[9]) == pep_residue:
                            contacts_list[i] += 1

    #print sum(contacts_list)
    contact_df = pd.DataFrame(
    {'Residue': resnum_list,
     'Peptide': pepnum_list,
     'CDR': CDR_list,
     'Contacts' : contacts_list
    })

    data = contact_df
    contact_df.to_csv("CDR_melt.txt", sep='\t')

    def facet(data, color):
        data = data.pivot(index="Residue", columns='Peptide', values='Contacts')
        g = sns.heatmap(data, cmap='Blues', cbar=False)

    with sns.plotting_context(font_scale=5.5):
        g = sns.FacetGrid(data, col="CDR", col_wrap=4, size=3, aspect=1)
        g = g.map_dataframe(facet)
        g.set_titles(col_template="{col_name}", fontweight='bold', fontsize=18)

    plt.savefig("CDR3_contacts.png", format="png", transparent=True)