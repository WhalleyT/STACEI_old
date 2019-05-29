import subprocess
import pandas as pd
def make_pmhc_pdb(infile, outfile, pmhc):
    writer = open(outfile, "w")

    with open(infile) as f:
        for line in f:
            if line.startswith("ATOM"):
                if line.split()[4] in pmhc:
                    writer.write(line)


def call_pisa(pdb, session):
    subprocess.call("pisa %s -analyse %s" % (session, pdb), shell=True)


def extract_pisa(session, chain, outfile):

    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    subprocess.call("pisa %s -list monomers > session.txt" % session, shell=True)

    with open("session.txt") as f:
        for line in f:
            if "    " + chain + "    " in line:
                chain = line.split()[0]
                break

    subprocess.call("pisa %s -detail monomers %s > %s" % (session, chain, outfile), shell=True)
    start = False
    data = []

    with open(outfile) as f:
        for line in f:

            if start:
                if len(line) > 1:
                    if line.strip()[0].isdigit():
                        line = line.split()
                        residue = d[line[2].split(":")[1]]
                        num = line[3]
                        asa = float(line[5])
                        bsa = float(line[6])

                        if asa > 0.0 and bsa > 0.0:
                            availabilty = (asa - bsa) / asa * 100
                        elif bsa == 0.0 and asa > 0.0:
                            availabilty = 100

                        data.append([residue, num, asa, bsa, availabilty])

            if "2. Residue Accessibility and Solvation Energy Effect" in line:
                start = True

            if "Total" in line and start:
                ASA = line.split()[2]
                BSA = line.split()[3]


    data = pd.DataFrame(data, columns = ["Residue", "No", "ASA", "BSA", "Availability"])
    data.to_csv(outfile, sep="\t")

    return BSA, ASA
