import subprocess
import pandas as pd
import os


def call_pisa(pdb, session):
    print("pisa %s -analyse %s" % (session, pdb))
    subprocess.call("pisa %s -analyse %s" % (session, pdb), shell=True)


def is_int(input):
    try:
        num = int(input)
    except ValueError:
        return False
    return True


def retrieve_interface(pairing):
    with open("session.txt") as f:
        for line in f:
            line = line.split()
            if line:
                if is_int(line[0]):
                    id = line[1]
                    from_chain = line[3]
                    to_chain = line[5]

                    if from_chain in pairing and to_chain in pairing:
                        if from_chain == pairing[1]:
                            structure = "Struct. 1"
                        else:
                            structure = "Struct. 2"

                        return id, structure

    return None, None


def extract_pisa(session, chain, all_chains, outfile, annotation_dictionary):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    subprocess.call("pisa %s -list interfaces > session.txt" % session, shell=True)
    print("pisa %s -list interfaces > session.txt" % session)

    # create a list of the other chains so we can collect interfaces
    from_chain = chain
    to_chain = list(set(all_chains) - set(chain))

    data = []

    for c in to_chain:
        tchain = c
        fchain = chain

        interface, structure = retrieve_interface([tchain, fchain])

        if interface is not None:
            donor = annotation_dictionary[fchain]
            acceptor = annotation_dictionary[tchain]

            subprocess.call("pisa %s -detail interfaces %s > interface.txt" % (session, interface), shell=True)

            with open("interface.txt") as f:
                chain_start, interface_start, end_line = False, False, False

                for line in f:

                    if chain_start and interface_start and not end_line:
                        if "-----'-'------------'--'" in line:
                            end_line = True
                        flat_line = line
                        line = line.split("|")

                        if len(line) > 1:
                            amino = d[line[2].split()[0].split(":")[-1]]
                            resi = line[2].split()[1]
                            asa = float(flat_line.split()[-3])
                            bsa = float(flat_line.split()[-2])

                            if asa > 0.0 and bsa > 0.0:
                                availabilty = (asa - bsa) / asa * 100
                            elif bsa == 0.0 and asa > 0.0:
                                availabilty = 100

                            data.append([donor, acceptor, resi, amino, asa, bsa, availabilty])

                    if structure in line:
                        chain_start = True

                    if "Interfacing Residues:" in line:
                        interface_start = True

    data = pd.DataFrame(data, columns=["Donor", "Acceptor", "No", "Residue", "ASA", "BSA", "Availability"])
    data.to_csv(outfile, sep="\t")

    os.remove("session.txt")
    os.remove("interface.txt")