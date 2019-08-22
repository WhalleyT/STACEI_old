import sys
import re
import os
import warnings

def write_SC_pipe(MHCa_res, MHCb_res, peptide_res, TCRa_res, TCRb_res):
    
    f = open('sc_in.txt','w')
    sc_MHCa_line = 'chain ' + MHCa_res + '\n'
    sc_MHCb_line = 'chain ' + MHCb_res + '\n'
    sc_pep_line  = 'chain ' + peptide_res + '\n'
    sc_TCRa_line = 'chain ' + TCRa_res + '\n'
    sc_TCRb_line = 'chain ' + TCRb_res + '\n'

    f.write('molecule 1\n')
    f.write(sc_MHCa_line)
    f.write(sc_MHCb_line)
    f.write(sc_pep_line)
    f.write('molecule 2\n')
    f.write(sc_TCRa_line)
    f.write(sc_TCRb_line)
    f.write('end\n')
    f.close()

    return

def remove_conformations(conf_list, line):
            for i in range(len(conf_list)):
                if conf_list[i] in line:
                    line = re.sub(conf_list[i], amino_acids[i], line)
            return line

def split_joined_conformation_residue(line, amino_acids):
    for i in range(len(amino_acids)):
        if amino_acids[i] in line and " " + amino_acids[i] not in line:
            line = re.sub(amino_acids[i], " " + amino_acids[i], line)
    return line

def parse_first_conformation(line, first_conformation):
    for i in range(len(first_conformation)):
        if first_conformation[i] in line:
            index = line.find(first_conformation[i])
            end = index + len(first_conformation[i]) - 1
            substitute = line[index:end]
            line = re.sub(substitute, substitute[:-1] + " ", line)
    return line

def clean_SC_PDB(pdb_file, out_name, sc_path):
    solvents_to_remove  = ["NAG ", "NA ", "EMC ", "GOL ", "SO4 ", "PGE ", "7PE ", "TAM", "EPE", "EDO", "PG4" ] 
    first_conformation  = ["CD1 ", "CE1 ", "OE1 ", "NE1 ", "OG1 ", "OD1 ", "NH1 ", "NE1 ", "ND1 ", "CG1 ","CZ1 ", "CZ1 ", "CH1 "]
    second_conformation = ["CD2 ", "CE2 ", "OE2 ", "NE2 ", "OG2 ", "OD2 ", "NH2 ", "NE2 ", "ND2 ", "CG2 ","CZ2 ", "CZ3 ", "CH2 "]



    amino_acids = [ 'ALA', 'ARG', 'ASN', 'ASP', 'ASX',
                    'CYS', 'GLU', 'GLN', 'GLX', 'GLY',
                    'HIS', 'ILE', 'LEU', 'LYS', 'MET',
                    'PHE', 'PRO', 'SER', 'THR', 'TRP',
                    'TYR', 'VAL'    ]

    a_conformations = []
    b_conformations = []

    for i in amino_acids:
        a_conf = "A"+i
        b_conf = "B"+i
        a_conformations.append(a_conf)
        b_conformations.append(b_conf)

    file = open(pdb_file, "r")
    PDB = []

    for line in file:
        PDB.append(line)
    file.close()
    PDB_no_solution = []
    filt_name = sc_path + "/" + out_name + ".pdb"
    with open(filt_name, 'w') as output:
        for line in PDB:
                if line.startswith("ATOM"):
                    line = split_joined_conformation_residue(line, amino_acids)

                    
                if any(solvent in line for solvent in solvents_to_remove) or any(confor in line for confor in second_conformation) == True:
                    dump = line # just a dump to get if else to work
                else:
                    output.write(line)
    return filt_name

def get_ATOM(file):
    out = open("sc.pdb", "w")

    with open(file) as pdb:
        for line in pdb:
            if line.startswith("ATOM"):
                out.write(line)


def run_SC(filtered_name):
    cmdstring = "sc XYZIN %s < sc_in.txt > sc_out.txt" % filtered_name
    os.system(cmdstring)

    found = False
    imaginary = False

    with open('sc_out.txt', 'r') as f:
        for line in f:
            if 'Shape complementarity statistic Sc' in line:
                SC = line
                found = True
            
            if "imaginary contain" in line:
                imaginary = True

    f = open('sc.txt', 'w')
    if found:
        f.write(SC)
    else:
        f.write("NA")
        warnings.warn("SC statistic not found, logging as NA")
        if imaginary:
            warnings.warn("This is because an 'imaginary contain' has been reported, this generally means H atoms at the interface, please remove them manually if you want to ignore this")

    f.close()
    return


def remove_double_conformations(pdb):
    outlist = []

    with open(pdb) as file:
        for line in file:
            if line.startswith("ATOM"):
                strsplit = re.split(r'(\s+)', line)
                res = strsplit[6]
                if len(res) == 4:
                    res = list(strsplit[6])
                    res[0] = " "

                    res = "".join(res)
                    strsplit[6] = res
                    line = "".join(strsplit)
            outlist.append(line)

    out = open(pdb, "w")

    for line in outlist:
        out.write(line)
