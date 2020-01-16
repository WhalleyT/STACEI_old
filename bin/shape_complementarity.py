import os
import warnings

def write_SC_pipe(MHCa_res, MHCb_res, peptide_res, TCRa_res, TCRb_res):
    
    f = open('sc_in.txt', 'w')
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
