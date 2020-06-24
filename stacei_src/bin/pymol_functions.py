#!/usr/bin/env python2
import time
import csv
import pymol
import os


def show_residues(PDB_path, residues):
    """Show the TCR-pMHC complex, no added ones"""
    pymol.cmd.load(PDB_path)
    pymol.cmd.hide('lines')
    for i in residues:
        res_chain = 'chain ' + i
        pymol.cmd.select(res_chain)
        pymol.cmd.show('cartoon', 'sele')

    pymol.cmd.color('chromium')
    return


def colour_CDR_loops(TCRa_res, TCRb_res):
    pymol.cmd.select('CDR1a', 'chain %s & (resi 23:33)' % TCRa_res)
    pymol.cmd.select('CDR2a', 'chain %s & (resi 48:58)' % TCRa_res)
    pymol.cmd.select('CDR3a', 'chain %s & (resi 90:110)' % TCRa_res)
    pymol.cmd.select('CDR3afw', 'chain %s & (resi 64:74)' % TCRa_res)
    pymol.cmd.select('CDR1b', 'chain %s & (resi 23:33)' % TCRb_res)
    pymol.cmd.select('CDR2b', 'chain %s & (resi 46:56)' % TCRb_res)
    pymol.cmd.select('CDR3b', 'chain %s & (resi 90:110)' % TCRb_res)
    pymol.cmd.select('CDR3bfw', 'chain %s & (resi 64:74)' % TCRb_res)

    colours = ['smudge', 'teal', 'lightpink', 'dirtyviolet', 'wheat', 'deepsalmon', 'aquamarine', 'paleyellow']
    loops = ['CDR1a', 'CDR2a', 'CDR3a', 'CDR3afw', 'CDR1b', 'CDR2b', 'CDR3b', 'CDR3bfw']
    for i in range(0, len(loops) - 1):
        j = colours[i]
        k = loops[i]
        pymol.cmd.color('%s' % j, '%s' % k)

    pymol.cmd.label('chain %s & (resi 28) and name CB' % TCRa_res, '"CDR1a" ')
    pymol.cmd.label('chain %s & (resi 100) and name CB' % TCRa_res, '"CDR2a" ')
    pymol.cmd.label('chain %s & (resi 69) and name CB' % TCRa_res, '"CDR3a" ')
    pymol.cmd.label('chain %s & (resi 100) and name CB' % TCRa_res, '"CDR3afw" ')
    pymol.cmd.label('chain %s & (resi 28) and name CB' % TCRb_res, '"CDR1b" ')
    pymol.cmd.label('chain %s & (resi 51) and name CB' % TCRb_res, '"CDR2b" ')
    pymol.cmd.label('chain %s & (resi 100) and name CB' % TCRb_res, '"CDR3b" ')
    pymol.cmd.label('chain %s & (resi 69) and name CB' % TCRb_res, '"CDR3bfw" ')

    return


def write_CDR_file(TCRa_res, TCRb_res, cwd):
    pymol.cmd.select('TCRs', 'chain %s+%s' % (TCRa_res, TCRb_res))
    pymol.cmd.orient('TCRs')
    pymol.cmd.zoom('TCRs')
    pymol.cmd.disable('TCRs')
    pymol.cmd.disable('dots')
    png_out = cwd + '/' + 'CDR_loop'
    pymol.cmd.bg_color(color="white")
    pymol.cmd.hide('nonbonded')
    pymol.cmd.png('cdr_loops', width=0, height=0, dpi=-1, ray=1, quiet=0)


def get_peptide(PDB_path, peptide_res):
    pymol.cmd.delete('all')
    pymol.cmd.load(PDB_path)
    pymol.cmd.hide('lines')
    pymol.cmd.select('peptide', 'chain %s' % peptide_res)
    pymol.cmd.zoom('peptide')
    pymol.cmd.orient('peptide')
    pymol.cmd.show('sticks', 'peptide')
    pymol.cmd.bg_color(color="white")


def colour_peptide_buriedness(peplen, path, home):
    BSA_values = []
    # open up BSA file
    f = open(home + "/" + path + '/peptide_BSA.txt', 'r')
    for row in f:
        BSA_values.append(float(row))

    peptide_position = list(range(1, peplen + 1))
    numpep = len(peptide_position)
    res = 0
    while (res < numpep):
        residue = res + 1
        pymol.cmd.select('res%i' % residue, ' peptide & resi %i' % residue)
        pymol.cmd.alter('res%i' % residue, 'b = %f' % BSA_values[res])
        res += 1

    minBSA = min(BSA_values);
    maxBSA = max(BSA_values)

    # pymol.cmd.label('res1', '1')
    # pymol.cmd.label('res%i' % numpep, '%i' % numpep)
    pymol.cmd.spectrum('b', 'red_white_blue', minimum=minBSA, maximum=maxBSA)
    pymol.cmd.disable('peptide')
    pymol.cmd.disable('res%i' % residue)
    pymol.cmd.bg_color(color="white")
    pymol.cmd.disable('dots')
    pymol.cmd.hide('nonbonded')
    pymol.cmd.select('het', 'hetatm')
    pymol.cmd.hide('everything', 'het')
    pymol.cmd.deselect()

    pymol.cmd.png('peptide_buried_pyMol', width=0, height=0, dpi=-1, ray=1, quiet=0)

    while not os.path.exists("peptide_buried_pyMol.png"):
        time.sleep(1)
    return numpep

def peptide_in_pocket(MHCa, MHCb):
    pymol.cmd.select("MHCa", 'chain %s' %MHCa)
    pymol.cmd.show("cartoon", "MHCa")
    pymol.cmd.color("grey60", "MHCa")
    pymol.cmd.select("MHCb", 'chain %s' %MHCa)
    pymol.cmd.show("cartoon", "MHCb")
    pymol.cmd.color("grey80", "MHCb")
    pymol.cmd.deselect()

    """ Being moody on some pymol?
    pymol.cmd.set_view(0.274721891, 0.908298671, -0.315469623,
                       -0.743769169, 0.408674866, 0.528954268,
                       0.609373331, 0.089321382, 0.787836611,
                       0.000000000, 0.000000000,  -75.816986084,
                       33.574874878, -7.009893417, 50.194580078,
                       59.774688721,   91.859283447,  -20.000000000)
    """
    #pymol.cmd.fog(0)
    #pymol.cmd.depth(0)
    pymol.cmd.ray()
    pymol.cmd.png("peptide_groove_pymol", width=0, height=0, dpi=-1, ray=0, quiet=0)

    while not os.path.exists("peptide_groove_pymol.png"):
        time.sleep(1)



def load_all_TCR_pMHC_complexes(PDB_path, peptide_res, MHCa_res, MHCb_res, TCRa_res, TCRb_res):
    pymol.cmd.delete('all')
    pymol.cmd.load(PDB_path)
    pymol.cmd.hide('lines')
    # pymol.cmd.load('ccp4_ncont.py')
    pymol.cmd.select('pMHC', 'chain %s+%s+%s' % (peptide_res, MHCa_res, MHCb_res))
    pymol.cmd.select('MHC', 'chain %s+%s' % (MHCa_res, MHCb_res))
    pymol.cmd.select('peptide', 'chain %s' % (peptide_res))
    pymol.cmd.select('TCR', 'chain %s+%s' % (TCRa_res, TCRb_res))
    pymol.cmd.select('complex', 'chain %s+%s+%s+%s+%s' % (peptide_res, MHCa_res, MHCb_res, TCRa_res, TCRb_res))
    pymol.cmd.bg_color(color="white")
    return


def orientate_pMHC():
    pymol.cmd.show('cartoon', 'pMHC')
    pymol.cmd.zoom('pMHC')
    pymol.cmd.orient('pMHC')
    pymol.cmd.color('chromium')
    pymol.cmd.disable('complex')
    return


def colour_contacts(MHC_to_pep_contacts, TCR_to_pMHC_contacts, MHCa_res, MHCb_res, TCRa_res, TCRb_res, peptide_res):
    for i in MHC_to_pep_contacts:
        # 0,1,5,6
        pymol.cmd.color('deepsalmon', 'chain %s & (resi %s)' % (i[0], i[1]))
        pymol.cmd.color('smudge', 'chain %s & (resi %s)' % (i[5], i[6]))

    pymol.cmd.bg_color(color="white")
    pymol.cmd.disable('dots')
    pymol.cmd.hide('nonbonded')
    pymol.cmd.rotate("y", 270)
    time.sleep(0.5)
    pymol.cmd.png('peptide_to_MHC', width=0, height=0, dpi=-1, ray=1, quiet=0)
    time.sleep(0.5)

    pymol.cmd.hide('cartoon')
    pymol.cmd.color('chromium')
    pymol.cmd.show('cartoon', 'complex')

    for i in TCR_to_pMHC_contacts:
        if i[0] == MHCa_res or MHCb_res:
            pymol.cmd.color('deepsalmon', 'chain %s & (resi %s)' % (i[0], i[1]))
        if i[0] == peptide_res:
            pymol.cmd.color('aquamarine', 'chain %s & (resi %s)' % (i[0], i[1]))
        if i[5] == TCRa_res:
            pymol.cmd.color('smudge', 'chain %s & (resi %s)' % (i[5], i[6]))
        if i[5] == TCRb_res:
            pymol.cmd.color('teal', 'chain %s & (resi %s)' % (i[5], i[6]))

    pymol.cmd.bg_color(color="white")
    pymol.cmd.disable('dots')
    pymol.cmd.hide('nonbonded')
    time.sleep(0.5)
    pymol.cmd.png('TCR_to_pMHC', width=0, height=0, dpi=-1, ray=1, quiet=0)
    time.sleep(0.5)

    pymol.cmd.hide('cartoon')
    pymol.cmd.color('grey')
    pymol.cmd.show('cartoon', 'TCR')

    for i in TCR_to_pMHC_contacts:
        if i[5] == TCRa_res:
            if 23 <= int(i[6]) <= 33:
                pymol.cmd.color('smudge', 'chain %s & (resi %s)' % (i[5], i[6]))
            if 48 <= int(i[6]) <= 58:
                pymol.cmd.color('teal', 'chain %s & (resi %s)' % (i[5], i[6]))
            if 90 <= int(i[6]) <= 110:
                pymol.cmd.color('lightpink', 'chain %s & (resi %s)' % (i[5], i[6]))
            if 64 <= int(i[6]) <= 74:
                pymol.cmd.color('dirtyviolet', 'chain %s & (resi %s)' % (i[5], i[6]))

        if i[5] == TCRb_res:
            if 23 <= int(i[6]) <= 33:
                pymol.cmd.color('wheat', 'chain %s & (resi %s)' % (i[5], i[6]))
            if 46 <= int(i[6]) <= 56:
                pymol.cmd.color('deepsalmon', 'chain %s & (resi %s)' % (i[5], i[6]))
            if 90 <= int(i[6]) <= 110:
                pymol.cmd.color('aquamarine', 'chain %s & (resi %s)' % (i[5], i[6]))
            if 64 <= int(i[6]) <= 74:
                pymol.cmd.color('paleyellow', 'chain %s & (resi %s)' % (i[5], i[6]))

        pymol.cmd.zoom('TCR')
        pymol.cmd.orient('TCR')
        pymol.cmd.bg_color(color="white")
        pymol.cmd.disable('dots')
        pymol.cmd.disable('TCR')
        pymol.cmd.hide('nonbonded')
        pymol.cmd.rotate("y", 270)
        time.sleep(0.5)
        pymol.cmd.png('CDR_contacts', width=0, height=0, dpi=-1, ray=1, quiet=0)
        time.sleep(0.5)
        return


def find_peptide_start(pdb_file, pep_residue):
    start = 0
    with open(pdb_file) as pdb:
        for line in pdb:
            if line.startswith("ATOM   "):
                split_lines = line.split()
                if split_lines[4] == pep_residue:
                    start = split_lines[5]
                    break
    return int(start)


def colour_peptide_contacts(peplen, TCR_to_pMHC, MHC_to_pep, peptide_res, PDB_path, numpep, start):
    pep_mhc_contacts = [0] * peplen
    pep_tcr_contacts = [0] * peplen

    for line in open(TCR_to_pMHC):
        columns = line.split("\t")
        if columns[0] == peptide_res:
            res = int(columns[1])
            idx = res - start
            pep_tcr_contacts[idx] += 1

    for line in open(MHC_to_pep):
        columns = line.split("\t")
        if columns[5] == peptide_res:
            res = int(columns[6])
            idx = res - start
            pep_mhc_contacts[idx] += 1

    # reload pymol
    pymol.cmd.delete('all')
    pymol.cmd.load(PDB_path)
    pymol.cmd.hide('lines')
    pymol.cmd.select('peptide', 'chain %s' % peptide_res)
    pymol.cmd.zoom('peptide')
    pymol.cmd.orient('peptide')
    pymol.cmd.show('sticks', 'peptide')

    res = 0
    while (res < numpep):
        residue = res + 1
        pymol.cmd.select('res%i' % residue, ' peptide & resi %i' % residue)
        pymol.cmd.alter('res%i' % residue, 'b = %f' % pep_mhc_contacts[res])  # raise here, but it works... ?
        res += 1

    minCont = min(pep_mhc_contacts)
    maxCont = max(pep_mhc_contacts)

    pymol.cmd.spectrum('b', 'blue_white_red', minimum=minCont, maximum=maxCont)
    pymol.cmd.disable('peptide')
    pymol.cmd.disable('res%i' % residue)
    pymol.cmd.bg_color(color="white")
    pymol.cmd.disable('dots')
    pymol.cmd.hide('nonbonded')
    time.sleep(0.5)
    pymol.cmd.png('peptide_heat_MHC_pyMol', width=0, height=0, dpi=-1, ray=1, quiet=0)
    time.sleep(0.5)

    pymol.cmd.delete('all')
    pymol.cmd.load(PDB_path)
    pymol.cmd.hide('lines')
    pymol.cmd.select('peptide', 'chain %s' % peptide_res)
    pymol.cmd.zoom('peptide')
    pymol.cmd.orient('peptide')
    pymol.cmd.show('sticks', 'peptide')

    res = 0
    while (res < numpep):
        residue = res + 1
        pymol.cmd.select('res%i' % residue, ' peptide & resi %i' % residue)
        pymol.cmd.alter('res%i' % residue, 'b = %i' % pep_tcr_contacts[res])
        res += 1

    minCont = min(pep_tcr_contacts);
    maxCont = max(pep_tcr_contacts)

    pymol.cmd.spectrum('b', 'blue_white_red', minimum=minCont, maximum=maxCont)
    pymol.cmd.disable('peptide')
    pymol.cmd.disable('res%i' % residue)
    pymol.cmd.bg_color(color="white")
    pymol.cmd.disable('dots')
    pymol.cmd.hide('nonbonded')
    time.sleep(0.5)
    pymol.cmd.png('peptide_heat_TCR_pyMol', width=0, height=0, dpi=-1, ray=1, quiet=0)
    time.sleep(0.5)
    return


def find_CDR_loop(TCR_chain, residue, TCRa_residue, TCRb_residue):
    residue = int(residue)
    loop = ""
    if TCR_chain == TCRa_residue:
        if 23 <= residue <= 33:
            loop = "CDR1a"
        if 48 <= residue <= 58:
            loop = "CDR2a"
        if 90 <= residue <= 110:
            loop = "CDR3a"
        if 64 <= residue <= 74:
            loop = "CDR3afw"

    elif TCR_chain == TCRb_residue:
        if 23 <= residue <= 33:
            loop = "CDR1b"
        if 46 <= residue <= 56:
            loop = "CDR2b"
        if 90 <= residue <= 110:
            loop = "CDR3b"
        if 64 <= residue <= 74:
            loop = "CDR3bfw"
    else:
        loop = "None"

    return (loop)


def color_MHC(MHC_residue, MHC_contacts, TCR_loops):
    CDR_colours = ['forest', 'lightteal', 'darksalmon', 'splitpea', 'raspberry', 'grey50', 'deepblue', 'brown',
                   'paleyellow', 'red']
    CDR_loops = ['CDR1a', 'CDR2a', 'CDR3a', 'CDR3afw', 'CDR1b', 'CDR2b', 'CDR3b', 'CDR3bfw', 'None']

    for i in range(len(MHC_contacts)):
        residue = int(MHC_contacts[i])
        loop = TCR_loops[i]
        for j in range(len(CDR_loops)):
            if loop == CDR_loops[j]:
                color = CDR_colours[j]

                # now have color and residue
                pymol.cmd.select('res%i' % residue, 'chain %s & resi %i' % (MHC_residue, residue))
                pymol.cmd.color('%s' % color, 'res%i' % residue)
                pymol.cmd.disable('res%i' % residue)
    return


def colour_CDR_wrapper(pdb_name, MHCa_res, MHCb_res, peptide_res, TCRa_res, TCRb_res, PDB_path):
    parsed_contact_file = pdb_name + '_TCR_to_pMHC_contacts_out.txt'

    with open(parsed_contact_file) as f:
        reader = csv.reader(f, delimiter="\t")
        contacts = list(reader)

    MHCa_contacts = []
    MHCb_contacts = []
    TCR_to_MHCa = []
    TCR_to_MHCb = []

    nrow = len(contacts)

    for row in range(nrow):
        if contacts[row][0] == MHCa_res:
            # append MHC residue
            MHCa_contacts.append(contacts[row][1])
            # append TCR CDR loop to corresponding index
            MHCa_CDR = find_CDR_loop(TCRa_res, contacts[row][6], TCRa_res, TCRb_res)
            TCR_to_MHCa.append(MHCa_CDR)

        if contacts[row][0] == MHCb_res and MHC_class == 2:
            MHCb_contacts.append(contacts[row][1])
            MHCb_CDR = find_CDR_loop(TCRb_res, contacts[row][6], TCRa_res, TCRb_res)
            TCR_to_MHCb.append(MHCb_CDR)
    """
    now to actually plot it into pymol
    at the moment it doesn't do it by
    most common contacting CDR, just 
    last one
    """

    CDR_colours = ['forest', 'lightteal', 'darksalmon', 'splitpea', 'raspberry', 'grey50', 'deepblue', 'brown',
                   'paleyellow', 'red']
    CDR_loops = ['CDR1a', 'CDR2a', 'CDR3a', 'CDR3afw', 'CDR1b', 'CDR2b', 'CDR3b', 'CDR3bfw', 'None']

    pymol.cmd.delete('all')
    pymol.cmd.load(PDB_path)
    pymol.cmd.hide('lines')
    pymol.cmd.bg_color(color="white")
    pymol.cmd.select('MHCs', 'chain %s+%s' % (MHCa_res, MHCb_res))
    pymol.cmd.show('cartoon', 'MHCs')
    pymol.cmd.zoom('MHCs')
    pymol.cmd.orient('MHCs')
    pymol.cmd.disable('MHCs')
    pymol.cmd.color('chromium')
    # add MHC centre of mass if needed?
    pymol.cmd.select("MHCa", 'chain %s' % MHCa_res)
    pymol.cmd.select("MHCb", 'chain %s' % MHCb_res)

    MHCa_COM = pymol.cmd.centerofmass("MHCa")
    MHCb_COM = pymol.cmd.centerofmass("MHCb")

    """
    def float_2_string(MHC_coords):
        floatstr = ['{:.14f}'.format(x) for x in MHC_coords]
        floatstr = "[ " + ', '.join(str(e) for e in MHC_coords) + " ]"
        return floatstr

    MHCa_COM = float_2_string(MHCa_COM)
    MHCa_COM = float_2_string(MHCb_COM)
    """

    pymol.cmd.pseudoatom("COM1", pos=MHCa_COM)
    pymol.cmd.show('spheres', "COM1")
    pymol.cmd.select("COM1")
    pymol.cmd.color("red", "COM1")

    pymol.cmd.pseudoatom("COM2", pos=MHCb_COM)
    pymol.cmd.show('spheres', "COM2")
    pymol.cmd.select("COM2")
    pymol.cmd.color("green", "COM2")

    color_MHC(MHCa_res, MHCa_contacts, TCR_to_MHCa)
    color_MHC(MHCb_res, MHCb_contacts, TCR_to_MHCb)
    pymol.cmd.disable('dots')
    pymol.cmd.hide('nonbonded')

    time.sleep(0.5)
    pymol.cmd.png('CDR_on_MHC', width=0, height=0, dpi=-1, ray=1, quiet=0)
    time.sleep(0.5)
    return
