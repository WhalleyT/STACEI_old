import os
import sys
import itertools
import re
from Bio.PDB import Selection as pdb_select
from Bio.PDB import PDBParser as pdb_parser
from Bio.PDB import is_aa as is_amino_acid


def read_file(infile, file_type):
    if infile is None:
        sys.exit('No file to read')
    if infile.split('.')[-1].lower() != str(file_type):
        sys.exit("\nFile extension must be of type ." + str(file_type) + "\n")
    else:
       #### print 'Reading file: ' + str(infile)
        return open(infile, "r")


def everything_parser(infile):
    everything_list = []
    for lines in infile.readlines():
        everything_list.append(lines[:-1])

   #### print '\n' + "Input file contains " + str(len(everything_list)) + ' lines'
    return everything_list


def ccp4_validation(infile):
    infile_yo = infile
    for lines in infile_yo:
        if "NCONT" in lines:
           #### print "Input file was validated as a ccp4 contact output file!"
            return True
    else:
        sys.exit("\nThis doesn't look like an NCONT output file. Please try again!")


def contact_parser(infile):
    infile_yo = infile
    contact_list = []
    for lines in infile_yo:
        if "]:" in lines:
            contact_list.append(lines)
   #### print str(len(contact_list)) + " contact lines detected"
    return contact_list


def do_contacts_match(list_everything, contact_list):
    everything_list_yo = list_everything
    contact_list_yo = contact_list
    total_contact_line = ''
    reference = 0
    analyte_len = len(contact_list_yo)
    for lines in everything_list_yo:
        if "Total" in lines and "contacts" in lines:
            total_contact_line += lines
    for s in total_contact_line.split():
        if s.isdigit():
            reference += int(s)
    if analyte_len == reference:
        #print "projectContact recognised the same number of contacts as reported by CCP4!"
        x = 1
    else:
        sys.exit("projectContacts did not find all the contacts in the contact file! Please report this error")


def fasta_parser(fasta):
    from Bio import SeqIO
    entries = []
    for seq_record in SeqIO.parse(fasta, "fasta"):
        entry = []
        seq_id = "".join(seq_record.id)
        sequence = "".join(seq_record.seq)
        entry.append(seq_id)
        entry.append(sequence)
        entries.append(entry)
    return entries


def depack_id(entries):
    new_entries = []
    for entry in entries:
        seq_id = entry[0]
        sequence = entry[1]
        new_entry = []
        id_terms = seq_id.split("|")
        #### print ID.split("|")
        for term in id_terms:
            if len(term) != 0:
                new_entry.append(term)
        new_entry.append(sequence)

        new_entries.append(new_entry)

    return new_entries


def find_locations(sub_entries):
    s = "=[]"
    locations = []
    for col in sub_entries:
        if s[0] in col and s[1] in col and s[2] in col:
            locations.append(col)
    return locations


def depack_locations(subentries):
    locations = []
    for col in subentries:
        location = [col.rsplit("=", 1)[0]]
        locationstring = (col.partition('[')[-1].rpartition(']')[0])
        location += map(int, locationstring.split(','))
        locations.append(location)
    return locations


def purge_cys_locations(locations):
    output = []
    for locs in locations:
        if "Cys" not in locs[0]:
            output.append(locs)
    return output


def line_needs_filled(line):
    if line.count(":") == 1:
        return True
    else:
        return False


def source_atom_grabber(line):
    split_line = line.split(':')
    return str(split_line[0] + ':')


def fill_line(contact_lines):
   #### print '\nFilling source atom lines...'
    new_lines = []
    position_counter = 0
    fill_counter = 0
    for lines in contact_lines:
        position_counter += 1
        if not line_needs_filled(lines):
            new_lines.append(lines)

        elif line_needs_filled(lines):
            fill_counter += 1
            track_back_counter = 2
            while str(source_atom_grabber(contact_lines[position_counter - track_back_counter])).count(' ') > 20:
                track_back_counter += 1
            prefix = source_atom_grabber(contact_lines[position_counter - track_back_counter])
            suffix = lines[28:]
            new_line = prefix + ' ' + suffix
            new_lines.append(new_line)
   #### print str(fill_counter) + " lines were filled with source atoms..."
    return new_lines


def contact_matrix_row(line, tcr_a_locations, tcr_b_locations, mhc_a_chain,
                       mhc_b_chain, peptide_chain, tcr_a_chain, tcr_b_chain):

    new_names = {mhc_a_chain: "MHCa",
                 mhc_b_chain: "MHCb",
                 peptide_chain: "peptide",
                 tcr_a_chain: "TCRa",
                 tcr_b_chain: "TCRb"
                 }
    s_new_name = new_names[line[4]]
    t_new_name = new_names[line[32]]

    annotation1 = ''
    if s_new_name == "TCRa":
        for loop in tcr_a_locations:
            for loc in loop[1:]:
                if loc == int(line[6:10]):
                    annotation1 = loop[0]

    if s_new_name == "TCRb":
        for loop in tcr_b_locations:
            for loc in loop[1:]:
                if loc == int(line[6:10]):
                    annotation1 = loop[0]
    annotation2 = ''
    contact_row = [line[4], s_new_name, int(line[6:10]), annotation1, line[11:14], line[19:22], line[24], line[32],
                   t_new_name, int(line[34:38]), annotation2, line[39:42], line[47:50], line[52], float(line[-4:])]

    return contact_row


def make_contact_matrix(contact_lines, tcr_a_locations, tcr_b_locations, mhc_a_chain,
                        mhc_b_chain, peptide_chain, tcr_a_chain, tcr_b_chain):
    omit_counter = 0
   #### print '\nCreating contact matrix...'
    contact_matrix = []
    for lines in contact_lines:
        if 'HOH' in lines:
            omit_counter += 1
        else:
            contact_matrix.append(contact_matrix_row(lines, tcr_a_locations, tcr_b_locations, mhc_a_chain,
                                                     mhc_b_chain, peptide_chain, tcr_a_chain, tcr_b_chain))
   #### print str(omit_counter) + ' contacts were omitted due to water "HOH"...'
    return sorted(contact_matrix, key=lambda items: (items[0], items[2]))


def is_hydrogen_bond(contact_row):
    hb_donors = ["ARGNE N", "ARGNH1N", "ARGNH2N", "ASNND2N", "GLNNE2N", "HISND1N", "HISNE2N", "LYSNZ N",
                 "SEROG O", "THROG1O", "TRPNE1N", "TYROH O", "ALAN  N", "ARGN  N", "ASNN  N", "ASPN  N", "CYSN  N"
                 "GLUN  N", "GLNN  N", "GLYN  N", "HISN  N", "ILEN  N", "LEUN  N", "LYSN  N", "METN  N", "PHEN  N",
                 "SERN  N", "THRN  N", "TRPN  N", "TYRN  N", "VALN  N"]

    hb_acceptors = ["ASNOD1O", "ASPOD1O", "ASPOD2O", "GLNOE1O", "GLUOE1O", "GLUOE2O", "HISND1O", "HISND2O",
                    "SEROG O", "THROG1O", "TYROH O", "ALAO  O", "ARGO  O", "ASNO  O", "ASPO  O", "CYSO  O", "GLUO  O",
                    "GLNO  O", "GLYO  O", "HISO  O", "ILEO  O", "LEUO  O", "LYSO  O", "METO  O", "PROO  O", "PHEO  O",
                    "SERO  O", "THRO  O", "TRPO  O", "TYRO  O", "VALO  O"]

    id_1 = contact_row[4] + contact_row[5] + contact_row[6]
    id_2 = contact_row[11] + contact_row[12] + contact_row[13]
    if id_1 in hb_donors and id_2 in hb_acceptors:
        if contact_row[14] <= 3.40:
            return True
    if id_2 in hb_donors and id_1 in hb_acceptors:
        if contact_row[14] <= 3.40:
            return True
    else:
        return False


def is_salt_bridge(contact_row):
    acid_atoms = ["GLUOE1O", "GLUOE2O", "ASPOD1O", "ASPOD2O"]
    base_atoms = ["LYSNZ N", "ARGNE N", "ARGNH1N", "ARGNH2N"]

    id_1 = contact_row[4] + contact_row[5] + contact_row[6]
    id_2 = contact_row[11] + contact_row[12] + contact_row[13]
    if id_1 in acid_atoms and id_2 in base_atoms:
        if contact_row[14] <= 3.40:
            return True
    if id_2 in acid_atoms and id_1 in base_atoms:
        if contact_row[14] <= 3.40:
            return True
    else:
        return False


def is_van_der_waals(contact_row):
    if contact_row[4] == contact_row[11]:
        if contact_row[5] == contact_row[12]:
            if contact_row[7] == contact_row[13]:
                # print "same atom hit"
                return False
    if contact_row[14] <= 4.0:
        return True


def bond_annotator(contacts_row):
    contacts_row.append('NO')

    if is_van_der_waals(contacts_row):
        contacts_row[15] = 'VW'
    if is_hydrogen_bond(contacts_row):
        contacts_row[15] = 'HB'
    if is_salt_bridge(contacts_row):
        contacts_row[15] = 'SB'


def annotate_all_wrapper(contact_matrix):
    contact_matrix_new = []
    omit_counter = 0
   #### print 'Anotating contacts...'
    for row in contact_matrix:
        bond_annotator(row)
        if row[15] == "NO":
            omit_counter += 1
        else:
            contact_matrix_new.append(row)

   #### print str(omit_counter) + ' contacts were omitted due to not meeting annotation criteria "NO"...'
    return contact_matrix_new


def clean_contacts(ncont, chains, fasta):
    in_file = read_file(ncont, "txt")

    in_file_name = ncont.rsplit('.', 1)[0]

    if type(in_file_name) != str:
        sys.exit('No file was loaded. Please view usage and provide a valid file to be processed')
    all_lines = everything_parser(in_file)
    contact_lines = contact_parser(all_lines)

    if fasta != "":
        fasta_entries = fasta_parser(fasta)
        fasta_entries = depack_id(fasta_entries)

        tcra = []
        for entry in fasta_entries:
            if "TCRA" in entry:
                tcra = entry
        tcr_a_locations = find_locations(tcra)
        tcr_a_locations = depack_locations(tcr_a_locations)
        tcr_a_locations = purge_cys_locations(tcr_a_locations)

        tcrb = []
        for entry in fasta_entries:
            if "TCRB" in entry:
                tcrb = entry
        tcr_b_locations = find_locations(tcrb)
        tcr_b_locations = depack_locations(tcr_b_locations)
        tcr_b_locations = purge_cys_locations(tcr_b_locations)

    else:
        tcr_a_locations = [""]
        tcr_b_locations = [""]

    # Sort chains
    mhc_a_chain, mhc_b_chain, peptide_chain, tcr_a_chain, tcr_b_chain = chains[0], chains[1], chains[2],\
                                                                        chains[3], chains[4]

    contact_lines = fill_line(contact_lines)
    contact_matrix = make_contact_matrix(contact_lines, tcr_a_locations, tcr_b_locations,
                                         mhc_a_chain, mhc_b_chain, peptide_chain, tcr_a_chain, tcr_b_chain)
    contact_matrix = annotate_all_wrapper(contact_matrix)

    out_file = open(str(in_file_name) + '_clean.txt', 'w')
    output2 = ''
    output2 += 'Donor_Chain_Letter \tDonor_Chain \tDonor_ResNum \tDonor_Annotation \tDonor_ResCode \tDonor_Position ' \
               '\tDonor_Atom \tAcceptor_Chain_Letter \tAcceptor_Chain \tAcceptor_ResNum \tAcceptor_Annotation ' \
               '\tAcceptor_ResCode \tAcceptor_Position \tAcceptor_Atom \tDistance \tType\n'

    for x in contact_matrix:
        for y in x:
            output2 += str(y) + '\t'
        output2 += '\n'

    out_file.write(output2)
   #### print '\nOutputted file: ' + str(in_file_name) + '_clean.txt'
   #### print('\n''     ~  End ProjectContacts.py v0.6      ~')


#def add_mhc_contact():


def run_ncont(pdb_name, mhca_res, mhcb_res, peptide_res, tcra_res, tcrb_res, filtered_name):
    pdb_mhc_pep_ncont = pdb_name + '_MHC_to_pep_contacts.txt'
    pdb_tcr_pmhc_ncont = pdb_name + '_TCR_to_pMHC_contacts.txt'


    pep_mhc_command = "ncont XYZIN %s <<eof > %s \n" \
                      "source /*/%s\n" \
                      "source /*/%s\n" \
                      "target /*/%s\n" \
                      "mindist 0.0\n" \
                      "maxdist 4.0" %(filtered_name, pdb_mhc_pep_ncont, mhca_res, mhcb_res, peptide_res)

    os.system(pep_mhc_command)

    tcr_pmhc_command = "ncont XYZIN %s <<eof > %s \n" \
                       "source /*/%s\n" \
                       "source /*/%s\n" \
                       "target /*/%s\n" \
                       "target /*/%s\n" \
                       "target /*/%s\n" \
                       "mindist 0.0\n" \
                       "maxdist 4.0" %(filtered_name, pdb_tcr_pmhc_ncont, tcra_res,
                                      tcrb_res, mhca_res, mhcb_res, peptide_res)
    os.system(tcr_pmhc_command)


def clean_pdb(pdb, name):
    out_name = open(name + "_clean_numbered_imgt.pdb", "w")

    with open(pdb) as infile:
        for line in infile:
            if line.startswith("ATOM"):
                out_name.write(line)


def data_stripper(contact_matrix):
    new_matrix = []
    for x in contact_matrix:
        y = []
        z = ''
        y.append(x[0])
        z += x[0]
        y.append(x[2])
        z += x[2]
        y.append(x[4])
        z += x[4]
        y.append(x[7])
        z += x[7]
        y.append(x[9])
        z += x[9]
        y.append(x[11])
        z += x[11]
        y.append(x[15])
        y.append(z)

        new_matrix.append(y)
    return new_matrix


def contact_grouper(contact_matrix):
    # Returns a three layer list: grouped contacts -> contacts -> contact parameters
   #### print "\nReducing contacts to single residue level.."

    result = []
    for key, group in itertools.groupby(contact_matrix, lambda x: x[7]):
        result.append(list(group))
    # print result
    return result


def force_bools(group):
    # Returns three bools: VW, HB, SB for each group
    vdw = False; HB = False; SB = False
    for x in group:
        if 'VW' in x[6]:
            vdw = True
        if 'HB' in x[6]:
            HB = True
        if 'SB' in x[6]:
            SB = True
    return vdw, HB, SB


def new_group(group):
    y = group[0]
    vdw, HB, SB = force_bools(group)
   #### print y[1]
    return [y[0], '', y[1], y[2], y[3], '', y[4], y[5], vdw, HB, SB]


def simple_contact_matrix(grouped_contact_matrix):
   #### print "\nGenerating single residue contact matrix.."

    y = []
    for x in grouped_contact_matrix:
        # print x
        y.append(new_group(x))
    return y


def three2one(seq):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    if len(seq) % 3 == 0:
        upper_seq = seq.upper()
        single_seq = ''
        for i in range(len(upper_seq) / 3):
            single_seq += d[upper_seq[3 * i:3 * i + 3]]
        return single_seq


def replace_3_to_1(aa):
    old = aa
    new = []
    for x in old:
        x[3] = three2one(x[3])
        x[7] = three2one(x[7])

        new.append(x)

    return new


def name_chains(line, where, chains):
    old_code = line[where]
    new_line = line
    if old_code == chains[0]:
        new_code = 'MHCA'
    if old_code == chains[1]:
        new_code = 'MHCB'
    if old_code == chains[2]:
        new_code = 'peptide'
    if old_code == chains[3]:
        new_code = 'TCRA'
    if old_code == chains[4]:
        new_code = 'TCRB'
    new_line[where] = new_code
    return new_line


def residue_only_contacts(infile, chains):
    in_file = read_file(infile, "txt")
    in_file_name = infile.rsplit('.', 1)[0]

    if type(in_file_name) != str:
        raise IOError('No file was loaded. Please view usage and provide a valid file to be processed')
    all_lines = everything_parser(in_file)
    contact_matrix = matrix_parser(all_lines)
    contact_matrix = contact_matrix[1:]
    contact_matrix = data_stripper(contact_matrix)

    grouped_contact_matrix = contact_grouper(contact_matrix)
    output_contact_matrix = simple_contact_matrix(grouped_contact_matrix)
    output_contact_matrix2 = output_contact_matrix[:]
    output_contact_matrix2 = replace_3_to_1(output_contact_matrix2)

    for x in output_contact_matrix:
        name_chains(x, 0, chains)
        name_chains(x, 4, chains)

    out_file = open(str(in_file_name) + '_residues.txt', 'w')

    output2 = ''
    output2 += 'Chain \tAnnotation \tResNum \tResCode \tChain \tAnnotation \tResNum \tResCode \tvdW \tHB \tSB \n'

    for x in output_contact_matrix:
       ### print x
        for y in x:
            output2 += str(y) + '\t'
        output2 += '\n'

    out_file.write(output2)
   ### print '\nOutputted file: ' + str(in_file_name) + '_residues.txt'
   ### print('\n''     ~  End ProjectContactsSimple.py v0.1      ~')


def matrix_parser(allLines):
    matrix = []
    for lines in allLines:
        split_line = lines.split('\t')
        matrix.append(split_line[:-1])
    return matrix


def add_index_code_seq(line):
    index = line[0] + line[2] + line[3]
    line.append(index)
    return line


def add_index_code_seq_all(seqMatrix):
    seq_matrix_new = []
    for x in seqMatrix:
        new_line = add_index_code_seq(x)
        seq_matrix_new.append(new_line)
    return seq_matrix_new


def add_index_code_contact(line):
    index = ''
    index += line[4] + line[6] + line[7]
   ### print line
    line.append(index)
    return line


def add_index_code_to_contact_all(con_matrix):
    conMatrixNew = []
    for x in con_matrix:
        new_line = add_index_code_seq(x)
        new_line = add_index_code_contact(x)
        conMatrixNew.append(new_line)
    return conMatrixNew


def pairContacts2parents2(contactMatrix,fullSeqMatrix):
    newFullSeqMatrix=[]
    for x in fullSeqMatrix:
        newSeqLine=[]
        hitCount=0
        for y in contactMatrix:
            if x[4] == y[11]:
                newSeqLine=[]
#                ### print "seq matrix in:"
#                ### print x
#                ### print "contact matrix goes in:"
#                ### print y
                hitCount+=1
                for terms in x:
                    newSeqLine.append(terms)
                newSeqLine.append(y[4])
                newSeqLine.append(y[5])
                newSeqLine.append(y[6])
                newSeqLine.append(y[7])
                newSeqLine.append(y[8])
                newSeqLine.append(y[9])
                newSeqLine.append(y[10])
                newSeqLine.append(y[11])
                newSeqLine.append(y[12])
                y[11],y[12] = "Done", "Done"
#                ### print "after done:"
#                ### print y
#                ### print "output line:"
#                ### print newSeqLine
#                ### print "\n\n"
                newFullSeqMatrix.append(newSeqLine)
        if hitCount == 0:
            newSeqLine=[]
            for terms in x:
                newSeqLine.append(terms)
            newSeqLine.append('')
            newSeqLine.append('')
            newSeqLine.append('')
            newSeqLine.append('')
            newSeqLine.append('')
            newSeqLine.append('')
            newSeqLine.append('')
            newSeqLine.append('')
            newSeqLine.append('')
            newFullSeqMatrix.append(newSeqLine)

    return newFullSeqMatrix


def fill_acceptor_annotation(seq_matrix):
    seq_matrix_search = seq_matrix
    seq_matrix_fix = seq_matrix
    for x in seq_matrix_fix:
        if len(x[5]) != 0:
            for y in seq_matrix_search:
                if y[4] == x[13]:
                    x[6] = y[1]
    return seq_matrix_fix


def remove_index(seq_matrix):
    new = []
    for line in seq_matrix:
        del line[4]
        del line[-2:]
        new.append(line)
    return new

def annotate_sequence_list(sequence_file, contact_file):
    seq_in_file = read_file(sequence_file, "txt")
    seq_in_file_name = sequence_file.rsplit('.', 1)[0]

    contact_in_file = read_file(contact_file, "txt")
    contactInFileName = contact_file.rsplit('.', 1)[0]
    # Openers #
    if type(seq_in_file_name) != str:
        raise IOError('No file was loaded. Please view usage and provide a valid file to be processed')
    if type(contactInFileName) != str:
        raise IOError('No file was loaded. Please view usage and provide a valid file to be processed')

    # Prep sequence file#
   ### print "\nParsing sequence information.."
    all_seq_lines = everything_parser(seq_in_file)
    sequence_list = matrix_parser(all_seq_lines)
    sequence_matrix = sequence_list[1:]
    sequence_matrix = add_index_code_seq_all(sequence_matrix)
   ### print "Done!\n"

    # Prep contact file#
   ### print "Parsing contact information.."
    allContactLines = everything_parser(contact_in_file)
    contactList = matrix_parser(allContactLines)
    contact_matrix = contactList[1:]
    contac_matrix = add_index_code_to_contact_all(contact_matrix)
   ### print "Done!\n"

    # Find pairs#
    paired_seq_matrix = pairContacts2parents2(contact_matrix, sequence_matrix)
    paired_seq_matrix = fill_acceptor_annotation(paired_seq_matrix)
    output_contact_matrix = remove_index(paired_seq_matrix)

    outFile = open(str(contactInFileName) + '_contacts_residues_full.txt', 'w')

    output2 = ''
    output2 += 'Chain \tAnnotation \tResNum \tResCode \tChain \tAnnotation \tResNum \tResCode \tvdW \tHB \tSB \n'

    for x in output_contact_matrix:
        for y in x:
            output2 += str(y) + '\t'
        output2 += '\n'

    outFile.write(output2)
   ### print '\nOutputted file: ' + str(contactInFileName) + '_contacts_residues_full.txt'
   ### print('\n''     ~  End fullResidueContacts.py v0.1      ~')

def pdb_to_contact_list():
    def readFile(FILE, fileType):
        if FILE == None:
            sys.exit('No file to read')
        if FILE.split('.')[-1].lower() != str(fileType):
            sys.exit("\nFile extension must be of type ." + str(fileType) + "\n")
        else:
           ### print 'Reading file: ' + str(FILE)
            return open(FILE, "r")


def lineFormat(line):
    line2 = []

    line2 = line.split(" ")
    line3 = line2[1:]
    line4 = []
    line4.append(line3[2])
    line4.append(line3[3])
    line4.append(line3[1])
    line4.append(three2one(line3[0]))
    return line4


def PDB_to_list(PDB_resi, chain_name, chain, MHCclass):
    outTxt = ''
    for i in range(0, len(PDB_resi)):
        if is_amino_acid(PDB_resi[i]):
            line = str(PDB_resi[i])
            line = line[8:]
            line = line.replace(' het=  resseq=', ' ')
            if line[-2].isalnum() == False:
                line = line.replace(' icode= >', ' ')
                line = line + chain_name + ' '
                if chain == 'MHCA' and MHCclass == 1:
                    number = int(re.sub("[^0-9]", "", line))
                    if 50 <= number <= 86:
                        line = line + 'MHCa1'
                    if 140 <= number <= 176:
                        line = line + 'MHCa2'
                if chain == 'MHCA' and MHCclass == 2:
                    number = int(re.sub("[^0-9]", "", line))
                    if 46 <= number <= 78:
                        line = line + 'MHCa1'
                if chain == 'MHCB' and MHCclass == 2:
                    number = int(re.sub("[^0-9]", "", line))
                    if 54 <= number <= 91:
                        line = line + 'MHCb1'
                lineF = lineFormat(line)

                for x in lineF:
                    outTxt += str(x) + "\t"
                outTxt += "\n"
            else:
               omitted = 1
    return outTxt


def pdb_to_sequence(pdb, chains, mhc_class, pdb_name):
    origPDB = read_file(pdb, "pdb")
    fileName = pdb.rsplit('.', 1)[0]

   ### print "PDB_name is ", pdb_name

    if not os.path.exists(pdb_name):
       ### print "Creating Directory " + pdb_name
        os.makedirs(fileName)

    if not os.path.exists(pdb_name + "/sequences"):
       ### print "Creating Directory " + pdb_name + "/sequences"
        os.makedirs(pdb_name + "/sequences")

    filtered_file = read_file(pdb, "pdb")

    # Sort chains
    mhca_a, mhc_b, pep, tcr_a, tcr_b = chains[0], chains[1], chains[2], chains[3], chains[4]

    # create a parser object
    parser = pdb_parser()

    # N.B. PDBParser is set to allow for buggy PDB structures by a defaut PERMISSIVE=1
    # can change this to =0 to throw an error outright. I have kept it to allow errors
    # as they often don't impact anything.

    struct = parser.get_structure(filtered_file, pdb)

   #### print "The warnings re: discontinous sequences are fine, ignore."

    mhca_chain = struct[0][mhca_a]
    mhcb_chain = struct[0][mhc_b]
    pep_chain = struct[0][pep]
    tcra_chain = struct[0][tcr_a]
    tcrb_chain = struct[0][tcr_b]

    # get residues
    mhca_residues = list(pdb_select.unfold_entities(mhca_chain, 'R'))
    mhcb_residues = list(pdb_select.unfold_entities(mhcb_chain, 'R'))
    pep_residues = list(pdb_select.unfold_entities(pep_chain, 'R'))
    tcra_residues = list(pdb_select.unfold_entities(tcra_chain, 'R'))
    tcrb_residues = list(pdb_select.unfold_entities(tcrb_chain, 'R'))

    out_txt = "Chain \tAnnotation \tResNum \tResCode\n"

    mhca_txt = PDB_to_list(mhca_residues, 'MHCA', 'MHCA', mhc_class)
    mchb_txt = PDB_to_list(mhcb_residues, 'MHCB', 'MHCB', mhc_class)
    peptide_txt = PDB_to_list(pep_residues, 'peptide', 'none', mhc_class)
    tcra_txt = PDB_to_list(tcra_residues, 'TCRA', 'alpha', mhc_class)
    tcrb_txt = PDB_to_list(tcrb_residues, 'TCRB', 'beta', mhc_class)

    out_txt += mhca_txt + mchb_txt + peptide_txt + tcra_txt + tcrb_txt

    outfile = open(pdb_name + '_sequence.txt', 'w')
    outfile.write(out_txt)
    outfile.close()


def find_term_from_id(data, chain, term):
    for entry in data:
        if chain in entry:
            for terms in entry:
                if term in terms:
                    return terms


def cdr_terms_extract_local(chain, CDRterm):
    name = (CDRterm.split("=")[0])
    alpha = '('
    bravo = ')'
    startpos = CDRterm.find(alpha) + len(alpha)
    endpos = CDRterm.find(bravo, startpos)
    locale = (CDRterm[startpos:endpos])
    locations = locale.split(",")
    positions = []
    positions.append(chain)
    positions.append(name)
    positions.append(locations)
    return positions


##################################################### BODY ##################################################
def add_cdr_to_sequence(fasta, seqTxt):
    fastaFile = read_file(fasta, "fasta")
    fastaFileName = fasta.rsplit('.', 1)[0]

    sequenceFile = read_file(seqTxt, "txt")
    sequenceFileName = seqTxt.rsplit('.', 1)[0]

    fastaEntries = fasta_parser(fastaFile)
    fastaEntries = depack_id(fastaEntries)

   #### print "############################################################################################"
    CDR1a = find_term_from_id(fastaEntries, "TCRA", "CDR1a")
    CDR1a = cdr_terms_extract_local("TCRA", CDR1a)
    CDR2a = find_term_from_id(fastaEntries, "TCRA", "CDR2a")
    CDR2a = cdr_terms_extract_local("TCRA", CDR2a)
    FWa = find_term_from_id(fastaEntries, "TCRA", "FWa")
    FWa = cdr_terms_extract_local("TCRA", FWa)

    CDR3a = find_term_from_id(fastaEntries, "TCRA", "CDR3a")
    CDR3a = cdr_terms_extract_local("TCRA", CDR3a)
    CDR1b = find_term_from_id(fastaEntries, "TCRB", "CDR1b")
    CDR1b = cdr_terms_extract_local("TCRB", CDR1b)
    CDR2b = find_term_from_id(fastaEntries, "TCRB", "CDR2b")
    CDR2b = cdr_terms_extract_local("TCRB", CDR2b)
    FWb = find_term_from_id(fastaEntries, "TCRB", "FWb")
    FWb = cdr_terms_extract_local("TCRB", FWb)

    CDR3b = find_term_from_id(fastaEntries, "TCRB", "CDR3b")
    CDR3b = cdr_terms_extract_local("TCRB", CDR3b)

    allLoops = [CDR1a, CDR2a, FWa, CDR3a, CDR1b, CDR2b, FWb, CDR3b]

    # Extracting sequence list #

    allSeqLines = everything_parser(sequenceFile)
    sequenceList = matrix_parser(allSeqLines)
    sequenceList = sequenceList[1:]

    # Pass fasta info into sequence list #

    for line in sequenceList:
        for loop in allLoops:
            if loop[0] == line[0]:
                #todo fix this mess below. clipping first element of string to accept first item of CDR loop
                if line[2] in [loop[2][0].split("[")[1]] + loop[2]:
                    line[1] = loop[1]

    # Make output #
    out = ""
    out += "Chain\tAnnotation\tResNum\tResCode\n"
    for line in sequenceList:
        for col in line:
            out += str(col) + "\t"
        out += "\n"

    # Write FASTA file #

   #### print "\nOutputting .txt file as: \n"
   #### print out

    txtOut = open(sequenceFileName + "_annot.txt", "w")

   #### print "Saving file as:" + sequenceFileName + "_annot.txt"

    txtOut.write(out)
    txtOut.close()

   #### print('     ~  End of annotateSequenceFile.py v0.2 BETA   ~')


def allLevels(contactMatrix, requestList):
    '''
    Request list must be a list containing an even number entries..
    odds are column you to be analysed
    evens are the value you are looking for in the previous odd entry
    '''
    counter = 0
    if len(requestList) <= 1:
       print "Request error"
    if len(requestList) % 2 != 0:
      print ("Request list error.. request list must be even")

    if len(requestList) == 2:
        for contactRow in contactMatrix:
            column = requestList[0]
            value = requestList[1]
            if str(contactRow[column]) == str(value):
                counter += 1

    if len(requestList) == 4:
        for contactRow in contactMatrix:
            if str(contactRow[requestList[0]]) == str(requestList[1]):
                if str(contactRow[requestList[2]]) == str(requestList[3]):
                    counter += 1

    if len(requestList) == 6:
        for contactRow in contactMatrix:
            if str(contactRow[requestList[0]]) == str(requestList[1]):
                if str(contactRow[requestList[2]]) == str(requestList[3]):
                    if str(contactRow[requestList[4]]) == str(requestList[5]):
                        counter += 1

    if len(requestList) == 8:
        for contactRow in contactMatrix:
            if str(contactRow[requestList[0]]) == str(requestList[1]):
                if str(contactRow[requestList[2]]) == str(requestList[3]):
                    if str(contactRow[requestList[4]]) == str(requestList[5]):
                        if str(contactRow[requestList[6]]) == str(requestList[7]):
                            counter += 1
    return counter


def stats(fileIn, inFileName):

    inFile = read_file(fileIn, "txt")

    if type(inFileName) != str:
        raise IOError('No file was loaded. Please view usage and provide a valid file to be processed')
    allLines = everything_parser(inFile)
    contactMatrix = matrix_parser(allLines)
    contactMatrix = contactMatrix[1:]

    output = ""
    outFile = open(inFileName + "_statistics.txt", 'w')
   #### print '\n Building statistics...\n'

    output += "\t\tpMHC\n"
    output += "TCRab\t\t" + str(len(contactMatrix)) + "\n"
    output += "TCRab\tvdW\t" + str(allLevels(contactMatrix, [15, "VW"])) + "\n"
    output += "TCRab\tHB\t" + str(allLevels(contactMatrix, [15, "HB"])) + "\n"
    output += "TCRab\tSB\t" + str(allLevels(contactMatrix, [15, "SB"])) + "\n"
    output += "--------------------" + "\n" + "\n"

    output += "\t\tpMHC" + "\n"

    output += "TCRab\tall\t" + str(len(contactMatrix)) + "\n"
    output += "\n"
    output += "TCRa\tall\t" + str(allLevels(contactMatrix, [1, "TCRa"])) + "\n"
    output += "TCRa\tCDR1a\t" + str(allLevels(contactMatrix, [1, "TCRa", 3, "CDR1a"])) + "\n"
    output += "TCRa\tCDR2a\t" + str(allLevels(contactMatrix, [1, "TCRa", 3, "CDR2a"])) + "\n"
    output += "TCRa\tCDR3a\t" + str(allLevels(contactMatrix, [1, "TCRa", 3, "CDR3a"])) + "\n"
    output += "TCRa\tFWa\t" + str(allLevels(contactMatrix, [1, "TCRa", 3, "FWa"])) + "\n"
    output += "TCRa\tNone\t" + str(allLevels(contactMatrix, [1, "TCRa", 3, ""])) + "\n"

    output += "\n"
    output += "TCRb\tall\t" + str(allLevels(contactMatrix, [1, "TCRb"])) + "\n"
    output += "TCRb\tCDR1b\t" + str(allLevels(contactMatrix, [1, "TCRb", 3, "CDR1b"])) + "\n"
    output += "TCRb\tCDR2b\t" + str(allLevels(contactMatrix, [1, "TCRb", 3, "CDR2b"])) + "\n"
    output += "TCRb\tCDR3b\t" + str(allLevels(contactMatrix, [1, "TCRb", 3, "CDR3b"])) + "\n"
    output += "TCRb\tFWb\t" + str(allLevels(contactMatrix, [1, "TCRb", 3, "FWb"])) + "\n"
    output += "TCRb\tNone\t" + str(allLevels(contactMatrix, [1, "TCRb", 3, ""])) + "\n"

    output += "--------------------" + "\n"
    output += "\t\tpMHCab\tMHCa\tMHCb\tpeptide" + "\n"

    output += "TCRa\tall\t" + str(allLevels(contactMatrix, [1, "TCRa"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 8, "MHCa"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 8, "MHCb"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 8, "peptide"])) + "\n"
    output += "TCRa\tCDR1a\t" + str(allLevels(contactMatrix, [1, "TCRa", 3, "CDR1a"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "CDR1a", 8, "MHCa"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "CDR1a", 8, "MHCb"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "CDR1a", 8, "peptide"])) + "\n"
    output += "TCRa\tCDR2a\t" + str(allLevels(contactMatrix, [1, "TCRa", 3, "CDR2a"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "CDR2a", 8, "MHCa"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "CDR2a", 8, "MHCb"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "CDR2a", 8, "peptide"])) + "\n"
    output += "TCRa\tCDR3a\t" + str(allLevels(contactMatrix, [1, "TCRa", 3, "CDR3a"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "CDR3a", 8, "MHCa"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "CDR3a", 8, "MHCb"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "CDR3a", 8, "peptide"])) + "\n"
    output += "TCRa\tFWa\t" + str(allLevels(contactMatrix, [1, "TCRa", 3, "FWa"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "FWa", 8, "MHCa"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "FWa", 8, "MHCb"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "Fwa", 8, "peptide"])) + "\n"
    output += "TCRa\tNone\t" + str(allLevels(contactMatrix, [1, "TCRa", 3, ""])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "", 8, "MHCa"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "", 8, "MHCb"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRa", 3, "", 8, "peptide"])) + "\n"
    output += "\n"
    output += "TCRb\tall\t" + str(allLevels(contactMatrix, [1, "TCRb"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 8, "MHCa"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 8, "MHCb"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 8, "peptide"])) + "\n"
    output += "TCRb\tCDR1b\t" + str(allLevels(contactMatrix, [1, "TCRb", 3, "CDR1b"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "CDR1b", 8, "MHCa"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "CDR1b", 8, "MHCb"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "CDR1b", 8, "peptide"])) + "\n"
    output += "TCRb\tCDR2b\t" + str(allLevels(contactMatrix, [1, "TCRb", 3, "CDR2b"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "CDR2b", 8, "MHCa"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "CDR2b", 8, "MHCb"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "CDR2b", 8, "peptide"])) + "\n"
    output += "TCRb\tCDR3b\t" + str(allLevels(contactMatrix, [1, "TCRb", 3, "CDR3b"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "CDR3b", 8, "MHCa"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "CDR3b", 8, "MHCb"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "CDR3b", 8, "peptide"])) + "\n"
    output += "TCRb\tFWb\t" + str(allLevels(contactMatrix, [1, "TCRb", 3, "FWb"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "FWb", 8, "MHCa"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "FWb", 8, "MHCb"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "Fwb", 8, "peptide"])) + "\n"
    output += "TCRb\tNone\t" + str(allLevels(contactMatrix, [1, "TCRb", 3, ""])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "", 8, "MHCa"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "", 8, "MHCb"])) + "\t" + str(
        allLevels(contactMatrix, [1, "TCRb", 3, "", 8, "peptide"])) + "\n"

    output += "--------------------"

    outFile.write(output)
    outFile.close()








