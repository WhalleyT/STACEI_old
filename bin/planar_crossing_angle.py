import os
import pymol
import sys
import time
import math

import numpy as np

import bin.data.colourSet as colour_set
import bin.data.viewSet as view_set


def read_file(infile, file_type):
    if infile is None:
        sys.exit('No file to read')
    if infile.split('.')[-1].lower() != str(file_type):
        sys.exit("\nFile extension must be of type ." + str(file_type) + "\n")
    else:
        print 'Reading file: ' + str(infile)
        return open(infile, "r")


def file_to_list(infile):
    all_lines = []
    for line in infile:
        all_lines.append(line)
    return all_lines


def initialise_pymol():
    print "\nInitialising pymol...\n"
    pymol.finish_launching(['pymol', '-qeim'])
    pymol.cmd.reinitialize()
    # set PyMOL parameters
    pymol.cmd.set("ray_shadows", "0")
    pymol.cmd.set("specular", "off")
    pymol.cmd.set("orthoscopic", "on")
    pymol.cmd.bg_color("white")
    pymol.cmd.set("ray_opaque_background", "0")
    return None


def align2template(fileName, MHCClass):
    print os.getcwd()
    pymol.cmd.delete("all")
    print "\nAligning file to template...\n"
    pymol.cmd.load(fileName + ".pdb")

    pymol.cmd.load("bin/data/" + MHCClass + "_template.pdb")

    pymol.cmd.align(fileName, MHCClass + "_template")
    pymol.cmd.save(fileName + "/crossingAngle/" + fileName + "_aligned_noMeta.pdb", fileName)

    print "\nAlignment to " + MHCClass + "_template.pdb complete!\n"
    return None


### SSBOND Metadata adder ###

def SSParser(inFile):
    print "Finding SSBOND Lines..."
    outSS = []
    for line in inFile:
        add = ''
        if "SSBOND" in line[0:6]:
            add += line
            outSS.append(add)
            add = ''
    return outSS


def headerParser(inFile):
    print "Finding header..."
    outHeader = []
    for line in inFile:
        add = ''
        if "HEADER" in line[0:6]:
            add += line
            outHeader.append(add)
            add = ''
    return outHeader


def titleParser(inFile):
    print "Finding title..."
    outTitle = []
    for line in inFile:
        add = ''
        if "TITLE" in line[0:6]:
            add += line
            outTitle.append(add)
            add = ''
    return outTitle


### Ray tools ###

def wait4ray(query):
    counter = 0
    while not os.path.exists(query):
        print ("=" * counter) + "| " + str(counter)
        time.sleep(1)
        counter += 1
    return None


def rayTime(saveas):
    print "Outputting image.. This may take a few seconds.."
    if os.path.exists(saveas):
        print "Removing " + saveas + " as it already exists!"
        os.remove(saveas)
    time.sleep(10)
    pymol.cmd.png(saveas, ray=1, width=3000, height=3000, dpi=300)
    wait4ray(saveas)
    print "Done! " + str(saveas) + " was outputted"


### PDB Atom Parser ###

def atomParser(allLines):
    print str(len(allLines)) + " went in"
    atomList = []
    elseList = []
    for line in allLines:
        if "ATOM" in line[0:5]:
            atomList.append(line)
        else:
            elseList.append(line)
    print str(len(atomList)) + " atom lines detected"
    print str(len(elseList)) + " else lines detected"
    return atomList, elseList


def atomRowParser(atomRow):
    atomRowParsed = []
    atomRowParsed.append(atomRow[0:6])  # 0 record name
    atomRowParsed.append(int(atomRow[6:11]))  # 1 atom serial numeber
    atomRowParsed.append(atomRow[12:16])  # 2 atom name
    atomRowParsed.append(atomRow[16])  # 3 alternate location indicator
    atomRowParsed.append(atomRow[17:20])  # 4 residue name
    atomRowParsed.append(atomRow[21])  # 5 chain identifier
    atomRowParsed.append(int(atomRow[22:26]))  # 6 residue sequence number
    atomRowParsed.append(atomRow[26])  # 7 code for insertion of residues
    atomRowParsed.append(float(atomRow[30:38]))  # 8 coordinate x
    atomRowParsed.append(float(atomRow[38:46]))  # 9 coordinate y
    atomRowParsed.append(float(atomRow[46:54]))  # 10 coordinate z
    atomRowParsed.append(float(atomRow[54:60]))  # 11 occupancy
    atomRowParsed.append(float(atomRow[60:66]))  # 12 temperature b factor
    atomRowParsed.append(atomRow[76:78])  # 13 element symbol
    atomRowParsed.append(atomRow[78:-1])  # 14 charge on atom
    return atomRowParsed


def make_atom_matrix(atomList):
    atomMatrix = []
    print '\nCreating atom matrix...'
    for row in atomList:
        atomMatrix.append(atomRowParser(row))
    print "Done!\n"
    return atomMatrix


def TCRpairFinder(SSBondList, TCRachain, TCRbchain):
    '''
    Takes in all the SSBonds and outputs two lists..

    the alpha chain cys and the beta chain cys
    '''
    print "Finding TCRa and TCRb cysteine pairs...\n "
    alphaCys = []
    betaCys = []
    for row in SSBondList:
        if row[0] == TCRachain:
            if 15 <= row[1] <= 35:
                if row[2] == TCRachain:
                    if 85 <= row[3] <= 110:
                        alphaCys = row
        if row[0] == TCRachain:
            if 15 <= row[3] <= 35:
                if row[2] == TCRachain:
                    if 85 <= row[1] <= 110:
                        alphaCys = row
    for row in SSBondList:
        if row[0] == TCRbchain:
            if 15 <= row[1] <= 35:
                if row[2] == TCRbchain:
                    if 85 <= row[3] <= 110:
                        betaCys = row
        if row[0] == TCRbchain:
            if 15 <= row[3] <= 35:
                if row[2] == TCRbchain:
                    if 85 <= row[1] <= 110:
                        betaCys = row
    print alphaCys
    print betaCys
    print "\n"
    if len(alphaCys) > 0 and len(betaCys) > 0:
        return alphaCys, betaCys
    else:
        return False, False


def TCRpairAtomFinder(TCRpairLocation, atomMatrix):
    print "Finding SG atoms..."

    print TCRpairLocation
    coordinate1 = []
    coordinate2 = []
    chain = ''
    cysResidue1 = 0
    cysResidue2 = 0
    chain = TCRpairLocation[0]
    cysResidue1 = TCRpairLocation[1]
    cysResidue2 = TCRpairLocation[3]
    for atoms in atomMatrix:
        if atoms[5] == chain:
            if atoms[6] == cysResidue1:
                if "SG" in atoms[2]:
                    if atoms[3] == '' or 'A':
                        coordinate1 = atoms
    for atoms in atomMatrix:
        if atoms[5] == chain:
            if atoms[6] == cysResidue2:
                if "SG" in atoms[2]:
                    if atoms[3] == '' or 'A':
                        coordinate2 = atoms
    return coordinate1, coordinate2


def fastaParser(fasta):
    from Bio import SeqIO
    entries = []
    for seq_record in SeqIO.parse(fasta, "fasta"):
        ID = ''
        sequence = ''
        entry = []
        ID = "".join(seq_record.id)
        sequence = "".join(seq_record.seq)
        entry.append(ID)
        entry.append(sequence)
        entries.append(entry)
    return entries


def depackID(entries):
    newEntries = []
    for entry in entries:
        ID = entry[0]
        sequence = entry[1]
        newEntry = []
        id_terms = ID.split("|")
        # print ID.split("|")
        for term in id_terms:
            if len(term) != 0:
                newEntry.append(term)
        newEntry.append(sequence)

        newEntries.append(newEntry)

    return newEntries


def findLocations(subentries):
    s = "=[]"
    locations = []
    for col in subentries:
        if s[0] in col and s[1] in col and s[2] in col:
            locations.append(col)
    return locations


def depackLocations(subentries):
    locations = []
    for col in subentries:
        location = []
        location.append(col.rsplit("=", 1)[0])
        locationstring = (col.partition('[')[-1].rpartition(']')[0])
        location += map(int, locationstring.split(','))
        locations.append(location)
    return locations


def findCysLocs(locations):
    for locs in locations:
        if "Cys1" in locs[0]:
            cys1 = locs
        if "Cys2" in locs[0]:
            cys2 = locs
    return cys1, cys2


### Cysteine pair finder from SSBOND list ###

def SSBondParser(inFile):
    print "\nFinding disulphide bridges..."
    inFile
    SSBondList = []
    for line in inFile:
        if "SSBOND" in line[0:6]:
            parsedLine = []
            parsedLine.append(line[15])  # Chain 1
            parsedLine.append(int(line[18:21]))  # Residue 1
            parsedLine.append(line[29])  # Chain 2
            parsedLine.append(int(line[32:35]))  # Residue 1
            SSBondList.append(parsedLine)
    print str(len(SSBondList)) + " SSBOND lines detected!\n"
    if len(SSBondList) > 0:
        return SSBondList
    else:
        return False


def SGcoords(row):
    x = ''
    y = ''
    z = ''
    SGcoords = []
    x = row[8]
    y = row[9]
    z = row[10]
    SGcoords = [x, y, z]
    return SGcoords


def cysFromSSBONDWrapper(inFile, atomMatrix, TCRachain, TCRbchain):
    SSLines = SSBondParser(inFile)
    if not SSLines:
        print "No SSBOND lines were detected in the PDB!\n"
        return False, False, False, False
    a, b = TCRpairFinder(SSLines, TCRachain, TCRbchain)
    if a == False or b == False:
        print "Could not find TCR cysteine pairs in the SSBOND list!\n"
        return False, False, False, False
    A1atoms, A2atoms = TCRpairAtomFinder(a, atomMatrix)
    B1atoms, B2atoms = TCRpairAtomFinder(b, atomMatrix)
    if A1atoms == False or A2atoms == False or B1atoms == False or B2atoms == False:
        print "Could not find TCR cysteine atoms from the cysteine pairs listed in SSBOND list!\n"
        return False, False, False, False
    return A1atoms, A2atoms, B1atoms, B2atoms


### MHC coordinates finder ###

def MHCaxisI(atomMatrix, MHCachain):
    print "\nMHC axis atoms are:\n"
    Calphas = []
    helixCalphas = []
    for row in atomMatrix:
        if row[5] == MHCachain and "CA" in row[2]:
            Calphas.append(row)
    for row in Calphas:
        if 50 <= row[6] <= 86 or 140 <= row[6] <= 176:
            if row[3] == '' or 'A':
                helixCalphas.append(row)
    for x in helixCalphas:
        print x
    print '\nNumber of Calpha atoms in MHC axis: ' + str(len(helixCalphas)) + "\n"
    return helixCalphas


def MHCaxisII(atomMatrix, MHCachain, MHCbchain):
    print "\nMHC axis atoms are:\n"
    Calphas = []
    helixCalphas = []
    for row in atomMatrix:
        if row[5] == MHCachain and "CA" in row[2]:
            Calphas.append(row)
    for row in Calphas:
        if 46 <= row[6] <= 78:
            if row[3] == '' or 'A':
                helixCalphas.append(row)
    Calphas = []
    for row in atomMatrix:
        if row[5] == MHCbchain and "CA" in row[2]:
            Calphas.append(row)
    for row in Calphas:
        if 54 <= row[6] <= 64:
            if row[3] == '' or 'A':
                helixCalphas.append(row)
        if 67 <= row[6] <= 91:
            if row[3] == '' or 'A':
                helixCalphas.append(row)
    for x in helixCalphas:
        print x
    print '\nNumber of Calpha atoms in MHC axis: ' + str(len(helixCalphas)) + "\n"
    return helixCalphas


def mhc_axis_coords(helix_C_alphas):
    print "Extracting MHC axis coordinates..."
    mhc_x = []
    mhc_y = []
    mhc_z = []
    for row in helix_C_alphas:
        mhc_x.append(row[8])
        mhc_y.append(row[9])
        mhc_z.append(row[10])
    mhc_coords = [mhc_x, mhc_y, mhc_z]
    return mhc_coords


def find_mhc_axis_coords(atom_matrix, mhc_class, a_chain, b_chain):
    if mhc_class == "I":
        print 'Establishing MHC I axis...'
        return mhc_axis_coords(MHCaxisI(atom_matrix, a_chain))
    if mhc_class == "II":
        print 'Establishing MHC II axis...'
        return mhc_axis_coords(MHCaxisII(atom_matrix, a_chain, b_chain))


def calc_plane_bis(x, y, z):
    a = np.column_stack((x, y, z))
    return np.linalg.lstsq(a, np.ones_like(x))[0]


def project_points(x, y, z, a, b, c):
    """
    Projects the points with coordinates x, y, z onto the plane
    defined by a*x + b*y + c*z = 1
    """
    vector_norm = a * a + b * b + c * c
    normal_vector = np.array([a, b, c]) / np.sqrt(vector_norm)
    point_in_plane = np.array([a, b, c]) / vector_norm

    points = np.column_stack((x, y, z))
    points_from_point_in_plane = points - point_in_plane
    proj_onto_normal_vector = np.dot(points_from_point_in_plane,
                                     normal_vector)
    proj_onto_plane = (points_from_point_in_plane -
                       proj_onto_normal_vector[:, None] * normal_vector)

    return point_in_plane + proj_onto_plane


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


# Quadrant finder #

def whichQuadrant(TCR):
    TCRAox = np.asscalar(TCR[0][0])
    TCRBox = np.asscalar(TCR[1][0])
    TCRAoy = np.asscalar(TCR[0][1])
    TCRBoy = np.asscalar(TCR[1][1])

    if TCRBox < TCRAox and TCRBoy < TCRAoy:
        print "Crossing angle is in the 0-90 quadrant\n"
        return 0, 90

    if TCRBox > TCRAox and TCRBoy < TCRAoy:
        print "Crossing angle is in the 90-180 quadrant\n"
        return 90, 180

    if TCRBox > TCRAox and TCRBoy > TCRAoy:
        print "Crossing angle is in the 180-270 quadrant\n"
        return 180, 270

    if TCRBox < TCRAox and TCRBoy > TCRAoy:
        print "Crossing angle is in the 180-270 quadrant\n"
        return 270, 360
    else:
        print "oh boy"
        return None


def whichCrossingAngle(crossingAngle, rangeLow):
    if rangeLow == 0 or rangeLow == 90:
        print "TCR binds with canonical polarity\n"
        return crossingAngle
    if rangeLow == 180 or rangeLow == 270:
        return 360 - crossingAngle
        print "TCR binds with reverse polarity\n"


########### BODY    ###############

### Initialiser ###

# Load input.pdb #
def calculate_and_print(pdb, fasta, mhc_class, ray, chains):

    if mhc_class is 1:
        mhc_class = "I"
    elif mhc_class is 2:
        mhc_class = "II"
    else:
        raise Exception("MHC class invalid")


    origPDB = read_file(pdb, "pdb")
    fileName = pdb.rsplit('.', 1)[0]

    if type(fileName) != str:
        sys.exit('No file to grab metadata was loaded. Please view usage and provide a valid .pdb to be processed')

    # Sort chains
    MHCachain, MHCbchain, peptidechain, TCRachain, TCRbchain = chains[0], chains[1], chains[2], chains[3], chains[4]

    # Make output folder #

    if not os.path.exists(fileName):
        print "Creating Directory " + fileName
        os.makedirs(fileName)

    if not os.path.exists(fileName + "/crossingAngle"):
        print "Creating Directory " + fileName + "/crossingAngle"
        os.makedirs(fileName + "/crossingAngle")

    # Align input.pdb to template #
    initialise_pymol()
    align2template(fileName, mhc_class)
    print "Opening aligned PDB file"
    alignedPDBnoMeta = read_file(fileName + "/crossingAngle/" + fileName + "_aligned_noMeta.pdb", "pdb")

    # Extract all lines from PDB file of orig and aligned PDB files #

    origallLines = file_to_list(origPDB)
    aligned_all_lines = file_to_list(alignedPDBnoMeta)
    alignedPDBnoMeta.close()

    ### Metadata Replacer ###

    print "Transfering metadata from " + fileName + ".pdb to " + fileName + "_aligned_noMeta\n"
    metatitle = titleParser(origallLines)
    metaheader = headerParser(origallLines)
    meta_ss = SSParser(origallLines)
    working_file = open(fileName + "/crossingAngle/" + fileName + '_aligned.pdb', 'w')
    outTxtMeta = ''
    outTxtMeta += "         " + str(1) + "         " + str(2) + "         " + str(3) + "         " + str(4) + \
                  "         " + str(5) + "         " + str(6) + "         " + str(7) + "         " + str(8) + "\n"
    outTxtMeta += str(12345678901234567890123456789012345678901234567890123456789012345678901234567890) + "\n"
    outTxtMeta += "REMARK   1                                                                      \n\
    REMARK   1  ~~~~ Metadata generated by crossingAngle_v2.0.py     ~~~~ \n\
    REMARK   1                                                                      \n\
    REMARK   1  This metadata contains :           \n\
    REMARK   1      HEADER                        \n\
    REMARK   1      TITLE         \n\
    REMARK   1      This REMARK \n\
    REMARK   1      SSBOND lines required for crossing angle calculations    \n\
    REMARK   1\n"
    
    for meta in metaheader:
        outTxtMeta += meta
        
    for meta in metatitle:
        outTxtMeta += meta
        
    for meta in meta_ss:
        outTxtMeta += meta
        
    for meta in aligned_all_lines:
        outTxtMeta += meta
        
    working_file.write(outTxtMeta)
    
    print "Metadata added to " + fileName + "/crossingAngle/" + fileName + "_aligned.pdb\n"

    # Load the working PDB file #
    print "Extracting all lines from " + fileName + "/crossingAngle/" + fileName + "_aligned.pdb\n"
    working_file = open(fileName + "/crossingAngle/" + fileName + '_aligned.pdb', 'r')
    
    all_lines = file_to_list(working_file)
    
    working_file.close()
    
    atom_list, else_list = atomParser(all_lines)
    atomMatrix = make_atom_matrix(atom_list)

    ############################# Cysteine axis #########################################

    if fasta is not None:
        print "Perfect! an annotated fasta file was provided to flag where the cysteine pair residues are!"
        # Unpack the fasta file with Cys information #
        fastaFile = read_file(fasta, "fasta")
        fastaEntries = fastaParser(fasta)
        fastaEntries = depackID(fastaEntries)

        TCRAflag = []
        for entry in fastaEntries:
            if "TCRA" in entry:
                TCRAflag = entry
        TCRAflaglocations = findLocations(TCRAflag)
        TCRAflaglocations = depackLocations(TCRAflaglocations)

        aCys1, aCys2 = findCysLocs(TCRAflaglocations)
        print aCys1
        TCRaCys = [TCRachain, aCys1[1], TCRachain, aCys2[1]]

        TCRBflag = []
        for entry in fastaEntries:
            if "TCRB" in entry:
                TCRBflag = entry
        TCRBflaglocations = findLocations(TCRBflag)
        TCRBflaglocations = depackLocations(TCRBflaglocations)
        TCRbCys1, TCRbCys2 = findCysLocs(TCRBflaglocations)

        bCys1, bCys2 = findCysLocs(TCRBflaglocations)
        TCRbCys = [TCRbchain, bCys1[1], TCRbchain, bCys2[1]]

        TCRaCys1, TCRaCys2 = TCRpairAtomFinder(TCRaCys, atomMatrix)
        TCRbCys1, TCRbCys2 = TCRpairAtomFinder(TCRbCys, atomMatrix)

    else:
        print "Determining the TCR axis via locating the cys pair atom coordinates...\n"
        TCRaCys1, TCRaCys2, TCRbCys1, TCRbCys2 = cysFromSSBONDWrapper(else_list, atomMatrix, TCRachain, TCRbchain)

    cysCheck = [TCRaCys1, TCRaCys2, TCRbCys1, TCRbCys2]
    if any(x == False for x in cysCheck):
        pymol.cmd.quit()
        sys.exit("Could not find the TCR disulphide bridge in the TCR residues. \
                 Please ensure the PDB file contains SSBOND list in header or provide an annotated fasta file generated by complexSequenceTools.py\
                 See Usage for more details.")

    print "\nTCRa and TCRb cysteine pair SG atoms are:\n"
    print TCRaCys1
    print TCRaCys2
    print TCRbCys1
    print TCRbCys2
    cysAlpha1 = SGcoords(TCRaCys1)
    cysAlpha2 = SGcoords(TCRaCys2)
    cysBeta1 = SGcoords(TCRbCys1)
    cysBeta2 = SGcoords(TCRbCys2)
    print "\nTCRa and TCRb cysteine pair SG coordinates are [x,y,z]:\n"
    print cysAlpha1
    print cysAlpha2
    print cysBeta1
    print cysBeta1
    TCRA = []
    TCRA.append((cysAlpha1[0] + cysAlpha2[0]) / 2)
    TCRA.append((cysAlpha1[1] + cysAlpha2[1]) / 2)
    TCRA.append((cysAlpha1[2] + cysAlpha2[2]) / 2)

    TCRA = np.array(TCRA)

    TCRB = []
    TCRB.append((cysBeta1[0] + cysBeta2[0]) / 2)
    TCRB.append((cysBeta1[1] + cysBeta2[1]) / 2)
    TCRB.append((cysBeta1[2] + cysBeta2[2]) / 2)

    TCRB = np.array(TCRB)

    print "\nTCRa cysteine centroid is [x,y,z]:"
    print TCRA
    print "\nTCRb cysteine centroid is [x,y,z]:"
    print TCRB
    print "\n\n"

    TCR = TCR = np.vstack((TCRA, TCRB))

    print "\nTCR is..."
    print type(TCR), "shape = ", TCR.shape
    print TCR, "\n"

    ################################## MHC axis #################################

    MHCcoords = find_mhc_axis_coords(atomMatrix, mhc_class, MHCachain, MHCbchain)
    MHCdata = np.array(MHCcoords)
    MHCdata = MHCdata.T
    print "Data is stored as a ", MHCdata.shape, "numpy array."

    # Calculate the mean of the points, i.e. the 'center' of the cloud
    print "Finding mean of MHC datapoints..."
    MHCdata2 = MHCdata

    MHCdata2mean = MHCdata2.mean(axis=0)

    ############## Got the data.. now the fun begins ########################
    print "\nGot the TCR and MHC coordinate data. Beginnning line and plane fitting..\n"

    ######### LOBF through MHC data ##################
    # Do an SVD on the mean-centered MHCdata.
    print "Performing SVD line of best fit (LOBF) on mean-centered MHC datapoints..."
    print "Generating a line of best fit through the MHC datapoints..."
    print "NB This uses the np.linalg.svd algorithm which is a more stable alternative to the largest eigenvalue method used previous."
    uu, dd, vv = np.linalg.svd(MHCdata2 - MHCdata2mean)

    # Now generate some points along this best fit line, for plotting.
    # A scale factor (sf) is used to determine the arbitary length of this line
    sf = 10
    # Straight line, so we only need 2 points.
    print "Generating two points along the best fit line. Using a scale factor of: ", str(sf)
    MHC = vv[0] * np.mgrid[-sf:sf:2j][:, np.newaxis]

    # shift by the mean to get the line in the right place
    print "Shifting best fit line to the mean.. "
    MHC += MHCdata2mean
    MHC = np.flipud(MHC)

    print "\nMHC is..."
    print type(MHC), "shape = ", MHC.shape
    print MHC

    print "Generating MHC LOBF points transformed to TCRB for better visualisation"

    mhc_to_tcrb = MHC[0] - TCR[1]
    mhc_at_tcrb = MHC - mhc_to_tcrb

    print "\nMHCatTCRB is..."
    print type(mhc_at_tcrb), "shape = ", mhc_to_tcrb.shape
    print mhc_at_tcrb

    # Plane of best fit through MHC data

    # Create plane data points in X,Y
    print "\nCreating a plane of best fit through the MHC to which we can project MHC LOBF and TCR centroids onto..\n"
    mn = np.min(MHCdata, axis=0)  # used to set range of the plane
    mx = np.max(MHCdata, axis=0)  # used to set range of the plane
    X, Y = np.meshgrid(np.linspace(mn[0] - 20, mx[0] + 20, 20),
                       np.linspace(mn[1] - 20, mx[1] + 20, 20))  # min/max +/- 20 seems to capture a good area
    print "Calculating plane.."
    print MHCdata[:, 0]
    print MHCdata[:, 1]
    print MHCdata[:,2]
    projectionPlane = calc_plane_bis(MHCdata[:, 0], MHCdata[:, 1], MHCdata[:,2])
    # this creates the coefficients of the plane in form ax + by +cz = d where d = 1
    # projectionPlane is a numpy array of shape (3,) i.e. 3 x 1 storing coeffcients a,b,c

    print "The coefficients of the MHC plane are.."
    print type(projectionPlane), "shape = ", projectionPlane.shape
    print projectionPlane, "\n"

    # evaluate it on grid
    print "Moving random generated plane coordinates to the plane described by the plane coefficients..\n"
    Z = (1 - projectionPlane[0] * X - projectionPlane[1] * Y) / projectionPlane[
        2]  # rearrangement of above plane equation for values of Z
    print "..Done!\n"

    # Project TRA and TRB points to plane
    print "Now we have the MHC plane.. performing an orthogonal projection of TCR and MHC LOBF points to this plane"
    print "Projecting TCR points to the generated MHC plane..\n"
    TRAproj = project_points(TCR[0, 0], TCR[0, 1], TCR[0, 2], projectionPlane[0], projectionPlane[1], projectionPlane[2])
    TRBproj = project_points(TCR[1, 0], TCR[1, 1], TCR[1, 2], projectionPlane[0], projectionPlane[1], projectionPlane[2])
    TCRproj = np.concatenate((TRAproj, TRBproj), axis=0)

    print "TCRproj is..."
    print type(TCRproj), "shape = ", TCRproj.shape
    print TCRproj, "\n"

    print "..Done!\n"

    ######### Projection of points to the plane ##############
    # Project MHC points to plane
    print "Projecting MHC line on to the generated MHC plane..\n"
    MHCAproj = project_points(MHC[0, 0], MHC[0, 1], MHC[0, 2], projectionPlane[0], projectionPlane[1], projectionPlane[2])
    MHCBproj = project_points(MHC[1, 0], MHC[1, 1], MHC[1, 2], projectionPlane[0], projectionPlane[1], projectionPlane[2])
    MHCproj = np.concatenate((MHCAproj, MHCBproj), axis=0)

    print "\nMHCproj is..."
    print type(MHCproj), "shape = ", MHCproj.shape
    print MHCproj, "\n"
    print "..Done!\n"

    print "Projecting MHCatTCRB line to the generated MHC plane..\n"
    MHCAproj2 = project_points(mhc_at_tcrb[0, 0], mhc_at_tcrb[0, 1], mhc_at_tcrb[0, 2], projectionPlane[0], projectionPlane[1],
                               projectionPlane[2])
    MHCBproj2 = project_points(mhc_at_tcrb[1, 0], mhc_at_tcrb[1, 1], mhc_at_tcrb[1, 2], projectionPlane[0], projectionPlane[1],
                               projectionPlane[2])
    MHCatTCRBP = np.concatenate((MHCAproj2, MHCBproj2), axis=0)

    print "\nMHCatTCRBProj is..."
    print type(MHCatTCRBP), "shape = ", MHCatTCRBP.shape
    print MHCatTCRBP, "\n"
    print "..Done!\n"

    # Vector between TCR and TCRproj

    print "Finally, the plane is transformed to TCRB to create an axis parallel to the MHC but at the TCRB location."
    print "This is performed by transforming TCRproj as only a line (not a plane) is required for this calculation.\n"

    TCR2TCRprojV = TCRproj - TCR
    print "Vector between TCR and TCRproj is.."
    print TCR2TCRprojV

    planeAtTCR = TCRproj - TCR2TCRprojV[1, :]

    print "\nRepresentation of the plane at the TCR is.."
    print planeAtTCR

    ####################### PDB file writer of all established points ################################

    print "We've established all the required coordinates.. let's save them first.."

    corner1txt = [X[0][0], Y[0][0], Z[0][0]]
    corner2txt = [X[0][-1], Y[0][-1], Z[0][-1]]
    corner3txt = [X[-1][0], Y[-1][0], Z[-1][0]]
    corner4txt = [X[-1][-1], Y[-1][-1], Z[-1][-1]]
    MHCatxt = MHC[0].tolist()
    MHCbtxt = MHC[1].tolist()
    TCRatxt = TCR[0].tolist()
    TCRbtxt = TCR[1].tolist()

    TCRaPtxt = TRAproj[0].tolist()
    TCRbPtxt = TRBproj[0].tolist()

    MHCaatTCRBtxt = mhc_at_tcrb[0].tolist()
    MHCbatTCRBtxt = mhc_at_tcrb[1].tolist()
    MHCaPtxt = MHCAproj[0].tolist()
    MHCbPtxt = MHCBproj[0].tolist()
    MHCaatTCRBPtxt = MHCatTCRBP[0].tolist()
    MHCbatTCRBPtxt = MHCatTCRBP[1].tolist()

    planeAtTCRatxt = planeAtTCR[0].tolist()
    planeAtTCRbtxt = planeAtTCR[1].tolist()

    outPDB = ''
    outPDBfile = open(fileName + "/crossingAngle/" + fileName + "_planeCorners.pdb", mode="w")

    # File header

    outPDB += "\
    REMARK   1                                                                \n\
    REMARK   1  ~~~~         Generated by crossingAngle v1.9.py          ~~~~ \n\
    REMARK   1                                                                \n\
    REMARK   1  This PDB File contains 18 arbitrary atoms that the            \n\
    REMARK   1  crossing angle calculates from:                               \n\
    REMARK   1      MHC plane corner 1                                        \n\
    REMARK   1      MHC plane corner 2                                        \n\
    REMARK   1      MHC plane corner 3                                        \n\
    REMARK   1      MHC plane corner 4                                        \n\
    REMARK   1      MHCa fit                                                  \n\
    REMARK   1      MHCb fit                                                  \n\
    REMARK   1      MHCa fit at TCRB                                          \n\
    REMARK   1      MHCb fit at TCRB                                          \n\
    REMARK   1      MHCa fit projected onto plane                             \n\
    REMARK   1      MHCb fit projected onto plane                             \n\
    REMARK   1      MHCa fit at TCRB projected onto plane                     \n\
    REMARK   1      MHCb fit at TCRB projected onto plane                     \n\
    REMARK   1      TCRa centroid                                             \n\
    REMARK   1      TCRb centroid                                             \n\
    REMARK   1      TCRa centroid projected onto plane                        \n\
    REMARK   1      TCRb centroid projected onto plane                        \n\
    REMARK   1      Plane at TCRa                                             \n\
    REMARK   1      Plane at TCRb                                             \n\
    REMARK   1                                                                \n"

    # Stuff to write to file

    a1, a2, a3 = str("%.3f" % corner1txt[0]), str("%.3f" % corner1txt[1]), str("%.3f" % corner1txt[2])
    a1, a2, a3 = a1.rjust(8), a2.rjust(8), a3.rjust(8)
    outPDB += "ATOM      1  SG  CYS P   1    " + a1 + a2 + a3 + "  1.00  1.00           S" + "\n"

    b1, b2, b3 = str("%.3f" % corner2txt[0]), str("%.3f" % corner2txt[1]), str("%.3f" % corner2txt[2])
    b1, b2, b3 = b1.rjust(8), b2.rjust(8), b3.rjust(8)
    outPDB += "ATOM      2  SG  CYS P   2    " + b1 + b2 + b3 + "  1.00  1.00           S" + "\n"

    c1, c2, c3 = str("%.3f" % corner3txt[0]), str("%.3f" % corner3txt[1]), str("%.3f" % corner3txt[2])
    c1, c2, c3 = c1.rjust(8), c2.rjust(8), c3.rjust(8)
    outPDB += "ATOM      3  SG  CYS P   3    " + c1 + c2 + c3 + "  1.00  1.00           S" + "\n"

    d1, d2, d3 = str("%.3f" % corner4txt[0]), str("%.3f" % corner4txt[1]), str("%.3f" % corner4txt[2])
    d1, d2, d3 = d1.rjust(8), d2.rjust(8), d3.rjust(8)
    outPDB += "ATOM      4  SG  CYS P   4    " + d1 + d2 + d3 + "  1.00  1.00           S" + "\n"

    e1, e2, e3 = str("%.3f" % MHCatxt[0]), str("%.3f" % MHCatxt[1]), str("%.3f" % MHCatxt[2])
    e1, e2, e3 = e1.rjust(8), e2.rjust(8), e3.rjust(8)
    outPDB += "ATOM      5  SG  CYS M   1    " + e1 + e2 + e3 + "  1.00  1.00           S" + "\n"

    f1, f2, f3 = str("%.3f" % MHCbtxt[0]), str("%.3f" % MHCbtxt[1]), str("%.3f" % MHCbtxt[2])
    f1, f2, f3 = f1.rjust(8), f2.rjust(8), f3.rjust(8)
    outPDB += "ATOM      6  SG  CYS M   2    " + f1 + f2 + f3 + "  1.00  1.00           S" + "\n"

    i1, i2, i3 = str("%.3f" % MHCaPtxt[0]), str("%.3f" % MHCaPtxt[1]), str("%.3f" % MHCaPtxt[2])
    i1, i2, i3 = i1.rjust(8), i2.rjust(8), i3.rjust(8)
    outPDB += "ATOM      7  SG  CYS N   1    " + i1 + i2 + i3 + "  1.00  1.00           S" + "\n"

    j1, j2, j3 = str("%.3f" % MHCbPtxt[0]), str("%.3f" % MHCbPtxt[1]), str("%.3f" % MHCbPtxt[2])
    j1, j2, j3 = j1.rjust(8), j2.rjust(8), j3.rjust(8)
    outPDB += "ATOM      8  SG  CYS N   2    " + j1 + j2 + j3 + "  1.00  1.00           S" + "\n"

    o1, o2, o3 = str("%.3f" % MHCaatTCRBtxt[0]), str("%.3f" % MHCaatTCRBtxt[1]), str("%.3f" % MHCaatTCRBtxt[2])
    o1, o2, o3 = o1.rjust(8), o2.rjust(8), o3.rjust(8)
    outPDB += "ATOM      9  CA  LEU O   1    " + o1 + o2 + o3 + "  1.00  1.00           S" + "\n"

    p1, p2, p3 = str("%.3f" % MHCbatTCRBtxt[0]), str("%.3f" % MHCbatTCRBtxt[1]), str("%.3f" % MHCbatTCRBtxt[2])
    p1, p2, p3 = p1.rjust(8), p2.rjust(8), p3.rjust(8)
    outPDB += "ATOM     10  N   ALA O   2    " + p1 + p2 + p3 + "  1.00  1.00           S" + "\n"

    q1, q2, q3 = str("%.3f" % MHCaatTCRBPtxt[0]), str("%.3f" % MHCaatTCRBPtxt[1]), str("%.3f" % MHCaatTCRBPtxt[2])
    q1, q2, q3 = q1.rjust(8), q2.rjust(8), q3.rjust(8)
    outPDB += "ATOM     11  CA  LEU Q   1    " + q1 + q2 + q3 + "  1.00  1.00           S" + "\n"

    r1, r2, r3 = str("%.3f" % MHCbatTCRBPtxt[0]), str("%.3f" % MHCbatTCRBPtxt[1]), str("%.3f" % MHCbatTCRBPtxt[2])
    r1, r2, r3 = r1.rjust(8), r2.rjust(8), r3.rjust(8)
    outPDB += "ATOM     12  N   ALA Q   2    " + r1 + r2 + r3 + "  1.00  1.00           S" + "\n"

    g1, g2, g3 = str("%.3f" % TCRatxt[0]), str("%.3f" % TCRatxt[1]), str("%.3f" % TCRatxt[2])
    g1, g2, g3 = g1.rjust(8), g2.rjust(8), g3.rjust(8)
    outPDB += "ATOM     13  SG  CYS T   1    " + g1 + g2 + g3 + "  1.00  1.00           S" + "\n"

    h1, h2, h3 = str("%.3f" % TCRbtxt[0]), str("%.3f" % TCRbtxt[1]), str("%.3f" % TCRbtxt[2])
    h1, h2, h3 = h1.rjust(8), h2.rjust(8), h3.rjust(8)
    outPDB += "ATOM     14  SG  CYS T   2    " + h1 + h2 + h3 + "  1.00  1.00           S" + "\n"

    k1, k2, k3 = str("%.3f" % TCRaPtxt[0]), str("%.3f" % TCRaPtxt[1]), str("%.3f" % TCRaPtxt[2])
    k1, k2, k3 = k1.rjust(8), k2.rjust(8), k3.rjust(8)
    outPDB += "ATOM     15  N   ALA U   1    " + k1 + k2 + k3 + "  1.00  1.00           S" + "\n"

    l1, l2, l3 = str("%.3f" % TCRbPtxt[0]), str("%.3f" % TCRbPtxt[1]), str("%.3f" % TCRbPtxt[2])
    l1, l2, l3 = l1.rjust(8), l2.rjust(8), l3.rjust(8)
    outPDB += "ATOM     16  CA  LEU U   2    " + l1 + l2 + l3 + "  1.00  1.00           S" + "\n"

    m1, m2, m3 = str("%.3f" % planeAtTCRatxt[0]), str("%.3f" % planeAtTCRatxt[1]), str("%.3f" % planeAtTCRatxt[2])
    m1, m2, m3 = m1.rjust(8), m2.rjust(8), m3.rjust(8)
    outPDB += "ATOM     17  O   VAL V   1    " + m1 + m2 + m3 + "  1.00  1.00           S" + "\n"

    n1, n2, n3 = str("%.3f" % planeAtTCRbtxt[0]), str("%.3f" % planeAtTCRbtxt[1]), str("%.3f" % planeAtTCRbtxt[2])
    n1, n2, n3 = n1.rjust(8), n2.rjust(8), n3.rjust(8)
    outPDB += "ATOM     18  N   ALA V   2    " + n1 + n2 + n3 + "  1.00  1.00           S" + "\n"

    outPDBfile.write(outPDB)
    outPDBfile.close()

    print "Done!", "\n", "Coordinate file written as ", fileName + "/crossingAngle/" + fileName + "_planeCorners.pdb\n"

    coordinateDict = {
        "corner1": "/corners//P/CYS`1/SG",
        "corner2": "/corners//P/CYS`2/SG",
        "corner3": "/corners//P/CYS`3/SG",
        "corner4": "/corners//P/CYS`4/SG",
        "TCRcentroida": "/corners//T/CYS`1/SG",
        "TCRcentroidb": "/corners//T/CYS`2/SG",
        "TCRcentroidaP": "/corners//U/ALA`1/N",
        "TCRcentroidbP": "/corners//U/LEU`2/CA",
        "MHCfita": "/corners//M/CYS`1/SG",
        "MHCfitb": "/corners//M/CYS`2/SG",
        "MHCafitP": "/corners//N/CYS`1/SG",
        "MHCbfitP": "/corners//N/CYS`2/SG",
        "MHCafitatTCR": "/corners//O/LEU`1/CA",
        "MHCbfitatTCR": "/corners//O/ALA`2/N",
        "MHCafitatTCRP": "/corners//Q/LEU`1/CA",
        "MHCbfitatTCRP": "/corners//Q/ALA`2/N",
        "planeAtTCRa": "/corners//V/VAL`1/O",
        "planeAtTCRb": "/corners//V/ALA`2/N"

    }

    complexDict = {
        "MHCa": MHCachain,
        "MHCb": MHCbchain,
        "p": peptidechain,
        "TCRa": TCRachain,
        "TCRb": TCRbchain
    }

    #################################### CALCULATIONS ##########################################
    print "\n ~~~~~~~~~~~ Bigging numerical calculations... ~~~~~~~~~~~\n"
    ######################  Calculate crossing angle (Rudolph 3D) ##############################

    print "By using MHC and TCR, we can calculate the crossing angle:"

    MHCV = MHC[0, :] - MHC[1, :]
    TCRV = TCR[1, :] - TCR[0, :]

    crossingAngleRudolph = math.degrees(angle_between(MHCV, TCRV))
    print crossingAngleRudolph

    print "\nIn order to determine the polarity of binding i.e. canonical (1 - 180 degrees) or reverse (180 - 360 degrees), we need to determine which quadrant the crossing angle lies in."
    print "\nThis is possible by moving all calculated points such that TCRB and MHCAatTCRB are both trannsformed to x,y,z = [0,0,0] and MHCBatTCRB sits along the x axis  to x,y,z = [20,0,0]"
    print "\nLet's check that stays the same if we move MHC points to MHCA at TCRB"

    MHCV2 = mhc_at_tcrb[0, :] - mhc_at_tcrb[1, :]
    TCRV2 = TCR[1, :] - TCR[0, :]

    print math.degrees(angle_between(MHCV2, TCRV2))

    ######### Determine polarity via alignment to origin using PyMOL ##########################
    print "\nMoving everything to the origin using PyMOL alignment to some target atoms..."

    pymol.cmd.delete("all")
    pymol.cmd.scene("*", "clear")
    # alignment
    pymol.cmd.load(fileName + "/crossingAngle/" + fileName + "_planeCorners.pdb", "corners")
    pymol.cmd.select("MHCmove", selection=coordinateDict["MHCafitatTCR"] + " + " + coordinateDict["MHCbfitatTCR"])
    pymol.cmd.load("bin/data/centroids_target.pdb", "centroids_target")
    pymol.cmd.align("MHCmove", "centroids_target")
    pymol.cmd.rotate("z", "180", camera=0)
    print "Aligned at origin!"

    # colour and show spheres and lines
    pymol.cmd.hide("everything")
    pymol.cmd.show("spheres", "MHCmove")
    pymol.cmd.select("TCRA", selection=coordinateDict["TCRcentroida"])
    pymol.cmd.color("red", "TCRA")
    pymol.cmd.select("TCRB", selection=coordinateDict["TCRcentroidb"])
    pymol.cmd.color("blue", "TCRB")
    pymol.cmd.distance("TCRline", "TCRA", "TCRB")
    pymol.cmd.select("MHCA", selection=coordinateDict["MHCafitatTCR"])
    pymol.cmd.color("green", "MHCA")
    pymol.cmd.select("MHCB", selection=coordinateDict["MHCbfitatTCR"])
    pymol.cmd.color("green", "MHCB")
    pymol.cmd.distance("MHCline", "MHCA", "MHCB")
    pymol.cmd.color("red", "TCRline")
    pymol.cmd.color("black", "MHCline")
    pymol.cmd.set("dash_gap", 0)

    ## Photo op here
    # set camera angle
    pymol.cmd.set_view(view_set.originView)
    # generate images
    print "Generating image.."
    if ray:
        centroidonxaxis = fileName + "/crossingAngle/" + fileName + "_crossingAngleOnAxis.png"
        rayTime(centroidonxaxis)

    # save the moved centroids

    pymol.cmd.save(fileName + "/crossingAngle/" + fileName + "_centroid_onxaxis.pdb")
    pymol.cmd.save(fileName + "/crossingAngle/" + fileName + "_centroid_onxaxis.pse")

    # Getting coords of TCRB
    TCRorig = np.vstack((pymol.cmd.get_coords("TCRA"), pymol.cmd.get_coords("TCRB")))
    MHCorig = np.vstack((pymol.cmd.get_coords("MHCA"), pymol.cmd.get_coords("MHCB")))

    pymol.cmd.delete("all")

    # Crossing angle calculator #
    print "By tracking the coordinates of TCRA, we can determine the polarity of binding. Originised points are.. \n"
    print "\nLine 1 (TCR) A:"
    print TCRorig[0]
    print "\nLine 1 (TCR) B:"
    print TCRorig[1]
    print "\nLine 2 (MHC) A:"
    print MHCorig[0]
    print "\nLine 2 (MHC) B:"
    print MHCorig[1]

    MHCV3 = MHCorig[0, :] - MHCorig[1, :]
    TCRV3 = TCRorig[1, :] - TCRorig[0, :]

    print "\nCrossing angle using originised points = ", math.degrees(angle_between(MHCV3, TCRV3))
    rangeLow, rangeHigh = whichQuadrant(TCRorig)
    smartCrossingAngle = whichCrossingAngle(crossingAngleRudolph, rangeLow)

    print "\n-----------------------\nCrossing angle (Rudolph 3D angle) = \n"
    print smartCrossingAngle, "\n-----------------------\n"

    ######################  Calculate 2D rotation angle  ##############################
    print "Now we can measure the angle between TCR and MHC on a 2D plane relative to the MHC surface..\n"
    MHCprojV = MHCproj[0, :] - MHCproj[1, :]
    TCRprojV = TCRproj[1, :] - TCRproj[0, :]

    print "MHCproj vector is..."
    print MHCprojV[0], "\n"
    print "TCRproj vector is..."
    print TCRprojV[0], "\n"

    print "Calculating angle.."
    crossingAngle2D = math.degrees(angle_between(MHCprojV, TCRprojV))
    print crossingAngle2D
    print "..Done!"

    print "Let's check that stays the same if we move MHCA point to TCRB.."

    MHCprojV2 = MHCatTCRBP[0, :] - MHCatTCRBP[1, :]
    TCRprojV2 = TCRproj[1, :] - TCRproj[0, :]
    print math.degrees(angle_between(MHCprojV2, TCRprojV2))

    print "\nAgain we must take into account the polarity of binding (which has already been determined):"

    smart2DAngle = whichCrossingAngle(crossingAngle2D, rangeLow)

    print "\n-----------------------\nRotation angle (2D crossing angle) = \n"
    print smart2DAngle, "\n-----------------------\n"

    ######################  Calculate tilt angle  ##############################
    print "Next we want to calculate the tilt!"
    print "To do this, we evaluate how far away TCR(beta) and TCRproj(beta) are \
    and then shift TCRproj coordinates (both alpha and beta) by this vector"

    planeAtTCRV = planeAtTCR[1, :] - planeAtTCR[0, :]
    print "Calculating angle.."
    tilt2D = math.degrees(angle_between(planeAtTCRV, TCRV))
    print "..Done!"

    print "Finally, we need to determine the direction of tilt: up (+) or down (-) from TCRB to TCRA."
    print "This is acheived by determining whether TCRA is closer to or further away to the plane than TCRB is.\n"

    BtoBP = np.asscalar(np.linalg.norm(TCR[1] - TCRproj[1]))
    AtoAP = np.asscalar(np.linalg.norm(TCR[0] - TCRproj[0]))

    print "TCRB distance from the plane is ", BtoBP
    print "TCRA distance from the plane is ", AtoAP

    if AtoAP >= BtoBP:
        print "Tilt angle is positive!"
        tilt2D = tilt2D * 1
    if AtoAP < BtoBP:
        print "Tilt angle is negative!"
        tilt2D = tilt2D * -1

    print "\n-----------------------\nTilt angle (2D elevation away from MHC) = \n"
    print tilt2D, "\n-----------------------\n"

    ########################## OUPUT WRITER ###########################################
    print "\nGenerating output file..\n"

    outTxt = ''
    outTxt += fileName + ".pdb\n\n"
    outTxt += "Crossing angle (Rudolph 3D angle)\t" + str("%.2f" % smartCrossingAngle) + "\n"
    outTxt += "Rotation angle (2D crossing angle)\t" + str("%.2f" % smart2DAngle) + "\n"
    outTxt += "Tilt angle (2D elevation away from MHC)\t" + str("%.2f" % tilt2D) + "\n"
    outTxt += "\n"

    outTxt += "\n\n"
    outTxt += "-------------\n\n"
    # Directionality parameters
    dxMHC = MHCbtxt[0] - MHCatxt[0]
    dyMHC = MHCbtxt[1] - MHCatxt[1]
    dzMHC = MHCbtxt[2] - MHCatxt[2]
    dxTCR = TCRbtxt[0] - TCRatxt[0]
    dyTCR = TCRbtxt[1] - TCRatxt[1]
    dzTCR = TCRbtxt[2] - TCRatxt[2]
    dxyzMHC = []
    dxyzTCR = []
    dxyzMHC.append(dxMHC)
    dxyzMHC.append(dyMHC)
    dxyzMHC.append(dzMHC)
    dxyzTCR.append(dxTCR)
    dxyzTCR.append(dyTCR)
    dxyzTCR.append(dzTCR)

    # Record coordinates with vectors in place #

    outTxt += "Coordinates of vector in position:\n"
    outTxt += "\tx\ty\tz\n"
    outTxt += "TCRA\t"
    for x in TCRatxt:
        outTxt += str("%.3f" % x)
        outTxt += "\t"
    outTxt += "\n"
    outTxt += "TCRB\t"
    for x in TCRbtxt:
        outTxt += str("%.3f" % x)
        outTxt += "\t"
    outTxt += "\n"
    outTxt += "MHCA\t"
    for x in MHCatxt:
        outTxt += str("%.3f" % x)
        outTxt += "\t"
    outTxt += "\n"
    outTxt += "MHCB\t"
    for x in MHCbtxt:
        outTxt += str("%.3f" % x)
        outTxt += "\t"
    outTxt += "\n\n"
    outTxt += "TCRB-A\t"
    for x in dxyzTCR:
        outTxt += str("%.3f" % x)
        outTxt += "\t"
    outTxt += "\n"
    outTxt += "MHCB-A\t"
    for x in dxyzMHC:
        outTxt += str("%.3f" % x)
        outTxt += "\t"
    outTxt += "\n\n"
    outTxt += "-------------\n\n"

    outTxt += "Crossing angle is in the " + str(rangeLow) + " to " + str(rangeHigh) + " degree range\n"
    if rangeHigh == 90 or rangeHigh == 180:
        outTxt += "TCR docks with canonical polarity\n\n"
    if rangeHigh == 270 or rangeHigh == 360:
        outTxt += "TCR docks with reverse polarity\n\n"
    outTxtFile = open(fileName + "/crossingAngle/" + fileName + '_crossingAngle.txt', 'w')
    outTxtFile.write(outTxt)

    print "Done!"

    ########## VISUALTISATION ####################################

    print "~~~~~~~~~~~ Creating visualisations of crossing angles in pymol... ~~~~~~~~~~~"

    # load and extract to objects
    pymol.cmd.load(fileName + "/crossingAngle/" + fileName + "_planeCorners.pdb", "corners")

    for key in coordinateDict:
        pymol.cmd.select(key + "s", selection=coordinateDict[key])
        pymol.cmd.extract(key, key + "s")
        pymol.cmd.delete(key + "s")

    pymol.cmd.load(fileName + "/crossingAngle/" + fileName + "_aligned.pdb", "complex")

    for key in complexDict:
        pymol.cmd.select(key + "s", selection="complex and chain " + complexDict[key])
        pymol.cmd.extract(key, key + "s")
        pymol.cmd.delete(key + "s")

    pymol.cmd.delete("corners")
    pymol.cmd.delete("complex")

    pymol.cmd.select("hetatoms", "hetatm")
    pymol.cmd.remove("hetatoms")

    pymol.cmd.hide("all")

    ######## Crossing Angle ############

    # show centroids
    pymol.cmd.show("spheres", "TCRcentroida")
    pymol.cmd.show("spheres", "TCRcentroidb")
    pymol.cmd.color(colour_set.generalcolour_set["TCRa"], "TCRcentroida")
    pymol.cmd.color(colour_set.generalcolour_set["TCRb"], "TCRcentroidb")

    # show line
    pymol.cmd.distance("TCRline", "TCRcentroida", "TCRcentroidb")
    pymol.cmd.hide("labels", "all")
    pymol.cmd.color("black", "TCRline")
    pymol.cmd.set("dash_gap", "0")

    # show MHC
    if mhc_class == "I":
        a1locs = range(50, 86)
        MHCa1h = ["MHCa"] + a1locs
        a2locs = range(140, 176)
        MHCa2h = ["MHCa"] + a2locs

    if mhc_class == "II":
        a1locs = range(46, 78)
        MHCa1h = ["MHCa"] + a1locs
        a2locs = range(54, 91)
        MHCa2h = ["MHCb"] + a2locs
    name = "MHCa1"
    locs = '+'.join(str(x) for x in MHCa1h[1:])
    pymol.cmd.select(name, selection=MHCa1h[0] + " and resi " + locs)
    name = "MHCa2"
    locs = '+'.join(str(x) for x in MHCa2h[1:])
    pymol.cmd.select(name, selection=MHCa2h[0] + " and resi " + locs)

    pymol.cmd.color(colour_set.generalcolour_set["MHCa"], "MHCa")
    pymol.cmd.color(colour_set.generalcolour_set["MHCb"], "MHCb")
    pymol.cmd.color(colour_set.generalcolour_set["p"], "p")
    pymol.cmd.set("transparency", 0.5)
    pymol.cmd.hide("lines", "all")
    pymol.cmd.show("surface", "MHCa")
    pymol.cmd.show("surface", "MHCb")
    pymol.cmd.show("surface", "p")
    pymol.cmd.show("cartoon", "MHCa1")
    pymol.cmd.show("cartoon", "MHCa2")
    pymol.cmd.set("cartoon_transparency", 0.5)

    # set view
    pymol.cmd.set_view(view_set.birdsEyeView)
    pymol.cmd.scene(key="crossing_angle", action="store")
    # generate images

    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_crossing_centroids_angle_pMHC.png"
        rayTime(scene_name)

    pymol.cmd.hide("cartoon", "all")
    pymol.cmd.hide("surface", "all")

    # Photo op
    pymol.cmd.set_view(view_set.birdsEyeView)
    pymol.cmd.scene(key="centroids_angle", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_crossing_centroids_angle.png"
        rayTime(scene_name)

    pymol.cmd.hide("dashes", "all")

    # Photo op
    pymol.cmd.set_view(view_set.birdsEyeView)
    pymol.cmd.scene(key="centroids", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_crossing_centroids.png"
        rayTime(scene_name)

    pymol.cmd.hide("spheres", "all")
    pymol.cmd.show("dashes", "TCRline")

    # Photo op
    pymol.cmd.set_view(view_set.birdsEyeView)
    pymol.cmd.scene(key="angle", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_crossing_angle.png"
        rayTime(scene_name)

    # save session
    pymol.cmd.save(fileName + "/crossingAngle/" + fileName + "_crossing_angle.pse")
    print "\nOutputted PyMOL session file: " + fileName + "/crossingAngle/" + fileName + "_crossing_angle.pse"

    ############# Rotation Angle ###########################

    # import plane
    # pymol.cmd.extend("plane", plane.plane)
    # pymol.cmd.extend("make_plane", plane.make_plane)
    pymol.cmd.scene("*", "clear")
    pymol.cmd.hide("all")

    pymol.cmd.distance("plane1", "corner1", "corner2")
    pymol.cmd.distance("plane2", "corner3", "corner4")
    pymol.cmd.distance("plane3", "corner2", "corner4")
    pymol.cmd.distance("plane4", "corner1", "corner3")

    pymol.cmd.distance("MHCfitPline", "MHCafitP", "MHCbfitP")
    pymol.cmd.distance("TCRPline", "TCRcentroidaP", "TCRcentroidbP")
    pymol.cmd.distance("MHCfitatTCRPline", "MHCafitatTCRP", "MHCbfitatTCRP")

    pymol.cmd.hide("labels", "all")
    pymol.cmd.show("spheres", "TCRcentroidaP")
    pymol.cmd.show("spheres", "TCRcentroidbP")

    pymol.cmd.color(colour_set.generalcolour_set["TCRa"], "TCRcentroidaP")
    pymol.cmd.color(colour_set.generalcolour_set["TCRb"], "TCRcentroidbP")

    pymol.cmd.color("yellow", "plane1")
    pymol.cmd.color("yellow", "plane2")
    pymol.cmd.color("yellow", "plane3")
    pymol.cmd.color("yellow", "plane4")

    pymol.cmd.color("gray70", "MHCfitPline")
    pymol.cmd.color("black", "TCRPline")
    pymol.cmd.color("gray30", "MHCfitatTCRPline")

    pymol.cmd.show("surface", "MHCa")
    pymol.cmd.show("surface", "MHCb")
    pymol.cmd.show("surface", "p")
    pymol.cmd.show("cartoon", "MHCa1")
    pymol.cmd.show("cartoon", "MHCa2")

    ## Photo op here
    pymol.cmd.set_view(view_set.birdsEyeView)
    pymol.cmd.zoom("center", 80)
    pymol.cmd.scene(key="rotation_angle_pMHC_1", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_rotation_angle_pMHC.png"
        rayTime(scene_name)

    ## Photo op here
    pymol.cmd.set_view(view_set.offAxis)
    pymol.cmd.scene(key="rotation_angle_pMHC_2", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_rotation_angle_pMHC_2.png"
        rayTime(scene_name)

    pymol.cmd.hide("cartoon", "all")
    pymol.cmd.hide("surface", "all")

    ## Photo op here
    pymol.cmd.set_view(view_set.birdsEyeView)
    pymol.cmd.zoom("center", 80)
    pymol.cmd.scene(key="rotation_angle_centroids_1", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_rotation_angle_centroids_1.png"
        rayTime(scene_name)

    ## Photo op here
    pymol.cmd.set_view(view_set.offAxis)
    pymol.cmd.scene(key="rotation_angle_centroids_2", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_rotation_angle_centroids_2.png"
        rayTime(scene_name)

    pymol.cmd.hide("dashes", "MHCfitPline")
    pymol.cmd.hide("dashes", "TCRPline")
    pymol.cmd.hide("dashes", "MHCfitatTCRPline")

    ## Photo op here
    pymol.cmd.set_view(view_set.birdsEyeView)
    pymol.cmd.zoom("center", 80)
    pymol.cmd.scene(key="rotation_centroids", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_rotation_centroids_1.png"
        rayTime(scene_name)

    ## Photo op here
    pymol.cmd.set_view(view_set.offAxis)
    pymol.cmd.scene(key="rotation_centroids_2", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_rotation_centroids_2.png"
        rayTime(scene_name)

    pymol.cmd.hide("spheres", "all")
    pymol.cmd.show("dashes", "MHCfitPline")
    pymol.cmd.show("dashes", "TCRPline")
    pymol.cmd.show("dashes", "MHCfitatTCRPline")

    ## Photo op here
    pymol.cmd.set_view(view_set.birdsEyeView)
    pymol.cmd.zoom("center", 80)
    pymol.cmd.scene(key="rotation_angle", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_rotation_angle_1.png"
        rayTime(scene_name)

    ## Photo op here
    pymol.cmd.set_view(view_set.offAxis)
    pymol.cmd.scene(key="rotation_angle_2", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_rotation_angle_2.png"
        rayTime(scene_name)
    # save session
    pymol.cmd.save(fileName + "/crossingAngle/" + fileName + "_rotation_angle.pse")
    print "\nOutputted PyMOL session file: " + fileName + "/crossingAngle/" + fileName + "_rotation_angle.pse"

    ############# Tilt Angle ###########################
    pymol.cmd.scene("*", "clear")
    pymol.cmd.hide("all")

    pymol.cmd.distance("planeTCRline", "planeAtTCRa", "planeAtTCRb")
    pymol.cmd.show("dashes", "TCRline")
    pymol.cmd.show("spheres", "TCRcentroida")
    pymol.cmd.show("spheres", "TCRcentroidb")
    pymol.cmd.color("gray50", "planeTCRline")
    pymol.cmd.set("transparency", 0.5)
    pymol.cmd.hide("labels", "all")

    pymol.cmd.show("cartoon", "TCRa")
    pymol.cmd.show("cartoon", "TCRb")

    pymol.cmd.color(colour_set.generalcolour_set["TCRa"], "TCRa")
    pymol.cmd.color(colour_set.generalcolour_set["TCRb"], "TCRb")

    pymol.cmd.show("surface", "MHCa")
    pymol.cmd.show("surface", "MHCb")
    pymol.cmd.show("surface", "p")
    pymol.cmd.show("cartoon", "MHCa1")
    pymol.cmd.show("cartoon", "MHCa2")

    pymol.cmd.show("dashes", "plane1")
    pymol.cmd.show("dashes", "plane2")
    pymol.cmd.show("dashes", "plane3")
    pymol.cmd.show("dashes", "plane4")
    ## Photo op here
    pymol.cmd.set_view(view_set.offAxis)
    pymol.cmd.scene(key="tilt_angle_TCR_pMHC_1", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_tilt_angle_TCR_pMHC_1.png"
        rayTime(scene_name)

    pymol.cmd.hide("cartoon", "TCRa")
    pymol.cmd.hide("cartoon", "TCRb")

    ## Photo op here
    pymol.cmd.set_view(view_set.offAxis)
    pymol.cmd.scene(key="tilt_angle_pMHC_1", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_tilt_angle_pMHC_1.png"
        rayTime(scene_name)

    pymol.cmd.hide("surface", "all")
    pymol.cmd.hide("cartoon", "all")

    ## Photo op here
    pymol.cmd.set_view(view_set.offAxis)
    pymol.cmd.scene(key="tilt_angle_centroids_1", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_tilt_angle_centroids_1.png"
        rayTime(scene_name)

    pymol.cmd.hide("dashes", "all")

    ## Photo op here
    pymol.cmd.set_view(view_set.offAxis)
    pymol.cmd.scene(key="tilt_centroids_1", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_tilt_centroids_1.png"
        rayTime(scene_name)

    pymol.cmd.hide("spheres", "all")
    pymol.cmd.show("dashes", "planeTCRline")
    pymol.cmd.show("dashes", "TCRline")
    pymol.cmd.show("dashes", "plane1")
    pymol.cmd.show("dashes", "plane2")
    pymol.cmd.show("dashes", "plane3")
    pymol.cmd.show("dashes", "plane4")

    ## Photo op here
    pymol.cmd.set_view(view_set.offAxis)
    pymol.cmd.scene(key="tilt_angle_1", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_tilt_angle_1.png"
        rayTime(scene_name)

    # save session
    pymol.cmd.save(fileName + "/crossingAngle/" + fileName + "_tilt_angle.pse")
    print "\nOutputted PyMOL session file: " + fileName + "/crossingAngle/" + fileName + "_tilt_angle.pse"

    ##### Tilt angle 2 (at aligned view such that TCR vector travels in x direction) #####

    pymol.cmd.scene("*", "clear")
    pymol.cmd.hide("all")
    # move to common location to take image with planeAtTCRb flat across image
    pymol.cmd.load("bin/data/centroids_target.pdb", "centroid_target")
    pymol.cmd.create("mobile", "TCRcentroidaP or TCRcentroidbP")
    pymol.cmd.align("mobile", "centroid_target")

    for key in coordinateDict:
        pymol.cmd.matrix_copy("mobile", key)
    for key in complexDict:
        pymol.cmd.matrix_copy("mobile", key)

    # pymol.cmd.remove("centroid_target")

    pymol.cmd.distance("planeTCRline", "planeAtTCRa", "planeAtTCRb")
    pymol.cmd.show("dashes", "TCRline")
    pymol.cmd.show("dashes", "planeTCRline")
    pymol.cmd.show("spheres", "TCRcentroida")
    pymol.cmd.show("spheres", "TCRcentroidb")
    pymol.cmd.color("gray50", "planeTCRline")
    pymol.cmd.set("transparency", 0.5)
    pymol.cmd.hide("labels", "all")

    pymol.cmd.show("cartoon", "TCRa")
    pymol.cmd.show("cartoon", "TCRb")

    pymol.cmd.color(colour_set.generalcolour_set["TCRa"], "TCRa")
    pymol.cmd.color(colour_set.generalcolour_set["TCRb"], "TCRb")

    pymol.cmd.show("surface", "MHCa")
    pymol.cmd.show("surface", "MHCb")
    pymol.cmd.show("surface", "p")
    pymol.cmd.show("cartoon", "MHCa1")
    pymol.cmd.show("cartoon", "MHCa2")

    pymol.cmd.show("dashes", "plane1")
    pymol.cmd.show("dashes", "plane2")
    pymol.cmd.show("dashes", "plane3")
    pymol.cmd.show("dashes", "plane4")

    ## Photo op here
    pymol.cmd.orient("mobile")
    pymol.cmd.zoom("center", 80)
    pymol.cmd.turn("x", -45)
    if rangeLow == 0 or rangeLow == 90:
        pymol.cmd.turn("x", -45)
    if rangeLow == 180 or rangeLow == 270:
        pymol.cmd.turn("x", +45)

    pymol.cmd.scene(key="tilt_angle_TCR_pMHC_2", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_tilt_angle_TCR_pMHC_2.png"
        rayTime(scene_name)

    pymol.cmd.hide("cartoon", "TCRa")
    pymol.cmd.hide("cartoon", "TCRb")

    ## Photo op here
    pymol.cmd.scene(key="tilt_angle_pMHC_2", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_tilt_angle_pMHC_2.png"
        rayTime(scene_name)

    pymol.cmd.hide("surface", "all")
    pymol.cmd.hide("cartoon", "all")

    ## Photo op here
    pymol.cmd.scene(key="tilt_angle_centroids_2", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_tilt_angle_centroids_2.png"
        rayTime(scene_name)

    pymol.cmd.hide("dashes", "all")

    ## Photo op here
    pymol.cmd.scene(key="tilt_centroids_2", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_tilt_centroids_2.png"
        rayTime(scene_name)

    pymol.cmd.hide("spheres", "all")
    pymol.cmd.show("dashes", "planeTCRline")
    pymol.cmd.show("dashes", "TCRline")
    pymol.cmd.show("dashes", "plane1")
    pymol.cmd.show("dashes", "plane2")
    pymol.cmd.show("dashes", "plane3")
    pymol.cmd.show("dashes", "plane4")

    ## Photo op here
    pymol.cmd.scene(key="tilt_angle_2", action="store")
    if ray:
        scene_name = fileName + "/crossingAngle/" + fileName + "_tilt_angle_2.png"
        rayTime(scene_name)

    # save session
    pymol.cmd.save(fileName + "/crossingAngle/" + fileName + "_tilt_angle_2.pse")
    print "\nOutputted PyMOL session file: " + fileName + "/crossingAngle/" + fileName + "_tilt_angle_2.pse"

    # planeProperties = {'ALPHA':0.6, 'COLOR':[0.55, 0.25, 0.60], 'INVERT':False}
    # plane.make_plane("MHCplane","corner2","corner3", "corner4", center = True, settings=planeProperties)

    # plane.plane("corner1","corner2", "corner3", "corner3", settings=planeProperties)

    print "\nDone!\n"

    # Clean up #
    outTxtFile.close()
    pymol.cmd.quit()
    os.remove(fileName + "/crossingAngle/" + fileName + "_aligned_noMeta.pdb")
    print "\ncrossingAngle.py created the following files in the directory " + fileName + "/crossingAngle" + ":"
    for files in os.listdir(os.getcwd() + "/" + fileName + "/crossingAngle"):
        print files

    print('     ~  End crossingAngle.py v1.9 BETA  ~')