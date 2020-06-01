import numpy
import os
import pymol
import sys
import time
import data.colourSet
import data.viewSet


def read_file(infile, file_type):
    if infile is None:
        sys.exit('No file to read')

    if infile.split('.')[-1].lower() != str(file_type):
        sys.exit("\nFile extension must be of type ." + str(file_type) + "\n")
    else:
        print('Reading file: ' + str(infile))
        return open(infile, "r")


def file_to_list(infile):
    all_lines = []
    for line in infile:
        all_lines.append(line)
    return all_lines


def align_to_template(file_name, mhc_class):
    print(os.getcwd()); input()
    print("\nAligning file to template...\n")
    pymol.finish_launching(['pymol', '-qeim'])
    pymol.cmd.load(file_name + ".pdb")
    pymol.cmd.load("bin/data/" + mhc_class + "_template.pdb")
    print("x")
    pymol.cmd.align(file_name, mhc_class + "_template")
    print("y")
    pymol.cmd.save(file_name + "/crossingAngle/" + file_name + "_aligned_noMeta.pdb", file_name)
    pymol.cmd.reinitialize()
    print("\nAlignment to " + mhc_class + "_template.pdb complete!\n")
    return None


def sulphide_parser(infile):
    print("Finding SSBOND Lines...")
    out_ss = []
    for line in infile:
        add = ''
        if "SSBOND" in line[0:6]:
            add += line
            out_ss.append(add)
    return out_ss


def header_parser(inFile):
    print("Finding header...")
    outHeader = []
    for line in inFile:
        add = ''
        if "HEADER" in line[0:6]:
            add += line
            outHeader.append(add)
    return outHeader


def title_parser(infile):
    print("Finding title...")
    outTitle = []
    for line in infile:
        add = ''
        if "TITLE" in line[0:6]:
            add += line
            outTitle.append(add)
    return outTitle


### Ray tools ###

def wait_for_ray(query):
    counter = 0
    while not os.path.exists(query):
        print(counter)
        time.sleep(1)
        counter += 1
    return None


def save_ray(save_as):
    print("Outputting image.. This may take a few seconds..")
    if os.path.exists(save_as):
        print("Removing " + save_as + " as it already exists!")
        os.remove(save_as)
    time.sleep(10)
    pymol.cmd.png(save_as, ray=1, width=3000, height=3000, dpi=300)
    wait_for_ray(save_as)
    print("Done! " + str(save_as) + " was outputted")


### PDB Atom Parser ###

def atom_parser(all_lines):
    print(str(len(all_lines)) + " went in")
    atom_list = []
    else_list = []
    for line in all_lines:
        if "ATOM" in line[0:5]:
            atom_list.append(line)
        else:
            else_list.append(line)
    print(str(len(atom_list)) + " atom lines detected")
    print(str(len(else_list)) + " else lines detected")
    return atom_list, else_list


def atomRowParser(atomRow):
    atomRowParsed = [atomRow[0:6], int(atomRow[6:11]), atomRow[12:16], atomRow[16], atomRow[17:20], atomRow[21],
                     int(atomRow[22:26]), atomRow[26], float(atomRow[30:38]), float(atomRow[38:46]),
                     float(atomRow[46:54]), float(atomRow[54:60]), float(atomRow[60:66]), atomRow[76:78],
                     atomRow[78:-1]]
    return atomRowParsed


def make_atom_matrix(atomList):
    atomMatrix = []
    print('\nCreating atom matrix...')
    for row in atomList:
        atomMatrix.append(atomRowParser(row))
    print("Done!\n")
    return atomMatrix


### Robust Cysteine pair finder ###

def sulphide_bond_parser(infile):
    print("\nFinding disulphide bridges...")
    ss_bond_list = []
    for line in infile:
        if "SSBOND" in line[0:6]:
            parsedLine = [line[15], int(line[18:21]), line[29], int(line[32:35])]
            ss_bond_list.append(parsedLine)
    print(str(len(ss_bond_list)) + " SSBOND lines detected!\n")
    if len(ss_bond_list) > 0:
        return ss_bond_list
    else:
        return False


def tcr_pair_finder(ss_bond_list):
    """
    Takes in all the SSBonds and outputs two lists.. 
    
    the alpha chain cys and the beta chain cys
    """
    print("Finding TCRa and TCRb cysteine pairs...\n ")
    alphaCys = []
    betaCys = []
    for row in ss_bond_list:
        if row[0] == 'D':
            if 15 <= row[1] <= 35:
                if row[2] == 'D':
                    if 85 <= row[3] <= 110:
                        alphaCys = row
        if row[0] == 'D':
            if 15 <= row[3] <= 35:
                if row[2] == 'D':
                    if 85 <= row[1] <= 110:
                        alphaCys = row
    for row in ss_bond_list:
        if row[0] == 'E':
            if 15 <= row[1] <= 35:
                if row[2] == 'E':
                    if 85 <= row[3] <= 110:
                        betaCys = row
        if row[0] == 'E':
            if 15 <= row[3] <= 35:
                if row[2] == 'E':
                    if 85 <= row[1] <= 110:
                        betaCys = row
    print(alphaCys)
    print(betaCys)
    print("\n")
    if len(alphaCys) > 0 and len(betaCys) > 0:
        return alphaCys, betaCys
    else:
        return False, False


def tcr_pair_atom_finder(pair_location, atom_matrix):
    print("Finding SG atoms...")
    coordinate1 = []
    coordinate2 = []
    chain = pair_location[0]
    cysResidue1 = pair_location[1]
    cysResidue2 = pair_location[3]
    for atoms in atom_matrix:
        if atoms[5] == chain:
            if atoms[6] == cysResidue1:
                if "SG" in atoms[2]:
                    if atoms[3] == '' or 'A':
                        coordinate1 = atoms
    for atoms in atom_matrix:
        if atoms[5] == chain:
            if atoms[6] == cysResidue2:
                if "SG" in atoms[2]:
                    if atoms[3] == '' or 'A':
                        coordinate2 = atoms
    return coordinate1, coordinate2


def get_sg_coords(row):
    x = row[8]
    y = row[9]
    z = row[10]
    SGcoords = [x, y, z]
    return SGcoords


def robust_cysteine(infile, atom_matrix):
    ss_lines = sulphide_bond_parser(infile)
    if not ss_lines:
        print("No SSBOND lines were detected in the PDB!\n")
        return False, False, False, False
    a, b = tcr_pair_finder(ss_lines)
    if a is False or b is False:
        print("Could not find TCR cysteine pairs in the SSBOND list!\n")
        return False, False, False, False
    a1_atoms, a2_atoms = tcr_pair_atom_finder(a, atom_matrix)
    b1_atoms, b2_atoms = tcr_pair_atom_finder(b, atom_matrix)
    if a1_atoms == False or a2_atoms == False or b1_atoms == False or b2_atoms == False:
        print("Could not find TCR cysteine atoms from the cysteine pairs listed in SSBOND list!\n")
        return False, False, False, False
    return a1_atoms, a2_atoms, b1_atoms, b2_atoms


# Weak Cys pair finder #    

def cysPairFinder(atomMatrix, chain):
    print("Finding CYS residues between residues 15 to 35 and between 85 to 105...")
    cysteines = []
    cyspairs = []
    cysPairsSG = []
    for row in atomMatrix:
        if row[5] == chain and "CYS" in row[4]:
            cysteines.append(row)
    for row in cysteines:
        if 15 <= row[6] <= 35 or 85 <= row[6] <= 105:
            if row[3] == '' or 'A':
                cyspairs.append(row)
    for x in cyspairs:
        print(x)
    print('\nNumber of cys atoms found in chain ' + str(chain) + ' around the S-S bridge region: ' + str(len(cyspairs)))
    for row in cyspairs:
        if 'SG' in row[2]:
            cysPairsSG.append(row)
    print('\nNumber of cys SG atoms found in chain ' + str(chain) + ' around the S-S bridge region: ' + str(
        len(cysPairsSG)) + "\n")
    for x in cysPairsSG:
        print(x)
    print("\n")
    return cysPairsSG


def cysWeak(atomMatrix):
    a = cysPairFinder(atomMatrix, 'D')
    a1, a2 = a[0], a[1]
    b = cysPairFinder(atomMatrix, 'E')
    b1, b2 = b[0], b[1]
    return a1, a2, b1, b2


### MHC coordinates finder ###

def MHCaxisI(atomMatrix):
    print("\nMHC axis atoms are:\n")
    Calphas = []
    helixCalphas = []
    for row in atomMatrix:
        if row[5] == 'A' and "CA" in row[2]:
            Calphas.append(row)
    for row in Calphas:
        if 50 <= row[6] <= 86 or 140 <= row[6] <= 176:
            if row[3] == '' or 'A':
                helixCalphas.append(row)
    for x in helixCalphas:
        print(x)
    print('\nNumber of Calpha atoms in MHC axis: ' + str(len(helixCalphas)) + "\n")
    return helixCalphas


def MHCaxisII(atomMatrix):
    print("\nMHC axis atoms are:\n")
    Calphas = []
    helixCalphas = []
    for row in atomMatrix:
        if row[5] == 'A' and "CA" in row[2]:
            Calphas.append(row)
    for row in Calphas:
        if 46 <= row[6] <= 78:
            if row[3] == '' or 'A':
                helixCalphas.append(row)
    Calphas = []
    for row in atomMatrix:
        if row[5] == 'B' and "CA" in row[2]:
            Calphas.append(row)
    for row in Calphas:
        if 54 <= row[6] <= 64:
            if row[3] == '' or 'A':
                helixCalphas.append(row)
        if 67 <= row[6] <= 91:
            if row[3] == '' or 'A':
                helixCalphas.append(row)
    for x in helixCalphas:
        print(x)
    print('\nNumber of Calpha atoms in MHC axis: ' + str(len(helixCalphas)) + "\n")
    return helixCalphas


def MHCaxisCoords(helixCalphas):
    print("Extracting MHC axis coordinates...")
    MHCx = []
    MHCy = []
    MHCz = []
    for row in helixCalphas:
        MHCx.append(row[8])
        MHCy.append(row[9])
        MHCz.append(row[10])
    MHCcoords = [MHCx, MHCy, MHCz]
    return MHCcoords


def findMHCaxisCoords(atomMatrix, MHCclass):
    if MHCclass == "I":
        print('Establishing MHC I axis...')
        return MHCaxisCoords(MHCaxisI(atomMatrix))
    if MHCclass == "II":
        print('Establishing MHC II axis...')
        return MHCaxisCoords(MHCaxisII(atomMatrix))


def generateMatrixM(mhc_3d_coords):
    mhc_x_bar = sum(mhc_3d_coords[0]) / float(len(mhc_3d_coords[0]))
    mhc_y_bar = sum(mhc_3d_coords[1]) / float(len(mhc_3d_coords[1]))
    mhc_z_bar = sum(mhc_3d_coords[2]) / float(len(mhc_3d_coords[2]))
    mhc_x_minus_xbar = []
    for values in mhc_3d_coords[0]:
        mhc_x_minus_xbar.append(values - mhc_x_bar)
    mhc_y_minus_ybar = []
    for values in mhc_3d_coords[1]:
        mhc_y_minus_ybar.append(values - mhc_y_bar)
    mhc_z_minus_zbar = []
    for values in mhc_3d_coords[2]:
        mhc_z_minus_zbar.append(values - mhc_z_bar)
    matrix_mlist = [mhc_x_minus_xbar, mhc_y_minus_ybar, mhc_z_minus_zbar]
    return numpy.array(matrix_mlist, dtype=float)


### Matrix A ###

def generateMatrixA(matrixM):
    print("\nGenerating matrix A...\n")
    matrixMt = matrixM.transpose()
    matrixMtt = matrixMt.transpose()
    print("Matrix A: \n")
    print(numpy.dot(matrixMtt, matrixMt))
    return numpy.dot(matrixMtt, matrixMt)


### MHC vector ###

def generateEigens(matrixA):
    eigenvecval = numpy.linalg.eigh(matrixA)
    eigenvalues = eigenvecval[0]
    eigenvectors = eigenvecval[1]
    print("\nLargest Eigenvalue:")
    print(eigenvalues[-1])
    print("\nLargest Eigenvector")
    print(eigenvectors[:, 2] * -1)
    print("\n")
    return eigenvalues[-1], eigenvectors[:, 2] * -1


def crossingAngleCalculator(TCRA, TCRB, MHCA, MHCB):
    print("\nCalculating crossing angle...\n")
    p = (TCRB[0] - TCRA[0]) * (MHCB[0] - MHCA[0]) + (TCRB[1] - TCRA[1]) * (MHCB[1] - MHCA[1]) + (TCRB[2] - TCRA[2]) * (
        MHCB[2] - MHCA[2])
    q = ((TCRB[0] - TCRA[0]) ** 2) + ((TCRB[1] - TCRA[1]) ** 2) + ((TCRB[2] - TCRA[2]) ** 2)
    r = ((MHCB[0] - MHCA[0]) ** 2) + ((MHCB[1] - MHCA[1]) ** 2) + ((MHCB[2] - MHCA[2]) ** 2)
    s = p / ((q ** 0.5) * (r ** 0.5))
    crossingAngle = (180 / numpy.pi) * numpy.arccos(s)
    return crossingAngle


def which_quadrant(TCRBx, TCRBy):
    if TCRBx > 0 and TCRBy > 0:
        print("Crossing angle is in the 0-90 quadrant\n")
        return 0, 90
    elif TCRBx < 0 < TCRBy:
        print("Crossing angle is in the 90-180 quadrant\n")
        return 90, 180
    elif TCRBx < 0 and TCRBy < 0:
        print("Crossing angle is in the 180-270 quadrant\n")
        return 180, 270
    elif TCRBx > 0 > TCRBy:
        print("Crossing angle is in the 270-360 quadrant\n")
        return 270, 360
    else:
        print("Crossing angle is not within any expected quadrant\n")
        sys.exit()


def which_crossing_angle(xing_angle, range_low):
    if range_low == 0 or range_low == 90:
        return xing_angle
    if range_low == 180 or range_low == 270:
        return 360 - xing_angle


def pyMOLparameters():
    # set PyMOL parameters
    pymol.cmd.set("ray_shadows", "0")
    pymol.cmd.set("specular", "off")
    pymol.cmd.set("orthoscopic", "on")
    pymol.cmd.bg_color("white")
    pymol.cmd.set("ray_opaque_background", "0")


def crossingAngleVisualisation(fileName, MHCclass, ray):
    pyMOLparameters()
    print("Generating crossing angle image")
    # load and extract to objects
    pymol.cmd.load(fileName + "/crossingAngle/" + fileName + "_centroids.pdb", "centroids")
    pymol.cmd.select("TCRaCentroid", selection="chain Y")
    pymol.cmd.extract("TCRaCentroido", "TCRaCentroid")
    pymol.cmd.select("TCRbCentroid", selection="chain Z")
    pymol.cmd.extract("TCRbCentroido", "TCRbCentroid")
    pymol.cmd.delete("centroids")
    pymol.cmd.load(fileName + "/crossingAngle/" + fileName + "_aligned.pdb", "complex")
    pymol.cmd.select("MHCa", selection="chain A")
    pymol.cmd.extract("MHCao", "MHCa")
    pymol.cmd.select("MHCb", selection="chain B")
    pymol.cmd.extract("MHCbo", "MHCb")
    pymol.cmd.select("p", selection="chain C")
    pymol.cmd.extract("po", "p")
    pymol.cmd.select("TCRa", selection="chain D")
    pymol.cmd.extract("TCRao", "TCRa")
    pymol.cmd.select("TCRb", selection="chain E")
    pymol.cmd.extract("TCRbo", "TCRb")
    # delete selections


    pymol.cmd.delete("complex")
    pymol.cmd.delete("MHCa")
    pymol.cmd.delete("MHCb")
    pymol.cmd.delete("p")
    pymol.cmd.delete("TCRaCentroid")
    pymol.cmd.delete("TCRbCentroid")
    pymol.cmd.delete("TCRa")
    pymol.cmd.delete("TCRb")

    pymol.cmd.select("hetatoms", "hetatm")
    pymol.cmd.remove("hetatoms")

    pymol.cmd.hide("cartoon", "all")
    # show centroids
    pymol.cmd.show("spheres", "TCRaCentroido")
    pymol.cmd.show("spheres", "TCRbCentroido")
    pymol.cmd.color(data.colourSet.generalColourSet["TCRa"], "TCRaCentroido")
    pymol.cmd.color(data.colourSet.generalColourSet["TCRb"], "TCRbCentroido")
    pymol.cmd.distance("TCRline", "TCRaCentroid", "TCRbCentroid")
    pymol.cmd.hide("labels", "all")
    pymol.cmd.color("black", "TCRline")
    pymol.cmd.set("dash_gap", "0")

    # Photo op. Let's image just the loops alone
    pymol.cmd.set_view(data.viewSet.birdsEyeView)
    pymol.cmd.scene(key="Centroids", action="store")
    if ray:
        crossingAngleCentroidsImage = fileName + "/crossingAngle/" + fileName + "_centroidsCrossingAngle.png"
        save_ray(crossingAngleCentroidsImage)

    # show MHC
    if MHCclass == "I":
        range_a1loc = list(range(50, 86))
        mhca_1 = ["MHCao"] + range_a1loc
        a2locs = list(range(140, 176))
        mhca_2 = ["MHCao"] + a2locs

    if MHCclass == "II":
        range_a1loc = list(range(46, 78))
        mhca_1 = ["MHCao"] + range_a1loc
        a2locs = list(range(54, 91))
        mhca_2 = ["MHCbo"] + a2locs

    name = "MHCa1"
    locs = '+'.join(str(x) for x in mhca_1[1:])
    pymol.cmd.select(name, selection=mhca_1[0] + " and resi " + locs)
    name = "MHCa2"
    locs = '+'.join(str(x) for x in mhca_2[1:])
    pymol.cmd.select(name, selection=mhca_2[0] + " and resi " + locs)

    pymol.cmd.color(data.colourSet.generalColourSet["MHCa"], "MHCao")
    pymol.cmd.color(data.colourSet.generalColourSet["MHCb"], "MHCbo")
    pymol.cmd.color(data.colourSet.generalColourSet["p"], "po")
    pymol.cmd.set("transparency", 0.5)
    pymol.cmd.hide("lines", "all")
    pymol.cmd.show("surface", "MHCao")
    pymol.cmd.show("surface", "MHCbo")
    pymol.cmd.show("surface", "po")
    pymol.cmd.show("cartoon", "MHCa1")
    pymol.cmd.show("cartoon", "MHCa2")
    pymol.cmd.set("cartoon_transparency", 0.5)

    # set view
    pymol.cmd.set_view(data.viewSet.birdsEyeView)
    pymol.cmd.scene(key="CentroidsMHC", action="store")
    # generate images

    if ray:
        crossingAngleImage = fileName + "/crossingAngle/" + fileName + "_crossingangle.png"
        save_ray(crossingAngleImage)

        # save sessiona and quit
    pymol.cmd.save(fileName + "/crossingAngle/" + fileName + "_aligned_crossingangle.pse")
    print("Outputted " + fileName + "_crossingangle.png image and " + fileName + \
          "_aligned_crossingangle.pse session file.\n")
    pymol.cmd.reinitialize()


def angleToTargetPyMOL(file_name, ray):
    pyMOLparameters()
    # alignment
    pymol.cmd.load(file_name + "/crossingAngle/" + file_name + "_centroids.pdb", "aligned_centroids")
    pymol.cmd.select("MHCmove", selection="chain M+N")
    pymol.cmd.load("bin/data/centroids_target.pdb", "centroids_target")
    pymol.cmd.align("MHCmove", "centroids_target")
    print("Aligned at origin!")
    # colour and show spheres and lines
    pymol.cmd.hide("everything")
    pymol.cmd.show("spheres")
    pymol.cmd.select("TCRA", selection="aligned_centroids and chain Y")
    pymol.cmd.color("red", "TCRA")
    pymol.cmd.select("TCRB", selection="aligned_centroids and chain Z")
    pymol.cmd.color("blue", "TCRB")
    pymol.cmd.distance("TCRline", "TCRA", "TCRB")
    pymol.cmd.select("MHCA", "aligned_centroids and chain M")
    pymol.cmd.color("green", "MHCA")
    pymol.cmd.select("MHCB", "aligned_centroids and chain N")
    pymol.cmd.color("green", "MHCB")
    pymol.cmd.distance("MHCline", "MHCA", "MHCB")
    pymol.cmd.color("red", "TCRline")
    pymol.cmd.color("black", "MHCline")
    pymol.cmd.set("dash_gap", 0)

    # set camera angle
    pymol.cmd.set_view("\
         0.997469127,    0.005865565,    0.070841491,\
        -0.057437014,    0.653679490,    0.754587770,\
        -0.041880410,   -0.756749690,    0.652363896,\
        -0.000028916,   -0.000041917, -166.161453247,\
        11.795754433,    8.656616211,    2.305841923,\
      -26788.419921875, 27120.732421875,    0.000000000")
    # generate images
    print("Generating image..")

    if ray:
        centroid_on_axis = file_name + "/crossingAngle/" + file_name + "_centroidonxaxis.png"
        save_ray(centroid_on_axis)

    # save the moved centroids

    pymol.cmd.save(file_name + "/crossingAngle/" + file_name + "_centroid_onxaxis.pdb")
    pymol.cmd.save(file_name + "/crossingAngle/" + file_name + "_centroid_onxaxis.pse")

    # Getting coords of TCRB
    tcra_array = pymol.cmd.get_coords("TCRA")
    tcra_list = (tcra_array.tolist()[0])

    tcrb_array = pymol.cmd.get_coords("TCRB")
    tcrb_list = (tcrb_array.tolist()[0])

    mhca_array = pymol.cmd.get_coords("MHCA")
    mhca_list = (mhca_array.tolist()[0])

    mhcb_array = pymol.cmd.get_coords("MHCB")
    mhcb_list = (mhcb_array.tolist()[0])
    pymol.cmd.reinitialize()
    return tcra_list, tcrb_list, mhca_list, mhcb_list


def calculate(pdb, mhc_class, ray, file_name):
    origPDB = read_file(pdb, "pdb")

    if mhc_class == 1:
        mhc_class = "I"
    if mhc_class == 2:
        mhc_class = "II"

    if type(file_name) != str:
        sys.exit('No file to grab metadata was loaded. Please view usage and provide a valid .pdb to be processed')

    # Align input.pdb to template #

    align_to_template(file_name, mhc_class)
    print("Opening aligned PDB file")
    alignedPDBnoMeta = read_file(file_name + "/crossingAngle/" + file_name + "_aligned_noMeta.pdb", "pdb")

    # Extract all lines from PDB file of orig and aligned PDB files #

    origall_lines = file_to_list(origPDB)
    alignedAllLines = file_to_list(alignedPDBnoMeta)
    alignedPDBnoMeta.close()

    ### Metadata Replacer ###

    print("Transfering metadata from " + file_name + ".pdb to " + file_name + "_aligned_noMeta\n")
    metatitle = title_parser(origall_lines)
    metaheader = header_parser(origall_lines)
    MetaSS = sulphide_parser(origall_lines)
    working_file = open(file_name + "/crossingAngle/" + file_name + '_aligned.pdb', 'w')
    outTxtMeta = ''
    outTxtMeta += "         " + str(1) + "         " + str(2) + "         " + str(3) + "         " + str(4) + \
                  "         " + str(5) + "         " + str(6) + "         " + str(7) + "         " + str(8) + "\n"
    outTxtMeta += str(12345678901234567890123456789012345678901234567890123456789012345678901234567890) + "\n"
    outTxtMeta += "REMARK   1                                                                      \n\
    REMARK   1  ~~~~ Metadata generated by crossingAngleFull_v1.0.py     ~~~~ \n\
    REMARK   1                                                                      \n\
    REMARK   1  This metadata contains :           \n\
    REMARK   1      HEADER                        \n\
    REMARK   1      TITLE         \n\
    REMARK   1      This REMARK \n\
    REMARK   1      SSBOND lines required for crossing angle calculations    \n\
    REMARK   1\n"
    for Metax in metaheader:
        outTxtMeta += Metax
    for Metax in metatitle:
        outTxtMeta += Metax
    for Metax in MetaSS:
        outTxtMeta += Metax
    for Metax in alignedAllLines:
        outTxtMeta += Metax
    working_file.write(outTxtMeta)
    print("Metadata added to " + file_name + "/crossingAngle/" + file_name + "_aligned.pdb\n")

    # Load the working PDB file
    print("Extracting all lines from " + file_name + "/crossingAngle/" + file_name + "_aligned.pdb\n")
    working_file = open(file_name + "/crossingAngle/" + file_name + '_aligned.pdb', 'r')
    all_lines = file_to_list(working_file)
    working_file.close()
    atomList, elseList = atom_parser(all_lines)
    atomM = make_atom_matrix(atomList)

    ### Cysteine axis ###

    cysWARNING = False
    TCRaCys1, TCRaCys2, TCRbCys1, TCRbCys2 = robust_cysteine(elseList, atomM)
    if TCRaCys1 == False or TCRaCys2 == False or TCRbCys1 == False or TCRbCys2 == False:
        cysWARNING = True
        TCRaCys1, TCRaCys2, TCRbCys1, TCRbCys2 = cysWeak(atomM)
    print("\nTCRa and TCRb cysteine pair SG atoms are:\n")
    print(TCRaCys1)
    print(TCRaCys2)
    print(TCRbCys1)
    print(TCRbCys2)
    cysAlpha1 = get_sg_coords(TCRaCys1)
    cysAlpha2 = get_sg_coords(TCRaCys2)
    cysBeta1 = get_sg_coords(TCRbCys1)
    cysBeta2 = get_sg_coords(TCRbCys2)
    print("\nTCRa and TCRb cysteine pair SG coordinates are [x,y,z]:\n")
    print(cysAlpha1)
    print(cysAlpha2)
    print(cysBeta1)
    print(cysBeta1)
    tcra = [(cysAlpha1[0] + cysAlpha2[0]) / 2, (cysAlpha1[1] + cysAlpha2[1]) / 2, (cysAlpha1[2] + cysAlpha2[2]) / 2]
    tcrb = [(cysBeta1[0] + cysBeta2[0]) / 2, (cysBeta1[1] + cysBeta2[1]) / 2, (cysBeta1[2] + cysBeta2[2]) / 2]
    print("\nTCRa cysteine centroid is [x,y,z]:")
    print(tcra)
    print("\nTCRb cysteine centroid is [x,y,z]:")
    print(tcrb)
    print("\n\n")

    # MHC axis #

    MHCcoords = findMHCaxisCoords(atomM, mhc_class)
    matrixM = generateMatrixM(MHCcoords)
    matrixA = generateMatrixA(matrixM)
    Evalue, Evector = generateEigens(matrixA)

    #################################### Visualisation ###############################################
    # These are visualisation parameters only #

    scalingFactor = 20
    MHCA = tcra
    MHCBx = MHCA[0] + scalingFactor * (Evector[0] * -1)
    MHCBy = MHCA[1] + scalingFactor * (Evector[1] * -1)
    MHCBz = MHCA[2] + scalingFactor * (Evector[2] * -1)
    MHCB = [MHCBx, MHCBy, MHCBz]

    # Round everything to 3dp #

    tcra = [round(elem, 3) for elem in tcra]
    tcrb = [round(elem, 3) for elem in tcrb]
    MHCA = [round(elem, 3) for elem in MHCA]
    MHCB = [round(elem, 3) for elem in MHCB]

    # Centroid PDB file #

    print("Generating centroid PDB file for visualisation...\n")
    centroidsPdbFile = open(file_name + "/crossingAngle/" + file_name + '_centroids.pdb', 'w')
    outPDB = ''
    outPDB += "         " + str(1) + "         " + str(2) + "         " + str(3) + "         " + str(4) + \
              "         " + str(5) + "         " + str(6) + "         " + str(7) + "         " + str(8) + "\n"
    outPDB += str(12345678901234567890123456789012345678901234567890123456789012345678901234567890) + "\n"

    outPDB += "REMARK   1                                                                      \n\
    REMARK   1  ~~~~         Generated by crossingAngleFull v1.0.py            ~~~~ \n\
    REMARK   1                                                                      \n\
    REMARK   1  This PDB File contains four arbitrary atoms:                         \n\
    REMARK   1      TCRa centroid in chain D                                        \n\
    REMARK   1      TCRb centroid in chain E                                        \n\
    REMARK   1      MHC axis point A fixed at TCRa                                  \n\
    REMARK   1      MHC axis point B extrapolated from MHC   A                      \n\
    REMARK   1                                                                      \n"

    a1, a2, a3 = str("%.3f" % tcra[0]), str("%.3f" % tcra[1]), str("%.3f" % tcra[2])
    a1, a2, a3 = a1.rjust(8), a2.rjust(8), a3.rjust(8)
    outPDB += "ATOM      1  SG  CYS Y   1    " + a1 + a2 + a3 + "  1.00  1.00           S" + "\n"

    b1, b2, b3 = str("%.3f" % tcrb[0]), str("%.3f" % tcrb[1]), str("%.3f" % tcrb[2])
    b1, b2, b3 = b1.rjust(8), b2.rjust(8), b3.rjust(8)
    outPDB += "ATOM      2  SG  CYS Z   1    " + b1 + b2 + b3 + "  1.00  1.00           S""\n"

    c1, c2, c3 = str("%.3f" % MHCA[0]), str("%.3f" % MHCA[1]), str("%.3f" % MHCA[2])
    c1, c2, c3 = c1.rjust(8), c2.rjust(8), c3.rjust(8)
    outPDB += "ATOM      3  CA  LEU M   1    " + c1 + c2 + c3 + "  1.00  1.00           S""\n"

    d1, d2, d3 = str("%.3f" % MHCB[0]), str("%.3f" % MHCB[1]), str("%.3f" % MHCB[2])
    d1, d2, d3 = d1.rjust(8), d2.rjust(8), d3.rjust(8)
    outPDB += "ATOM      4  N   ALA N   1    " + d1 + d2 + d3 + "  1.00  1.00           S"
    centroidsPdbFile.write(outPDB)
    centroidsPdbFile.close()
    print("Done!\n")
    # Make me that image #

    crossingAngleVisualisation(file_name, mhc_class, ray)

    # Directionality parameters
    print("Angle calculator output: ")
    angleInposition = crossingAngleCalculator(tcra, tcrb, MHCA, MHCB)
    print(("%.2f" % angleInposition))
    dxMHC = MHCB[0] - MHCA[0]
    dyMHC = MHCB[1] - MHCA[1]
    dzMHC = MHCB[2] - MHCA[2]
    dxTCR = tcrb[0] - tcra[0]
    dyTCR = tcrb[1] - tcra[1]
    dzTCR = tcrb[2] - tcra[2]
    dxyzMHC = []
    dxyzTCR = []
    dxyzMHC.append(dxMHC)
    dxyzMHC.append(dyMHC)
    dxyzMHC.append(dzMHC)
    dxyzTCR.append(dxTCR)
    dxyzTCR.append(dyTCR)
    dxyzTCR.append(dzTCR)

    # Record coordinates with vectors in place #
    outTxt = ''
    outTxt += file_name + ".pdb\nCoordinates of vector in position:\n"
    outTxt += "\tx\ty\tz\n"
    outTxt += "TCRA\t"
    for x in tcra:
        outTxt += str("%.3f" % x)
        outTxt += "\t"
    outTxt += "\n"
    outTxt += "TCRB\t"
    for x in tcrb:
        outTxt += str("%.3f" % x)
        outTxt += "\t"
    outTxt += "\n"
    outTxt += "MHCA\t"
    for x in MHCA:
        outTxt += str("%.3f" % x)
        outTxt += "\t"
    outTxt += "\n"
    outTxt += "MHCB\t"
    for x in MHCB:
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
    outTxt += "Angle calculator output: "
    outTxt += ("%.3f" % angleInposition)
    outTxt += "\n\n-------------"

    #################################### Calculation ###############################################
    # These are calculation parameters only #

    print("Moving everything to the origin...")
    MHCA = [0.0, 0.0, 0.0]
    tcra, tcrb, MHCA, MHCB = angleToTargetPyMOL(file_name, ray)
    tcra = [round(elem, 3) for elem in tcra]
    tcrb = [round(elem, 3) for elem in tcrb]
    print(tcrb)
    MHCA = [round(elem, 3) for elem in MHCA]
    MHCB = [round(elem, 3) for elem in MHCB]

    # Crossing angle calculator #

    print("\nLine 1 (TCR) A:")
    print(tcra)
    print("\nLine 1 (TCR) B:")
    print(tcrb)
    print("\nLine 2 (MHC) A:")
    print(MHCA)
    print("\nLine 2 (MHC) B:")
    print(MHCB)
    crossingAngle = crossingAngleCalculator(tcra, tcrb, MHCA, MHCB)

    # Directionality parameters

    dxMHC = MHCB[0] - MHCA[0]
    dyMHC = MHCB[1] - MHCA[1]
    dzMHC = MHCB[2] - MHCA[2]
    dxTCR = tcrb[0] - tcra[0]
    dyTCR = tcrb[1] - tcra[1]
    dzTCR = tcrb[2] - tcra[2]
    dxyzMHC = []
    dxyzTCR = []
    dxyzMHC.append(dxMHC)
    dxyzMHC.append(dyMHC)
    dxyzMHC.append(dzMHC)
    dxyzTCR.append(dxTCR)
    dxyzTCR.append(dyTCR)
    dxyzTCR.append(dzTCR)

    # Record coordinates with vectors at origin #
    outTxt += "\n\nCoordinates of vector at origin and x axis:\n"
    outTxt += "\tx\ty\tz\n"
    outTxt += "TCRA\t"
    for x in tcra:
        outTxt += ("%.3f" % x)
        outTxt += "\t"
    outTxt += "\n"
    outTxt += "TCRB\t"
    for x in tcrb:
        outTxt += str("%.3f" % x)
        outTxt += "\t"
    outTxt += "\n"
    outTxt += "MHCA\t"
    for x in MHCA:
        outTxt += str("%.3f" % x)
        outTxt += "\t"
    outTxt += "\n"
    outTxt += "MHCB\t"
    for x in MHCB:
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

    # Quadrant finder #

    TCRBx = tcrb[0]
    TCRBy = tcrb[1]
    print("TCRBx is " + str(TCRBx))
    print("TCRBy is " + str(TCRBy) + "\n")
    rangeLow, rangeHigh = which_quadrant(TCRBx, TCRBy)
    print("Output from angle calculator = " + ("%.2f" % crossingAngle))
    smartCrossingAngle = which_crossing_angle(crossingAngle, rangeLow)

    ########################## ANSWER ###########################################
    print("\n#########    ANSWER    ############\n")
    print("The crossing angle is:")
    print(("%.2f" % smartCrossingAngle))

    cysWarningmessage = "WARNING!!!!!.. this value was calculated by attempting to find the disulphide bridge in " \
                        "each TCR chain and may cause errors if there are CYS residues within the CDR3 region.\n \
                        Please check the output to confirm that the 'TCRa and TCRb cysteine pair SG atoms' " \
                        "are indeed forming a disulphide bridge in your structure.\n"

    cysOkaymessage = "Calculated using the robust Cys detection method\n"
    if cysWARNING:
        print(cysWarningmessage)
        outTxt += cysWarningmessage

    if not cysWARNING:
        print(cysOkaymessage)
        outTxt += cysOkaymessage

    ### OUTPUT ANSWER ###

    outTxt += "Crossing angle is in the " + str(rangeLow) + " to " + str(rangeHigh) + " degree range\n"
    if rangeHigh == 90 or rangeHigh == 180:
        outTxt += "TCR docks with canonical polarity\n\n"
    if rangeHigh == 270 or rangeHigh == 360:
        outTxt += "TCR docks with reverse polarity\n\n"
    outTxt += "Crossing angle = " + ("%.2f" % smartCrossingAngle)
    outTxtFile = open(file_name + "/crossingAngle/" + file_name + '_crossingAngle.txt', 'w')
    outTxtFile.write(outTxt)

    # Clean uo #
    #pymol.cmd.quit()
    outTxtFile.close()
    os.remove(file_name + "/crossingAngle/" + file_name + "_aligned_noMeta.pdb")
    print("\ncrossingAngle.py created the following files in the directory " + file_name + "/crossingAngle" + ":")
    for files in os.listdir(os.getcwd() + "/" + file_name + "/crossingAngle"):
        print(files)
    print("Quitting")
    #pymol.cmd.quit()
    print("done")