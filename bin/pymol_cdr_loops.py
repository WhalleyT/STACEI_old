import os
import pymol
import sys
import argparse
import time
import bin.data.colourSet as colourSet
import bin.data.viewSet as viewSet
import numpy

description = \
    "foo"


### File loader ###

def parse_args():
    """parse arguments for command line"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--PDB', dest='pdb', type=str, required=True, help='The PDB file to be analysed')
    parser.add_argument('--MHCclass', dest='MHCclass', type=str, required=True,
                        help='The MHC class i.e. class I or II of the complex structure')
    parser.add_argument('--Fasta', dest='fasta', type=str, required=True, default=None,
                        help='The fasta file that contains information about where the CDR loops are in the structure to be analysed')
    parser.add_argument('--Chains', dest='chains', type=str, required=False,
                        help='Chains of TCR-pMHC complex in order MHCa,MHCb,peptide,TCRa,TCRb', default="ABCDE")
    parser.add_argument('--Ray', dest='ray', type=str, required=False,
                        help='Do you want to render images. Set to False for speed and dev. Default = True',
                        default=True)

    args = parser.parse_args()
    return args


### File loader ###

def readFile(FILE, fileType):
    if FILE == None:
        sys.exit('No file to read')
    if FILE.split('.')[-1].lower() != str(fileType):
        sys.exit("\nFile extension must be of type ." + str(fileType) + "\n")
    else:
        print 'Reading file: ' + str(FILE)
        return open(FILE, "r")


def fileToList(inFile):
    allLines = []
    for line in inFile:
        allLines.append(line)
    return allLines


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


def purgeCysLocs(locations):
    output = []
    for locs in locations:
        if "Cys" not in locs[0]:
            output.append(locs)
    return output


def initialisePymol():
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


############ BODY #########################

### Initialiser ###

# Load input.pdb #
def generate(pdb, fasta, MHCclass, chains, ray, fileName):
    print('     ~  Running autoCDRloops.py v0.1 BETA  ~')

    PDBfile = readFile(pdb, "pdb")
    FASTAfile = readFile(fasta, "fasta")

    if MHCclass == 1:
        MHCclass = "I"
    else:
        MHCclass = "II"

    if not os.path.exists(fileName):
        print "Creating Directory " + fileName
        os.makedirs(fileName)

    if not os.path.exists(fileName + "/visualisation"):
        print "Creating Directory " + fileName + "/visualisation"
        os.makedirs(fileName + "/visualisation")

    if not os.path.exists(fileName + "/pdbs"):
        print "Creating Directory " + fileName + "/pdbs"
        os.makedirs(fileName + "/pdbs")

    # Unpack the fasta file with CDR information #

    fastaEntries = fastaParser(fasta)
    fastaEntries = depackID(fastaEntries)

    TCRA = []
    for entry in fastaEntries:
        if "TCRA" in entry:
            TCRA = entry
    TCRAlocations = findLocations(TCRA)
    TCRAlocations = depackLocations(TCRAlocations)
    TCRAlocations = purgeCysLocs(TCRAlocations)

    TCRB = []
    for entry in fastaEntries:
        if "TCRB" in entry:
            TCRB = entry
    TCRBlocations = findLocations(TCRB)
    TCRBlocations = depackLocations(TCRBlocations)
    TCRBlocations = purgeCysLocs(TCRBlocations)

    # Sort chains
    MHCachain, MHCbchain, peptidechain, TCRachain, TCRbchain = chains[0], chains[1], chains[2], chains[3], chains[4]

    # Find the MHC helices
    if MHCclass == "I":
        a1locs = range(50, 86)
        MHCa1 = ["MHCa"] + a1locs
        a2locs = range(140, 176)
        MHCa2 = ["MHCa"] + a2locs

    if MHCclass == "II":
        a1locs = range(46, 78)
        MHCa1 = ["MHCa"] + a1locs
        a2locs = range(54, 91)
        MHCa2 = ["MHCb"] + a2locs

    # Let's get started

    initialisePymol()
    pymol.cmd.load(fileName + "/crossingAngle/" + fileName + "_aligned.pdb", "complex")

    # Make chains objects
    pymol.cmd.select("MHCas", selection="chain " + MHCachain)
    pymol.cmd.extract("MHCa", "MHCas")
    pymol.cmd.select("MHCbs", selection="chain " + MHCbchain)
    pymol.cmd.extract("MHCb", "MHCbs")
    pymol.cmd.select("ps", selection="chain " + peptidechain)
    pymol.cmd.extract("p", "ps")
    pymol.cmd.select("TCRas", selection="chain " + TCRachain)
    pymol.cmd.extract("TCRa", "TCRas")
    pymol.cmd.select("TCRbs", selection="chain " + TCRbchain)
    pymol.cmd.extract("TCRb", "TCRbs")

    pymol.cmd.select("hetatoms", "hetatm")
    pymol.cmd.remove("hetatoms")

    # delete selections
    pymol.cmd.delete("complex")
    pymol.cmd.delete("MHCas")
    pymol.cmd.delete("MHCbs")
    pymol.cmd.delete("ps")
    pymol.cmd.delete("TCRas")
    pymol.cmd.delete("TCRbs")

    # General colours

    complexObjects = ["MHCa", "MHCb", "p", "TCRa", "TCRb"]

    for objs in complexObjects:
        pymol.cmd.color(colourSet.generalColourSet[objs], objs)

    pymol.cmd.hide("cartoon", "p")
    pymol.cmd.show("sticks", "p")
    pymol.cmd.util.cnc("p")

    # Colour the CDR loops
    for loop in TCRAlocations:
        name = loop[0]
        locs = '+'.join(str(x) for x in loop[1:])
        pymol.cmd.select(name + "s", selection="TCRa and resi " + locs)
        pymol.cmd.color(colourSet.CDRcolourSet[name], name + "s")

    for loop in TCRBlocations:
        name = loop[0]
        locs = '+'.join(str(x) for x in loop[1:])
        pymol.cmd.select(name + "s", selection="TCRb and resi " + locs)
        pymol.cmd.color(colourSet.CDRcolourSet[name], name + "s")

    # Select the MHCa helices

    name = "MHCa1"
    locs = '+'.join(str(x) for x in MHCa1[1:])
    pymol.cmd.select(name, selection=MHCa1[0] + " and resi " + locs)
    name = "MHCa2"
    locs = '+'.join(str(x) for x in MHCa2[1:])
    pymol.cmd.select(name, selection=MHCa2[0] + " and resi " + locs)

    # General front scene
    pymol.cmd.set_view(viewSet.newFrontView)
    pymol.cmd.scene(key="front", action="store")
    if ray == True:
        frontViewImage = fileName + "/visualisation/" + "front.png"
        rayTime(frontViewImage)

    # General side scene
    pymol.cmd.set_view(viewSet.sideView)
    pymol.cmd.scene(key="side", action="store")
    if ray == True:
        sideViewImage = fileName + "/visualisation/" + "side.png"
        rayTime(sideViewImage)

    # Birds eye with CDR loops over MHC

    pymol.cmd.hide("all")

    for loop in TCRAlocations:
        name = loop[0]
        pymol.cmd.create(name, name + "s")
        pymol.cmd.show("ribbon", name)

    for loop in TCRBlocations:
        name = loop[0]
        pymol.cmd.create(name, name + "s")
        pymol.cmd.show("ribbon", name)

    # Photo op. Let's image just the loops alone
    pymol.cmd.set_view(viewSet.newBirdsEyeView)
    pymol.cmd.scene(key="CDRloops", action="store")
    if ray == True:
        CDRloopsImage = fileName + "/visualisation/" + "CDRloops.png"
        rayTime(CDRloopsImage)

    # pMHC surface
    pymol.cmd.show("surface", "MHCa")
    pymol.cmd.show("surface", "MHCb")
    pymol.cmd.show("surface", "p")
    pymol.cmd.set("transparency", 0.5)
    pymol.cmd.color("gray40", "p")

    # pMHC helices

    pymol.cmd.show("cartoon", "MHCa1")
    pymol.cmd.show("cartoon", "MHCa2")
    pymol.cmd.set("cartoon_transparency", 0.5)

    ## Photo op here. CDR birds eye with CDR ribbons on top of surface pMHC shoing cartoon helices

    pymol.cmd.set_view(viewSet.newBirdsEyeView)
    pymol.cmd.scene(key="CDRbirdseye", action="store")
    if ray == True:
        CDRbirdseyeout = fileName + "/visualisation/" + "CDRloopsBirdsEye.png"
        rayTime(CDRbirdseyeout)

    ## Now let's add some surface to the CDR loops up the birdseye
    for loop in TCRAlocations:
        name = loop[0]
        if "FW" not in name:
            pymol.cmd.show("surface", name)

    for loop in TCRBlocations:
        name = loop[0]
        if "FW" not in name:
            pymol.cmd.show("surface", name)

    ## Photo op here. CDR birds eye with CDR ribbons and surface on top of surface pMHC showing cartoon helices
    pymol.cmd.set_view(viewSet.newBirdsEyeView)
    pymol.cmd.scene(key="CDRsurfaceBirdseye", action="store")
    if ray == True:
        CDRsurfaceBirdseye = fileName + "/visualisation/" + "CDRsurfaceBirdseye.png"
        rayTime(CDRsurfaceBirdseye)
    # Now let's show the surface of pMHC with coloured contacts by the CDR loops

    # Hide the loop surfaces
    for loop in TCRAlocations:
        name = loop[0]
        pymol.cmd.hide("surface", name)
    for loop in TCRBlocations:
        name = loop[0]
        pymol.cmd.hide("surface", name)

    # first select them
    for loop in TCRAlocations:
        name = loop[0]
        pymol.cmd.select(name + "_near", selection="MHCa or MHCb or p within 4.0 of " + name)
        pymol.cmd.color(colourSet.CDRcolourSet[name], name + "_near")

    for loop in TCRBlocations:
        name = loop[0]
        pymol.cmd.select(name + "_near", selection="MHCa or MHCb or p within 4.0 of " + name)
        pymol.cmd.color(colourSet.CDRcolourSet[name], name + "_near")

    ## Photo op here
    pymol.cmd.set_view(viewSet.newBirdsEyeView)
    pymol.cmd.scene(key="CDRfootprint", action="store")
    if ray == True:
        CDRfootprint = fileName + "/visualisation/" + "CDRfootprint.png"
        rayTime(CDRfootprint)

    # Now let's just view the pMHC surface as this might be useful for overlays

    # Hide the loop ribbons
    for loop in TCRAlocations:
        name = loop[0]
        pymol.cmd.hide("ribbon", name)
    for loop in TCRBlocations:
        name = loop[0]
        pymol.cmd.hide("ribbon", name)

    # Return the colours back to norman
    for objs in complexObjects:
        pymol.cmd.color(colourSet.generalColourSet[objs], objs)
    pymol.cmd.copy("p2", "p")
    pymol.cmd.show("sticks", "p2")
    pymol.cmd.util.cnc("p2")

    ## Photo op here
    pymol.cmd.set_view(viewSet.newBirdsEyeView)
    pymol.cmd.scene(key="pMHCsurface", action="store")
    if ray == True:
        pMHCsurface = fileName + "/visualisation/" + "pMHCsurface.png"
        rayTime(pMHCsurface)

    # MHC helices cos pymol is being a pain
    pymol.cmd.hide("all")
    pymol.cmd.show("cartoon", "MHCa1")
    pymol.cmd.show("cartoon", "MHCa2")
    pymol.cmd.set("cartoon_transparency", 0.5)

    ## Photo op here
    pymol.cmd.set_view(viewSet.newBirdsEyeView)
    pymol.cmd.scene(key="MHChelices", action="store")
    if ray == True:
        MHChelices = fileName + "/visualisation/" + "MHChelices.png"
        rayTime(MHChelices)

    # CDR (backbone) centre of mass

    pymol.cmd.hide("all")

    COMs = []
    for loop in TCRAlocations:
        COM = []
        name = loop[0]
        COM.append(name)
        loopCOM = pymol.cmd.centerofmass(name + " and n. c+o+n+ca")
        COM.append(loopCOM)
        COMs.append(COM)
        pymol.cmd.pseudoatom(name + "_COM", pos=loopCOM)
        pymol.cmd.show("spheres", name + "_COM")
        pymol.cmd.color(colourSet.CDRcolourSet[name], name + "_COM")
        pymol.cmd.set("sphere_scale", 1.5, name + "_COM")
        pymol.cmd.save(fileName + "/pdbs/" + fileName + "_" + name + "_COM.pdb", name + "_COM")

    for loop in TCRBlocations:
        COM = []
        name = loop[0]
        COM.append(name)
        loopCOM = pymol.cmd.centerofmass(name + " and n. c+o+n+ca")
        COM.append(loopCOM)
        COMs.append(COM)
        pymol.cmd.pseudoatom(name + "_COM", pos=loopCOM)
        pymol.cmd.show("spheres", name + "_COM")
        pymol.cmd.color(colourSet.CDRcolourSet[name], name + "_COM")
        pymol.cmd.set("sphere_scale", 1.5, name + "_COM")
        pymol.cmd.save(fileName + "/pdbs/" + fileName + "_" + name + "_COM.pdb", name + "_COM")

    for COM in COMs:
        print COM

    ## Photo op here
    pymol.cmd.set_view(viewSet.newBirdsEyeView)
    pymol.cmd.scene(key="CDRCOM", action="store")
    if ray == True:
        CDRCOMimage = fileName + "/visualisation/" + "CDR_centre_of_mass.png"
        rayTime(CDRCOMimage)

    # Save the session
    pymol.cmd.save(fileName + "/sessions/" + fileName + "_autoCDRloops.pse")


    # for mean in loopMeans:
    #    print mean

    # Save the loops

    for loop in TCRAlocations:
        name = loop[0]
        pymol.cmd.save(fileName + "/pdbs/" + fileName + "_" + name + ".pdb", name)

    for loop in TCRBlocations:
        name = loop[0]
        pymol.cmd.save(fileName + "/pdbs/" + fileName + "_" + name + ".pdb", name)

    # Quit pymol
    #pymol.cmd.quit()
    print('     ~  End autoCDRloops.py v0.1 BETA  ~')


