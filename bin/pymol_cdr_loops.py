import os
import pymol
import sys
import time
import bin.data.colourSet as colourSet
import bin.data.viewSet as viewSet

def readFile(FILE, fileType):
    if FILE == None:
        sys.exit('No file to read')
    if FILE.split('.')[-1].lower() != str(fileType):
        sys.exit("\nFile extension must be of type ." + str(fileType) + "\n")
    else:
        print('Reading file: ' + str(FILE))
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
        location += locationstring.split(',')
        locations.append(location)
    return locations


def purgeCysLocs(locations):
    output = []
    for locs in locations:
        if "Cys" not in locs[0]:
            output.append(locs)
    return output


def initialisePymol():
    print("\nInitialising pymol...\n")
    pymol.finish_launching(['pymol', '-qeimc'])
    pymol.cmd.reinitialize()
    # set PyMOL parameters
    pymol.cmd.set("ray_shadows", "0")
    pymol.cmd.set("specular", "off")
    pymol.cmd.set("orthoscopic", "on")
    pymol.cmd.bg_color("white")
    pymol.cmd.set("ray_opaque_background", "0")
    return None


def wait4ray(query):
    print("Waiting for image to render...")
    while not os.path.exists(query):
        time.sleep(1)
    return None

def rayTime(saveas, tracing):
    print("Outputting image.. This may take a few seconds..")
    if os.path.exists(saveas):
        print("Removing " + saveas + " as it already exists!")
        os.remove(saveas)
    time.sleep(2)
    pymol.cmd.png(saveas, ray=tracing, width=3000, height=3000, dpi=300)
    wait4ray(saveas)
    print("Done! " + str(saveas) + " was outputted")


############ BODY #########################

### Initialiser ###

# Load input.pdb #
def generate(pdb, fasta, MHCclass, chains, ray, fileName, id):
    print('     ~  Running autoCDRloops.py v0.1 BETA  ~')

    PDBfile = readFile(pdb, "pdb")
    FASTAfile = readFile(fasta, "fasta")

    if MHCclass == 1:
        MHCclass = "I"
    else:
        MHCclass = "II"

    if not os.path.exists(fileName):
        print("Creating Directory " + fileName)
        os.makedirs(fileName)

    if not os.path.exists(fileName + "/visualisation"):
        print("Creating Directory " + fileName + "/visualisation")
        os.makedirs(fileName + "/visualisation")

    if not os.path.exists(fileName + "/pdbs"):
        print("Creating Directory " + fileName + "/pdbs")
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

    # Find the MHC helices and groove
    if MHCclass == "I":
        a1locs = list(range(50, 86))
        MHCa1 = ["chain " + MHCachain ] + a1locs
        a2locs = list(range(140, 176))
        MHCa2 = ["chain " + MHCachain ] + a2locs

    if MHCclass == "II":
        a1locs = list(range(46, 78))
        MHCa1 = ["chain " + MHCachain ] + a1locs
        a2locs = list(range(54, 91))
        MHCa2 = ["chain " + MHCbchain ] + a2locs

    # Let's get started

    initialisePymol()
    pymol.cmd.load(fileName + "/crossingAngle/" + id + "_aligned.pdb", "complex")

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


    locs = '+'.join(str(x) for x in MHCa1[1:])
    sele = MHCa1[0] + " and resi " + locs
    pymol.cmd.select("MHCa1", sele)

    locs = '+'.join(str(x) for x in MHCa2[1:])
    sele = MHCa2[0] + " and resi " + locs
    pymol.cmd.select("MHCa2", sele)


    pymol.cmd.show("cartoon", "MHCa1")
    pymol.cmd.show("cartoon", "MHCa2")

    # General front scene
    pymol.cmd.set_view(viewSet.frontView)
    pymol.cmd.scene(key="front", action="store")

    frontViewImage = fileName + "/visualisation/" + id + "_front.png"

    if ray:
        rayTime(frontViewImage, 1)
    else:
        rayTime(frontViewImage, 0)

    # General side scene
    pymol.cmd.set_view(viewSet.sideView)
    pymol.cmd.scene(key="side", action="store")

    sideViewImage = fileName + "/visualisation/" + id + "_side.png"

    if ray:
        rayTime(sideViewImage, 1)
    else:
        rayTime(sideViewImage, 0)

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
    pymol.cmd.set_view(viewSet.birdsEyeView)
    pymol.cmd.scene(key="CDRloops", action="store")

    CDRloopsImage = fileName + "/visualisation/" + id + "_CDRloops.png"
    
    if ray:
        rayTime(CDRloopsImage, 1)
    else:
        rayTime(CDRloopsImage, 0)

    # Bruce edited here 16/06/19
    # pMHC surface
    if MHCclass == "I":
        pymol.cmd.select("MHCgrooves", selection="MHCa and resi " + '+'.join(str(x) for x in range(1,176)))
    if MHCclass == "II":
        pymol.cmd.select("MHCgrooves", selection="MHCa and resi "+ '+'.join(str(x) for x in range(1,78))+" or MHCb and resi "+ '+'.join(str(x) for x in range(1,91)))
    
    pymol.cmd.extract("MHCgroove", "MHCgrooves")
    pymol.cmd.delete("MHCgrooves")
    pymol.cmd.show("surface", "MHCgroove")
    
    if MHCclass == "I":
        pymol.cmd.color(colourSet.generalColourSet["MHCa"], "MHCgroove and chain "+MHCachain)
    
    if MHCclass == "II": 
        pymol.cmd.color(colourSet.generalColourSet["MHCa"], "MHCgroove and chain "+MHCachain)
        pymol.cmd.color(colourSet.generalColourSet["MHCb"], "MHCgroove and chain "+MHCbchain)
        
    # End Bruce edit 16/06/19
    
    pymol.cmd.show("surface", "p")
    pymol.cmd.set("transparency", 0.5)
    pymol.cmd.color("gray40", "p")

    # pMHC helices

    pymol.cmd.show("cartoon", "MHCa1")
    pymol.cmd.show("cartoon", "MHCa2")
    pymol.cmd.set("cartoon_transparency", 0.5)

    ## Photo op here. CDR birds eye with CDR ribbons on top of surface pMHC shoing cartoon helices

    pymol.cmd.set_view(viewSet.birdsEyeView)
    pymol.cmd.scene(key="CDRbirdseye", action="store")

    CDRbirdseyeout = fileName + "/visualisation/" + id  + "_CDRloopsBirdsEye.png"
    
    if ray:
        rayTime(CDRbirdseyeout, 1)
    else:
        rayTime(CDRbirdseyeout, 0)    

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
    pymol.cmd.set_view(viewSet.birdsEyeView)
    pymol.cmd.scene(key="CDRsurfaceBirdseye", action="store")

    CDRsurfaceBirdseye = fileName + "/visualisation/" + id + "_CDRsurfaceBirdseye.png"

    if ray:
        rayTime(CDRsurfaceBirdseye, 1)
    else:
        rayTime(CDRsurfaceBirdseye, 0)   
    
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
    pymol.cmd.set_view(viewSet.birdsEyeView)
    pymol.cmd.scene(key="CDRfootprint", action="store")
    
    CDRfootprint = fileName + "/visualisation/" + id + "_CDRfootprint.png"
    
    if ray:
        rayTime(CDRfootprint, 1)
    else:
        rayTime(CDRfootprint, 0)

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
    pymol.cmd.set_view(viewSet.birdsEyeView)
    pymol.cmd.scene(key="pMHCsurface", action="store")

    pMHCsurface = fileName + "/visualisation/" + id + "_pMHCsurface.png"
    
    if ray:
        rayTime(pMHCsurface, 1)
    else:
        rayTime(pMHCsurface, 0)       

    # MHC helices cos pymol is being a pain

    locs = '+'.join(str(x) for x in MHCa1[1:])
    sele = MHCa1[0] + " and resi " + locs
    pymol.cmd.select("MHCa1", sele)
    

    locs = '+'.join(str(x) for x in MHCa2[1:])
    sele = MHCa2[0] + " and resi " + locs
    pymol.cmd.select("MHCa2", sele)
    
    pymol.cmd.hide("all")
    pymol.cmd.show("cartoon", "MHCa1")
    pymol.cmd.show("cartoon", "MHCa2")
    pymol.cmd.set("cartoon_transparency", 0.5)

    ## Photo op here
    pymol.cmd.set_view(viewSet.birdsEyeView)
    pymol.cmd.scene(key="MHChelices", action="store")

    MHChelices = fileName + "/visualisation/" + id + "_MHChelices.png"

    if ray:
        rayTime(MHChelices, 1)
    else:
        rayTime(MHChelices, 0)   

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
        pymol.cmd.save(fileName + "/pdbs/" + id + "_" + name + "_COM.pdb", name + "_COM")

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
        pymol.cmd.save(fileName + "/pdbs/" + id + "_" + name + "_COM.pdb", name + "_COM")

    for COM in COMs:
        print(COM)

    ## Photo op here
    pymol.cmd.set_view(viewSet.birdsEyeView)
    pymol.cmd.scene(key="CDRCOM", action="store")

    CDRCOMimage = fileName + "/visualisation/" + id + "_CDR_centre_of_mass.png"
    
    if ray:
        rayTime(CDRCOMimage, 1)
    else:
        rayTime(CDRCOMimage, 0)

    # Save the session
    pymol.cmd.save(fileName + "/sessions/" + id + "_autoCDRloops.pse")


    # for mean in loopMeans:
    #    print mean

    # Save the loops

    for loop in TCRAlocations:
        name = loop[0]
        pymol.cmd.save(fileName + "/pdbs/" + id + "_" + name + ".pdb", name)

    for loop in TCRBlocations:
        name = loop[0]
        pymol.cmd.save(fileName + "/pdbs/" + id + "_" + name + ".pdb", name)

    # Quit pymol
    #pymol.cmd.quit()
    print('     ~  End autoCDRloops.py v0.1 BETA  ~')


