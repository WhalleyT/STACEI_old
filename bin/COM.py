import pymol
import numpy as np

def init_pymol():
    print("\nInitialising pymol...\n")
    pymol.finish_launching(['pymol', '-nqc'])
    pymol.cmd.reinitialize()

    return np.array([float(x), float(y), float(z)])


def dist(residue_one, residue_two):
    diff = residue_one - residue_two
    return np.sqrt(np.sum(diff * diff))


def COM(pdb, pairs_of_chains):
    pymol.cmd.delete("all")
    print("\nAligning file to template...\n")
    pymol.cmd.load(pdb)

    distances = []
    for pair in pairs_of_chains:
        print(pair[0])

        pymol.cmd.select("to", "chain %s" % pair[0])
        pymol.cmd.select("from", "chain %s" % pair[1])

        to_COM = np.asarray(pymol.cmd.centerofmass("to"), dtype=np.float32)
        from_COM = np.asarray(pymol.cmd.centerofmass("from"), dtype=np.float32)

        distance = dist(to_COM, from_COM)

        distances.append(distance)

    return distances


COM("pdbs/1zgl.pdb", [("A", "B")])