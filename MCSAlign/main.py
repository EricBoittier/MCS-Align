from src.MCSAlign import MCSAlign, np_to_xyz, set_global_smiles, set_global_mols

N_GEN = 20
TANIMOTO_CUTOFF = 0.1

global_smiles = set_global_smiles(5,10)
global_mols = set_global_mols(global_smiles)

def generate(mcsa, name):
    mcsa.find_interaction_mol()
    mcsa.set_smiles(global_smiles)
    mcsa.set_mols(global_mols)
    mcsa.set_similarities()
    mcsa.find_matches(TANIMOTO_CUTOFF, n_values=N_GEN)
    mcsa.find_MCSs()
    print("generating structures")
    
    for index in range(N_GEN):
        try:
            m1 = mcsa.all_mcs_structures[index]
            atoms, names = (mcsa.rotate_align(index, mcsa.motif))

            str_out = np_to_xyz(atoms, names)

            with open(f"gens_{name}/gen{index}.xyz", "w") as f:
                f.write(str_out)
                
        except Exception as e:
            print(e)
        
# name = "62"
# mcsa = MCSAlign("/data/unibas/boittier/hydra/train/62-53-3.out", "cN")
# generate(mcsa, name)

# name = "132"
# mcsa = MCSAlign("/data/unibas/boittier/hydra/train/132-64-9.out", "cOc")
# generate(mcsa, name)

# name = "327"
# mcsa = MCSAlign("/data/unibas/boittier/hydra/train/327-54-8.out", "cF")
# generate(mcsa, name)

# name = "67"
# mcsa = MCSAlign("/data/unibas/boittier/hydra/train/67-64-1.out", "C=O")
# generate(mcsa, name)

# name = "98-86"
# mcsa = MCSAlign("/data/unibas/boittier/hydra/train/98-86-2.out", "C=O")
# generate(mcsa, name)

#  New start here

# name = "109"
# mcsa = MCSAlign("/data/unibas/boittier/hydra/test/109-99-9.out", "O")
# generate(mcsa, name)

# name = "434"
# mcsa = MCSAlign("/data/unibas/boittier/hydra/test/434-45-7.out", "cC=O")
# generate(mcsa, name)

# name = "50"
# mcsa = MCSAlign("/data/unibas/boittier/hydra/test/50-0-0.out", "C=O")
# generate(mcsa, name)

# name = "502"
# mcsa = MCSAlign("/data/unibas/boittier/hydra/test/502-49-8.out", "C=O")
# generate(mcsa, name)

# name = "547"
# mcsa = MCSAlign("/data/unibas/boittier/hydra/test/547-64-8.out", "O=CCO")
# generate(mcsa, name)

# name = "75"
# mcsa = MCSAlign("/data/unibas/boittier/hydra/test/75-89-8.out", "CO")
# generate(mcsa, name)

# name = "80"
# mcsa = MCSAlign("/data/unibas/boittier/hydra/test/80-73-9.out", "C=O")
# generate(mcsa, name)


name = "98-85"
mcsa = MCSAlign("/data/unibas/boittier/hydra/train/98-85-1.out", "CO")
generate(mcsa, name)

name = "110"
print(name)
mcsa = MCSAlign("/data/unibas/boittier/hydra/test/110-86-1.out", "na")
generate(mcsa, name)

name = "125132"
print(name)
mcsa = MCSAlign("/data/unibas/boittier/hydra/test/125132-75-4.out", "OC")
generate(mcsa, name)


name = "288"
print(name)
mcsa = MCSAlign("/data/unibas/boittier/hydra/train/288-32-4.out", "NC", using_pdb=True)
generate(mcsa, name)
