from MCSAlign import *

N_GEN = 20

mcsa = MCSAlign("/data/unibas/boittier/hydra/train/67-64-1.out", "C=O")

mcsa.find_interaction_mol()

mcsa.set_smiles(5,11)
mcsa.set_similarities()

print("finding matches structures")

mcsa.find_matches(0.16, n_values=20)

mcsa.find_MCSs()

print("generating structures")

for index in range(N_GEN):
    m1 = mcsa.all_mcs_structures[index]
    
    try:

        atoms, names = (mcsa.rotate_align(index, mcsa.motif))
        
        str_out = np_to_xyz(atoms, names)
        
        with open(f"gens_67/gen{index}.xyz", "w") as f:
            f.write(str_out)
    except Exception as e:
        print(e)        

