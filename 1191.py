from MCSAlign import *

N_GEN = 20

mcsa = MCSAlign("/data/unibas/boittier/hydra/train/1191-95-3.out", "C=O")
mcsa.find_interaction_mol()

mcsa.set_smiles(6,10)
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
        
        with open(f"gens_1191/gen{index}.xyz", "w") as f:
            f.write(str_out)
            
    except Exception as e:
        print(e)


