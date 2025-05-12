from src.MCSAlign import *

N_GEN = 20

mcsa = MCSAlign("hydra/train/2406-25-9.out", "N-O")
mcsa.find_interaction_mol()

mcsa.set_smiles(5,11)
mcsa.set_similarities()

print("finding matches structures")

mcsa.find_matches(0.16, n_values=20)

mcsa.find_MCSs()

print("generating structures")

for index in range(N_GEN):
    m1 = mcsa.all_mcs_structures[index]
    atoms, names = (mcsa.rotate_align(index, mcsa.motif))
    
    str_out = np_to_xyz(atoms, names)
    
    with open(f"gens_2406/gen{index}.xyz", "w") as f:
        f.write(str_out)
        

