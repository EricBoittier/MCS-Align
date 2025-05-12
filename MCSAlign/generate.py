from src.MCSAlign import *

N_GEN = 20

def generate(outfile, motif, smiles, mols):
    mcsa = MCSAlign(outfile, motif)

    mcsa.find_interaction_mol()

    mcsa.set_smiles(smiles)
    mcsa.set_mols(mols)
    
    mcsa.set_similarities()

    print("finding matches structures")

    mcsa.find_matches(0.16, n_values=20)

    mcsa.find_MCSs()

    print("-"*100)
    print(f"Generating {N_GEN} structures")

    for index in range(N_GEN):
        m1 = mcsa.all_mcs_structures[index]

        try:
            atoms, names = (mcsa.rotate_align(index, mcsa.motif))

            str_out = np_to_xyz(atoms, names)

            with open(f"gens_67/gen{index}.xyz", "w") as f:
                f.write(str_out)
            print(f"generated {index} of {N_GEN}")
        except Exception as e:
            print(e)        
