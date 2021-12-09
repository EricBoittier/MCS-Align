from MCSAlign import MCSAlign

mcsa = MCSAlign("/home/boittier/hydra/train/62-53-3.out")
mcsa.set_smiles(9,9)
mcsa.set_similarities()
mcsa.find_matches(0.1)

