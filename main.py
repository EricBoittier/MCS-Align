from MCSAlign import MCSAlign

mcsa = MCSAlign("/data/unibas/boittier/hydra/train/62-53-3.out")
mcsa.set_smiles(8,9)
print("finding similarities")
mcsa.set_similarities()
#print(mcsa.similarities)
mcsa.find_matches(0.16)
print(mcsa.matches)
