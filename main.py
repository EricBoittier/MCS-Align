from MCSAlign import MCSAlign

mcsa = MCSAlign("/data/unibas/boittier/hydra/train/62-53-3.out")
mcsa.set_smiles(8,9)
print("finding similarities")
mcsa.set_similarities()
#print(mcsa.similarities)
mcsa.find_matches(0.16)
print(mcsa.matches)
mcsa.find_MCS()
m1 = mcsa.matches[0][1]
print(mcsa.rotate_align(m1))


