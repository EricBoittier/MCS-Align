import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import IPythonConsole
import py3Dmol
import scipy.spatial.transform as Rotation
import numpy as np
import os

IPythonConsole.drawOptions.addAtomIndices = True
#RDLogger.DisableLog('rdApp.*')


GDB11_PATH = "/data/unibas/boittier/GDB"

global_smiles = []
global_mols = []

def set_global_smiles(GDB_min, GDB_max):
    global_smiles = []
    for i in range(GDB_min, GDB_max + 1):
        if i < 10:
            i = f"0{i}"
        elif i > 11:
            #  Too big for GDB11
            raise Exception
        else:
            pass

        lines = open(os.path.join(GDB11_PATH, f"gdb11_size{i}.smi")).readlines()
        global_smiles.extend(lines)
    return global_smiles
        
def set_global_mols(global_smiles):
    global_mols = []
    for sml in global_smiles:
        mol = Chem.MolFromSmiles(sml.split()[0])
        mol = Chem.AddHs(mol)
        global_mols.append(mol)
    return global_mols


def mol_to_xyz_np(mol):
    lines = Chem.MolToMolBlock(mol).splitlines()[4:]
    lines = [x for x in lines if len(x.split()) > 15]
    xyz = [[float(x.split()[0]), float(x.split()[1]), float(x.split()[2])] for x in lines]
    atoms = [str(x.split()[3]) for x in lines]
    return np.array(xyz), atoms


def np_to_xyz(np_array, atoms):
    s = f"{len(atoms)} \n \n"
    for xyz, a in zip(np_array, atoms):
        s += f"{a} {xyz[0]} {xyz[1]} {xyz[2]}\n"
    return s

def keep_indices(xyz, indices, atom_names):
    
    _ = [index for index in indices if not atom_names[index] == "H"]

    output = [None] * len(_)
    
    for i, index in enumerate(_):
        output[i] = xyz[index]

    return np.array(output)

class MCSAlign:
    def __init__(self, target, motif, using_pdb=False):
        self.target_complex = target + ".complex.mol2"
        self.target_molecule = target + ".molecule.mol2"
        
        self.motif = Chem.MolFromSmarts(motif)

        if not using_pdb:
        
            self.target_complex_mol = Chem.MolFromMol2File(self.target_complex)

            if self.target_complex_mol is not None:
    #             print("self.target_complex_mol is None")
                self.target_complex_molHs = Chem.MolFromMol2File(self.target_complex, removeHs=False)

            self.target_mol = Chem.MolFromMol2File(self.target_molecule)
            
            if self.target_mol is not None:
#             print("self.target_mol is None")
                self.target_molHs = Chem.MolFromMol2File(self.target_molecule, removeHs=False)
            
        else:
            self.target_complex += ".pdb"
            self.target_molecule += ".pdb"
            
            self.target_complex_mol = Chem.MolFromPDBFile(self.target_complex)
            print(self.target_complex_mol)

            if self.target_complex_mol is not None:
    #             print("self.target_complex_mol is None")
                self.target_complex_molHs = Chem.MolFromPDBFile(self.target_complex, removeHs=False)

            self.target_mol = Chem.MolFromPDBFile(self.target_molecule)
            if self.target_mol is not None:
#             print("self.target_mol is None")
                self.target_molHs = Chem.MolFromPDBFile(self.target_molecule, removeHs=False)
        

        if self.target_molHs is not None:
            self.target_FP = Chem.RDKFingerprint(self.target_molHs)
            
        self.smiles = None
        self.mols = None
        
        self.matches = None
        self.matches_mols = None
        self.similarities = None
        self.smile_mols = None
        self.match_mols_Hs = None

        self.mcs_smarts = None
        self.mcs_mol = None
        
        self.all_mcs_smarts = None
        self.all_mcs_mols = None
        self.all_mcs_structures = None
        
        self.interaction_mol_xyz = None
        self.interaction_mol_names = None
        
    def find_interaction_mol(self):
        complex_xyz, complex_names = mol_to_xyz_np(self.target_complex_molHs)
        target_xyz, target_names = mol_to_xyz_np(self.target_molHs)
        self.interaction_mol_xyz = []
        self.interaction_mol_names = []
        for xyz, name in zip(complex_xyz, complex_names):
            if xyz not in target_xyz:
                self.interaction_mol_names.append(name)
                self.interaction_mol_xyz.append(xyz)
            pass
        

    def set_smiles(self, smiles):
        self.smiles = smiles
        
    def set_mols(self, mols):
        self.mols = mols


    def set_target(self, target):
        self.target_complex = target + ".complex.mol2"
        self.target_molecule = target + ".molecule.mol2"

        self.target_mol = Chem.MolFromMol2File(self.target_molecule)
        self.target_molHs = Chem.MolFromMol2File(self.target_molecule, removeHs=False)
        self.target_FP = Chem.RDKFingerprint(self.target_molHs)

    def set_similarities(self):
        self.similarities = []
        self.smile_mols = []
        for i, mol in enumerate(self.mols):
            # Check for match
            if mol.HasSubstructMatch(self.motif):
                mol_fp = Chem.RDKFingerprint(mol)
                # calculate similarity
                sim = DataStructs.FingerprintSimilarity(mol_fp, self.target_FP)
                # ignore identical molecules
                if sim == 1:
                    sim = 0

                self.similarities.append(sim)
                self.smile_mols.append(mol)

    def find_matches(self, cutoff, n_values=False):
        self.matches = []
        self.matches_mols = []
        
        if n_values:
            keys = sorted(range(len(self.similarities)), key=lambda i: self.similarities[i])[-n_values:]
            for i, key in enumerate(keys):
                mol = self.smile_mols[key]
                AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
                smi = Chem.MolToSmiles(mol)
                self.matches.append((smi, mol))
                self.matches_mols.append(mol)
        else:
            for i, sim in enumerate(self.similarities):
                if sim > cutoff:
                    mol = self.smile_mols[i]
                    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
                    smi = Chem.MolToSmiles(mol)
                    self.matches.append((smi, mol))
                    self.matches_mols.append(mol)

    def find_indices(self, mol, mcs):
        return mol.GetSubstructMatches(mcs)

    def find_MCS(self):
        mcs_mols = self.matches_mols
        mcs_mols.append(self.target_molHs)
                
        res = rdFMCS.FindMCS(mcs_mols)  #, bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                                        #ringMatchesRingOnly=True)
        self.mcs_smarts = res.smartsString
        self.mcs_mol = Chem.MolFromSmarts(self.mcs_smarts)
        
    def find_MCSs(self):
        self.all_mcs_smarts = []
        self.all_mcs_mols = []
        self.all_mcs_structures = []
                
        for mol in self.matches_mols:
            _ = [mol]
            _.append(self.target_molHs)
                
            res = rdFMCS.FindMCS(_)  #, bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                                        #ringMatchesRingOnly=True)
            self.all_mcs_smarts.append(res.smartsString)
            self.all_mcs_mols.append(Chem.MolFromSmarts(res.smartsString))
            self.all_mcs_structures.append(mol)

    def rotate_align(self, match_index, mcs):
        #  Load match data
        mol = self.all_mcs_structures[match_index]
        mol_mcs = self.all_mcs_mols[match_index]
        m1_xyz, m1_atoms = mol_to_xyz_np(self.target_molHs)
        m2_xyz, m2_atoms = mol_to_xyz_np(mol)

        #  Match indices from complete MCS
        indices_match1 = self.find_indices(self.target_molHs, mol_mcs)
        indices_match2 = self.find_indices(mol, mol_mcs)
        align_1 = keep_indices(m1_xyz, indices_match1[0], m1_atoms)
        align_2 = keep_indices(m2_xyz, indices_match2[0], m2_atoms)
        #  Perform rotation
        rotation, rmsd = Rotation.Rotation.align_vectors(align_1, align_2)
        rotated = rotation.apply(m2_xyz)
        #  Perform translation
        com1 = np.mean(align_1.T, axis=1)
        align_2 = keep_indices(rotated, indices_match2[0], m2_atoms)
        com2 = np.mean(align_2.T, axis=1)
        trans = com1 - com2
        rotated = rotated + trans
        
        #  Match indices from the substructure given
        indices_match1 = self.find_indices(self.target_molHs, mcs)
        indices_match2 = self.find_indices(mol, mcs)
        align_1 = keep_indices(m1_xyz, indices_match1[0], m1_atoms)
        align_2 = keep_indices(rotated, indices_match2[0], m2_atoms)
        #  Perform rotation
        rotation, rmsd = Rotation.Rotation.align_vectors(align_1, align_2)
        rotated = rotation.apply(rotated)
        #  Perform translation
        com1 = np.mean(align_1.T, axis=1)
        align_2 = keep_indices(rotated, indices_match2[0], m2_atoms)
        com2 = np.mean(align_2.T, axis=1)
        trans = com1 - com2
        rotated = rotated + trans
        
        
        rotated = np.append(rotated, self.interaction_mol_xyz, axis=0)
        m2_atoms.extend(self.interaction_mol_names)
        
        #  save output for debugging
        xyz_start_str = np_to_xyz(m1_xyz, m1_atoms)
#         xyz_final_str = np_to_xyz(m2_xyz_trans, m2_atoms)
        xyz_final1_str = np_to_xyz(rotated, m2_atoms)
        
        with open("test1.xyz", "w") as f:
            f.write(xyz_start_str)
#         with open("test2.xyz", "w") as f:
#             f.write(xyz_final_str)
        with open("test3.xyz", "w") as f:
            f.write(xyz_final1_str)
            
        return rotated, m2_atoms





