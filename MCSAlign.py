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
RDLogger.DisableLog('rdApp.*')


GDB11_PATH = "/data/unibas/boittier/GDB"


def mol_to_xyz_np(mol):
    lines = Chem.MolToMolBlock(mol).splitlines()[4:]
    lines = [x for x in lines if len(x.split()) > 5]
    xyz = [[float(x.split()[0]), float(x.split()[1]), float(x.split()[2])] for x in lines]
    atoms = [str(x.split()[3]) for x in lines]
    return np.array(xyz), atoms


def np_to_xyz(np_array, atoms):
    s = f"{len(atoms)} \n \n"
    for xyz, a in zip(np_array, atoms):
        s += f"{a} {xyz[0]} {xyz[1]} {xyz[2]}\n"
    return s

def keep_indices(xyz, indices):
    output = np.zeros((len(indices), 3))
    for i, index in enumerate(indices):
        output[i][0] = xyz[index][0]
        output[i][1] = xyz[index][1]
        output[i][2] = xyz[index][2]

    return output

class MCSAlign:
    def __init__(self, target):
        self.target_complex = target + ".complex.mol2"
        self.target_molecule = target + ".molecule.mol2"

        self.target_mol = Chem.MolFromMol2File(self.target_molecule)
        self.target_molHs = Chem.MolFromMol2File(self.target_molecule, removeHs=False)
        self.target_FP = Chem.RDKFingerprint(self.target_molHs)
        self.smiles = None

        self.matches = None
        self.matches_mols = None
        self.similarities = None
        self.smile_mols = None
        self.match_mols_Hs = None

        self.mcs_smarts = None
        self.mcs_mol = None

    def set_smiles(self, GDB_min, GDB_max):
        self.smiles = []
        for i in range(GDB_min, GDB_max + 1):
            if i < 10:
                i = f"0{i}"
            elif i > 11:
                #  Too big for GDB11
                raise Exception
            else:
                pass

            lines = open(os.path.join(GDB11_PATH, f"gdb11_size{i}.smi")).readlines()
            self.smiles.extend(lines)

    def set_target(self, target):
        self.target_complex = target + ".complex.mol2"
        self.target_molecule = target + ".molecule.mol2"

        self.target_mol = Chem.MolFromMol2File(self.target_molecule)
        self.target_molHs = Chem.MolFromMol2File(self.target_molecule, removeHs=False)
        self.target_FP = Chem.RDKFingerprint(self.target_molHs)

    def set_similarities(self):
        self.similarities = []
        self.smile_mols = []
        for i, sml in enumerate(self.smiles):
            mol = Chem.MolFromSmiles(sml.split()[0])
            mol_fp = Chem.RDKFingerprint(mol)

            # calculate similarity
            sim = DataStructs.FingerprintSimilarity(mol_fp, self.target_FP)
            # ignore identical molecules
            if sim == 1:
                sim = 0

            self.similarities.append(sim)
            self.smile_mols.append(mol)

    def find_matches(self, cutoff):
        self.matches = []
        self.matches_mols = []
        for i, sim in enumerate(self.similarities):
            if sim > cutoff:
                mol = self.smile_mols[i]
                m_Hs = Chem.AddHs(mol)
                AllChem.EmbedMolecule(m_Hs, randomSeed=0xf00d)
                mol = m_Hs
                smi = Chem.MolToSmiles(mol)
                self.matches.append((smi, mol))
                self.matches_mols.append(mol)

    def find_indices(self, mol):
        return mol.GetSubstructMatches(self.mcs_mol)

    def find_MCS(self):
        mcs_mols = self.matches_mols
        mcs_mols.append(self.target_molHs)
        res = rdFMCS.FindMCS(mcs_mols, bondCompare=rdFMCS.BondCompare.CompareOrderExact,
                             ringMatchesRingOnly=True)
        self.mcs_smarts = res.smartsString
        self.mcs_mol = Chem.MolFromSmarts(self.mcs_smarts)

    def rotate_align(self, mol):
        m1_xyz, m1_atoms = mol_to_xyz_np(self.target_molHs)
        m2_xyz, m2_atoms = mol_to_xyz_np(mol)

        com1 = np.mean(m1_xyz.T, axis=1)
        com2 = np.mean(m2_xyz.T, axis=1)
        trans = com1 - com2

        indices_match1 = self.find_indices(self.target_molHs)
        indices_match2 = self.find_indices(mol)
        align_1 = keep_indices(m1_xyz, indices_match1[0])
        align_2 = keep_indices(m2_xyz, indices_match2[0])

        rotation, rmsd = Rotation.Rotation.align_vectors(align_1, align_2+trans)

        rotated = rotation.apply(m2_xyz)
        translated = rotated + trans

        xyz_final_str = np_to_xyz(translated, m2_atoms)
        xyz_start_str = np_to_xyz(m1_xyz, m1_atoms)

        with open("test1.xyz", "w") as f:
            f.write(xyz_final_str)
        with open("test2.xyz", "w") as f:
            f.write(xyz_start_str)

        return rotated, m2_atoms





