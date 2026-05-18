"""
mutate.py
Applies point mutations to a PDB file using PyMOL and saves the mutant.

Usage:
    pymol -cq mutate.py
"""

import os
import sys
from pymol import cmd

# ── Edit these values ──────────────────────────────────────────────────────────
PDB_FILE = "9NOR_4HE7.pdb"

# Add all your mutations here as (CHAIN, RES_NUM, ORIG_AA, MUT_AA)
MUTATIONS = [
    ("C", 31, "H", "R"),
    ("C", 36, "E", "D"),
    ("C", 41, "E", "A"),
]

# ──────────────────────────────────────────────────────────────────────────────

MUT_AA3 = {
    "A": "ALA", "G": "GLY", "V": "VAL", "L": "LEU", "I": "ILE",
    "P": "PRO", "F": "PHE", "W": "TRP", "M": "MET", "S": "SER",
    "T": "THR", "C": "CYS", "Y": "TYR", "H": "HIS", "D": "ASP",
    "E": "GLU", "N": "ASN", "Q": "GLN", "K": "LYS", "R": "ARG",
}


mutation_id = "_".join(f"{orig}{res}{mut}" for (chain, res, orig, mut) in MUTATIONS)
mutant_pdb  = f"{mutation_id}.pdb"

cmd.load(PDB_FILE, "complex")
print(f"Chains: {cmd.get_chains('complex')}")
print(f"Atoms:  {cmd.count_atoms('complex')}")

for (chain, res_num, orig_aa, mut_aa) in MUTATIONS:
    target_sel = f"chain {chain} and resi {res_num}"

    if cmd.count_atoms(target_sel) == 0:
        sys.exit(f"ERROR: Residue not found — chain {chain}, resi {res_num}. Aborting.")

    print(f"Applying {orig_aa}{res_num}{mut_aa} on chain {chain}...")
    cmd.wizard("mutagenesis")
    cmd.get_wizard().do_select(target_sel)
    cmd.get_wizard().set_mode(MUT_AA3[mut_aa])
    cmd.get_wizard().apply()
    cmd.set_wizard()
    print(f"  ✓ {orig_aa}{res_num}{mut_aa} applied")

cmd.save(mutant_pdb, "complex")
print(f"\n✓ Mutant PDB saved: {mutant_pdb}")

# Write the mutant PDB path to a file so pipeline.sh can pick it up
with open("mutant_pdb.txt", "w") as f:
    f.write(mutant_pdb)
