#!/bin/bash
# pipeline.sh
# Runs PyMOL mutagenesis then a full GROMACS MD pipeline.
#
# Usage:
#   bash pipeline.sh

set -euo pipefail

# ── 1. PyMOL mutagenesis ──────────────────────────────────────────────────────
echo "==> Applying mutations with PyMOL..."
pymol -cq mutate.py

MUTANT_PDB=$(cat mutant_pdb.txt)
echo "==> Using mutant PDB: $MUTANT_PDB"

# ── 2. pdb2gmx ────────────────────────────────────────────────────────────────
# Selections: 1 = CHARMM-36m force field, 1 = TIP3P water model
echo "==> Running pdb2gmx..."
printf "1\n1\n" | gmx pdb2gmx -f "$MUTANT_PDB" -o brazzein_processed.gro -ignh

# ── 3. Box & solvation ────────────────────────────────────────────────────────
echo "==> Setting up simulation box..."
gmx editconf -f brazzein_processed.gro -o brazzein_newbox.gro -bt dodecahedron -d 1.5 -c

echo "==> Solvating..."
gmx solvate -cp brazzein_newbox.gro -cs spc216.gro -o brazzein_solv.gro -p topol.top

# ── 4. Ions ───────────────────────────────────────────────────────────────────
echo "==> Adding ions..."
gmx grompp -f ions.mdp -c brazzein_solv.gro -p topol.top -o ions.tpr
# Selection 13 = SOL (solvent group) for ion replacement
printf "13\n" | gmx genion -s ions.tpr -o brazzein_solv_ions.gro \
    -p topol.top -pname NA -nname CL -neutral -conc 0.15

# ── 5. Energy minimisation ────────────────────────────────────────────────────
echo "==> Running energy minimisation..."
gmx grompp -f em.mdp -c brazzein_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# ── 6. NVT equilibration ──────────────────────────────────────────────────────
echo "==> Running NVT equilibration..."
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

# ── 7. NPT equilibration ──────────────────────────────────────────────────────
echo "==> Running NPT equilibration..."
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt

# ── 8. Index file ─────────────────────────────────────────────────────────────
echo "==> Building index file..."
python find_number.py

{
    read receptor_start
    read receptor_end
    read ligand_start
    read ligand_end
} < atom_number.txt

gmx make_ndx -f npt.gro -o index.ndx <<EOF
a ${receptor_start}-${receptor_end}
name 17 receptor
a ${ligand_start}-${ligand_end}
name 18 ligand
q
EOF

# ── 9. Production MD ──────────────────────────────────────────────────────────
echo "==> Running production MD..."
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr
gmx mdrun -v -deffnm md

# ── 10. Trajectory processing ─────────────────────────────────────────────────
echo "==> Processing trajectory..."
printf "1\n1\n0\n" | gmx trjconv -s md.tpr -f md.xtc \
    -o md_cluster.xtc -pbc cluster -center -n index.ndx

printf "17\n1\n" | gmx trjconv -s md.tpr -f md_cluster.xtc \
    -fit rot+trans -o md_trajectory.xtc -n index.ndx

printf "1\n" | gmx trjconv -s md.tpr -f md_trajectory.xtc \
    -o last_coordinate.pdb -dump 10000

# ── 11. Interaction energy ────────────────────────────────────────────────────
# Energy terms 51 & 52: Coul-SR:receptor-ligand and LJ-SR:receptor-ligand
# Verify indices match your topology: gmx energy -f md.edr (run interactively)
echo "==> Extracting interaction energy..."
printf "51\n52\n" | gmx energy -f md.edr -o interaction_energy.xvg -xvg none

echo ""
echo "✓ Pipeline complete."
