#!/bin/bash
'''
#------------------------------------------------------------------------------------------------#
# USE: 
#  ./prepare_1mol_AA_system.sh  solute.gro  box_solvent.gro  NAME  #
#
# DESCRIPTION:
#   Generates a box with a molecule solvated in a solvent of choice. The script needs 4 inputs:
#   1. SOLUTE.gro, 
#   2. BOX_SOLVENT.gro
#   3. solvent name
#   4. solvent number of atoms
#
# BEFORE:  
#   Make sure your topology file is correctly set up. 
#
# COMMENTS:
#   -                                                                                              
#------------------------------------------------------------------------------------------------#


# Load GROMACS 2016.5 or later
#source /usr/local/gromacs-2019.5/bin/GMXRC
#source /Users/goncalojustino/miniconda3/envs/martini/bin.ARM_NEON_ASIMD/GMXRC
source /opt/gromacs_20251_cuda/bin/GMXRC

#  Check if the files have been passed to the script
solute="$1"
solvent_box="$2"
solvent_name="$3"
solvent_atoms="$4"
size_1=${#solute}
size_2=${#solvent_box}

#if [ $size_1 = 0 ] || [$size_2 = 0 ] ; then
if [ $size_1 = 0 ] || [ $size_2 = 0 ] ; then
        echo ""
        echo "Missing file of solute or solvent box. Check that. e.g.:"
        echo ""
        echo "(WATER) ./prepare_1mol_AA_system.sh  benz.gro  spc216.gro       SOL    3"
        exit
fi


gmx_20251cuda insert-molecules -ci $solute -nmol 1 -box 3.8 3.8 3.8 -rot xyz -seed 0 -o out.gro

gmx_20251cuda solvate -cp out.gro -cs $solvent_box -o initial.gro

cp system_EMPTY.top system.top
solvent_lines=$(grep $solvent_name initial.gro | wc -l)
solvent_molecules=$(expr $solvent_lines / $solvent_atoms )
echo "$solvent_name               $solvent_molecules" >> system.top

gmx_20251cuda grompp -p system.top -c initial.gro -f em_mod.mdp  -o 1-min.tpr -po 1-min.mdp 
gmx_20251cuda mdrun -v -deffnm 1-min
'''
gmx_20251cuda grompp -p system.top -c 1-min.gro   -f eq_mod.mdp  -o 2-eq.tpr  -po 2-eq.mdp 
gmx_20251cuda mdrun -v -deffnm 2-eq -nb gpu -nt 8 -pin on

gmx_20251cuda grompp -p system.top -c 2-eq.gro    -f run_mod.mdp -o 3-run.tpr -po 3-run.mdp 
gmx_20251cuda mdrun -v -deffnm 3-run -nb gpu -nt 8 -pin on

'''