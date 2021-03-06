###########################
# script for mutation residues on GltPh trimer
###########################


# load a molecule
mol new start.psf
mol addfile start.pdb

# parse segments
animate write psf non_protein.psf sel [atomselect top "not protein"]
animate write pdb non_protein.pdb sel [atomselect top "not protein"]

animate write pdb prot_A.pdb sel [atomselect top "chain A"]
animate write pdb prot_B.pdb sel [atomselect top "chain B"]
animate write pdb prot_C.pdb sel [atomselect top "chain C"]

# mutate protein segments
package require psfgen

topology top_all36_prot.rtf
topology top_all36_lipid.rtf
topology toppar_water_ions.str

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
pdbalias residue NA SOD
pdbalias atom NA NA SOD
pdbalias residue HOH TIP3
pdbalias atom TIP3 OW OH2

resetpsf

foreach S { A B C} {
   segment $S {
      pdb prot_$S.pdb
      mutate 345 VAL
      mutate 366 ALA
   }
   coordpdb prot_$S.pdb $S
}

guesscoord
regenerate angles dihedrals 

writepdb protein.pdb
writepsf protein.psf


# combine protein and non-protein segments
resetpsf

readpsf  protein.psf
coordpdb protein.pdb

readpsf  non_protein.psf
coordpdb non_protein.pdb

writepdb mol.pdb
writepsf mol.psf

quit
