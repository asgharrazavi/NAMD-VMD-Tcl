package require psfgen
topology top_all36_prot.rtf
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD

segment PRO1 {
        pdb pro1.pdb
        first NTER
        last CTER
        mutate 434 GLU
}
coordpdb pro1.pdb PRO1

segment PRO2 {
        pdb pro2.pdb
        first NTER
        last CTER
        mutate 434 GLU
}
coordpdb pro2.pdb PRO2


patch DISU PRO1:123 PRO1:132
patch GLUP PRO1:434
patch DISU PRO2:123 PRO2:132
patch GLUP PRO2:434
regenerate dihedrals angles

guesscoord
writepdb prot2.pdb
writepsf prot2.psf
