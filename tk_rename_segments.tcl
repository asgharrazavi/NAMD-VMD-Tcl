mol new  58_pop2.psf  type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile 58_pop2.pdb type pdb
set pop2 [atomselect 0 "all"]
$pop2 set segname PIP2
animate write psf pip2.psf sel [atomselect top "resname SAPI24"]
animate write pdb pip2.pdb sel [atomselect top "resname SAPI24"]

mol new  58_pops.psf  type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile 58_pops.pdb type pdb
set pops [atomselect 1 "all"]
$pops set segname POPS
animate write psf pops.psf sel [atomselect top "resname POPS"]
animate write pdb pops.pdb sel [atomselect top "resname POPS"]

mol new  58_protein.psf  type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile 58_protein.pdb type pdb
animate write psf protein.psf sel [atomselect top "protein"]
animate write pdb protein.pdb sel [atomselect top "protein"]

mol new  58_popc.psf  type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile 58_popc.pdb type pdb
set popc [atomselect 3 "all"]
$popc set segname POPC
animate write psf popc.psf sel [atomselect top "resname POPC"]
animate write pdb popc.pdb sel [atomselect top "resname POPC"]

mol new  58_pope.psf  type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile 58_pope.pdb type pdb
set pope [atomselect 4 "all"]
$pope set segname POPE
animate write psf pope.psf sel [atomselect top "resname POPE"]
animate write pdb pope.pdb sel [atomselect top "resname POPE"]

mol new  58_dpsm.psf  type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile 58_dpsm.pdb type pdb
set dpsm [atomselect 5 "all"]
$dpsm set segname SSM
animate write psf dpsm.psf sel [atomselect top "resname SSM"]
animate write pdb dpsm.pdb sel [atomselect top "resname SSM"]

mol new  58_chol.psf  type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile 58_chol.pdb type pdb
set chl [atomselect 6 "all"]
$chl set segname CHOL
animate write psf chol.psf sel [atomselect top "resname CHL1"]
animate write pdb chol.pdb sel [atomselect top "resname CHL1"]


quit
