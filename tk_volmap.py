import os,sys,glob
import numpy as np
import multiprocessing

def job(idd):
    os.system('mkdir %d' %idd)
    os.chdir('%s' %idd)
    for j in ['A','B','C']:
        os.system(''' echo 10 11 |  gmx trjconv -n index.ndx -s protein_A_com_removed.pdb -f protein_200ps_%s_%d.xtc -fit rot+trans  -o t_.pdb -sep  -skip 5 ''' %(j,idd))
  	n_pdb = glob.glob('t_*pdb')
  	for k in range(len(n_pdb)):
	    os.system(''' /Applications/VMD_1.9.2.app/Contents/vmd/vmd_MACOSXX86 -dispdev none  t_%d.pdb <<EOF 
volmap density [atomselect 0 "protein and resid 350 to 358"] -res 1.0 -weight none -allframes -combine avg -o  hp2_%s_%d.dx -minmax {{-30.0 -30.0 -20.0} {30.0 30.0 20.0}} -radscale 1.0 -checkpoint 0
volmap density [atomselect 0 "protein and resid 50 to 70   140 to 155"] -res 1.0 -weight none -allframes -combine avg -o  tri_%s_%d.dx -minmax {{-30.0 -30.0 -20.0} {30.0 30.0 20.0}} -radscale 1.0 -checkpoint 0

quit
EOF
''' %(k,j,k,j,k))
        os.system('rm t*pdb')
    os.chdir('../')

for i in range(2,13):
    p = multiprocessing.Process(target=job,args=([i]))
    p.start()


