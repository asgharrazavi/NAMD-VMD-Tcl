import os,sys

# generate input text file
text = ''
for i in range(25,32):
    text += '3kdp_eq_nod804_000%d.coor.dcd  ' %i
for i in range(10,25):
    text += '3kdp_eq_000%d.coor.dcd  ' %i



os.system('./catdcd -i ind.ind -o all-protein-pot.dcd %s' %text)

