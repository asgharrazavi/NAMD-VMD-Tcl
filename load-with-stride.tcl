
# add a molecule
mol new  NKA_4hqj.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

# add first 10
for {set i 1} {$i < 9} {incr i} {
    mol addfile NKA_4hqj_0000$i.coor.dcd type dcd first 0 last -1 step 10 filebonds 1 autobonds 1 waitfor all
}

# add first 10
for {set i 10} {$i < 60} {incr i} {
    mol addfile NKA_4hqj_000$i.coor.dcd type dcd first 0 last -1 step 10 filebonds 1 autobonds 1 waitfor all
}
