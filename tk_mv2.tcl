# script for calculating velocities (MV**2) of atoms based on *.vel file
# load *.vel in frame 1 as moladdfile *.vel type namdbin

set aa [atomselect top "protein and resid 1 to 620 and name CA" frame 1]
set v2 [vecadd [vecmul [$aa get x] [$aa get x] ] [vecmul [$aa get y] [$aa get y]] [vecmul [$aa get z] [$aa get z]] ]
set mv2 [vecmul [$aa get mass] [exp {$v2}]]
$aa set beta $mv2
set outfile [open "mv2_beta.txt" "w"]
puts $outfile $mv2
