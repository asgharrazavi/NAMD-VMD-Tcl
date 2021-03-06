set sel1 [atomselect top "protein and resid 60 and name CZ"]
set sel2 [atomselect top "protein and resid 446 and name CD"]
set nf [molinfo top get numframes]
set outfile [open "r60_e446.txt" "w"]
for {set i 0} {$i < $nf} {incr i} {
    $sel1 frame $i
    $sel2 frame $i
    set com1 [measure center $sel1 weight mass]
    set com2 [measure center $sel2 weight mass]
    set simdata($i.r) [veclength [vecsub $com1 $com2]]
    puts $outfile "$simdata($i.r)"
}

close $outfile

