#package require bigdcd
source bigdcd.tcl

proc pucker5 {frame} {
global mol selt1 selt2 selt3 selt4 selt5  f_r_out all outfile puckeroutfile
set size 5
set sel1 [atomselect $mol "$selt1"]
set sel2 [atomselect $mol "$selt2"]
set sel3 [atomselect $mol "$selt3"]
set sel4 [atomselect $mol "$selt4"]
set sel5 [atomselect $mol "$selt5"]

#. frames are looped over by bigdcd
 puts "frame $frame "

#. get coords from selections
  set com1 [measure center $sel1 weight mass]
  set com2 [measure center $sel2 weight mass]
  set com3 [measure center $sel3 weight mass]
  set com4 [measure center $sel4 weight mass]
  set com5 [measure center $sel5 weight mass]

  puts $outfile "$frame $com1 $com2 $com3 $com4 $com5 "

}



# global
set mol [mol new example.pdb type pdb waitfor all]
#. selections
set selt1 "resname RNGQ and name C2'"
set selt2 "resname RNGQ and name C3'"
set selt3 "resname RNGQ and name C4'"
set selt4 "resname RNGQ and name O4'"
set selt5 "resname RNGQ and name C1'"

set all [atomselect $mol all]
set f_r_out "out"
set outfile [open $f_r_out w]

set files [glob dcd/*.dcd] 
#. which pucker to call 
bigdcd pucker5 auto {*}$files
bigdcd_wait
close $outfile
quit
