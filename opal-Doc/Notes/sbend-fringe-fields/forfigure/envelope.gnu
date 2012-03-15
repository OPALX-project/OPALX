set terminal postscript landscape enhanced color colortext  "Helvetica" 20
#set term postscript eps enhanced


#set time
set grid

# lw 3 "Helvetica" 14 

#set output "path.eps"
#set xlabel "z (m)"  0.000000,0.000000  font ""
#set ylabel "x (m)"  0.000000,0.000000  font ""
#plot 'test_position.txt' u 3:1 title 'path'

#set output "| ps2pdf - envelopefit.pdf"
set output "envelopefit.eps"
#set output "envelope.eps"
#set title 'envelope y'
set xrange [0:75]

set multiplot  
set origin 0.0,0.0 
set size 1.0,0.5 
set label "*SM4" at 1.894, 0.0081  font "Helvetica, 15"
set label "*AXC" at 11.4969, 0.0051  font "Helvetica, 15"
set label "*AXD" at 43.3515, 0.005  font "Helvetica, 15"
set label "*QXA11" at 25.2009, 0.005  font "Helvetica, 15"
set label "*ZS3" at 54.0121, 0.006  font "Helvetica, 15"
set label "*MIC" at 59.9499, 0.004  font "Helvetica, 15"
#set label "KX9I" at 34.002, 0.005  font "Helvetica, 15"

set xlabel "pos (m)"  0.000000,0.000000   font ""
set ylabel "2 {/Symbol s} envelope_x (m)"  0.000000,0.000000  font ""
plot 'opalfit.stat' u ($2):6 w l lw 5  title 'opal-t' , 'measurex' using 1:2:3:4 title 'measure'  with yerrorbars lw 5, 'transfit' u 1:($3/1000) title 'transport' w l lw 5
#,'transferline.stat' u ($2):6 w l lt 5 lw 5
#plot  'wholenpx' u 1:2 title 'transport' w l lw 3,'transferline.stat' u 2:6 title 'opal-t'  w l lw 3


set origin 0.0,0.5
set size 1.0,0.5 
set xlabel "pos (m)" 
set ylabel "2 {/Symbol s} envelope_y (m)"  
plot 'opalfit.stat' u ($2):7 w l lw 5 title 'opal-t', 'measurey' using 1:2:3:4 title 'measure'  with yerrorbars lw 5, 'transfit' u 1:($5/1000) title 'transport' w l lw 5
#, 'transferline.stat' u ($2):7 w l lt 5 lw 5
#plot 'wholenpy' u 1:2 title 'transport'  w l lw 3,'transferline.stat' u 2:7 title 'opal-t'  w l  lw 3
#set label "SM4" at 1.894,.003
unset multiplot



