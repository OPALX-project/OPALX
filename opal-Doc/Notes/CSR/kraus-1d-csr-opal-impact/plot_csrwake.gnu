set line style 1 lt 1 lw 2
set line style 2 lt 2 lw 2 
set line style 3 lt 3 lw 2 
set line style 4 lt 6 lw 2 
set line style 5 lt 1 lw 2 
set line style 6 lt 2 lw 2 
set line style 7 lt 3 lw 2 
set line style 8 lt 6 lw 2 
plot 'BC1/OPAL/CSRWake96.txt' u 1:2 w l ls 1 t 'OPAL'
replot 'BC1/ImpactT/CSRPote001.txt' u ($1+2.89819e-5):2 w l ls 2 t 'Impact-T'
replot 'BC1/OPAL/CSRWake96.txt' u 1:4 w l axis x1y2 ls 3 t 'OPAL'
replot 'BC1/ImpactT/CSRPote001.txt' u ($1+2.89819e-5):3 w li axis x1y2 ls 4 t 'Impact-T'
set xlabel 's [m]' font 'Helvetica-Bold,25'
set ylabel 'E_z [V/m]' font 'Helvetica-Bold,25'
set y2label 'line density [C/m]' font 'Helvetica-Bold,25'
set ytics nomirror
set y2tics
set key 0.0032,-5000
set term post enh col 'Helvetica,22' lw 2
set out 'comparison.ps'
replot
