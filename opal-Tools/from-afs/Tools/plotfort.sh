# plotfort.sh
# <number> <title> <row number> <collumn number> <xtitle> <ytitle>
#
fn=fort.$1

gnuplot <<EOF
set grid
set time
set title ' $2'
set xlabel '$5'  font "Times,24"
set ylabel '$6'  font "Times,24"
set pointsize 2.0
set terminal postscript  enhanced color "Helvetica" 20
set output "$fn.ps"
plot '$fn' u $3:$4  w l
EOF
