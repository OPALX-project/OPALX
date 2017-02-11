# Plot of Enge function.
reset
set output "Enge-func-plot.png"
set term pngcairo
#set termoption enhanced

# Constants
z0 = 5.0 # Enge function origin (cm).
D = 2.0  # Magnet gap (cm).
delta1 = 10.0 # Fringe field region in front of z0.
delta2 = 10.0 # Fringe field region after of z0.

# Enge coefficients.
c0 = 0.478959
c1 = 1.911289
c2 = -1.185953
c3 = 1.630554
c4 = -1.082657
c5 = 0.318111

# Define Enge function.
f(x) = ((x > (z0 - delta1)) && (x < (z0 + delta2))) ? \
1.0 / (1.0 + exp(c0 \
+ c1 * (z0 - x) / D \
+ c2 * ((z0 - x) / D)**2 \
+ c3 * ((z0 - x) / D)**3 \
+ c4 * ((z0 - x) / D)**4 \
+ c5 * ((z0 - x) / D)**5)) : 0.0

# Set plot parameters.
unset key
set xrange[0.0:10.0]
set yrange[0.0:1.2]

set xlabel "z (cm)"
set ylabel "B_z (Normalized)"

# Draw points of interest.
set arrow from z0,0.0 to z0,1.2 nohead
set object 1 rect from z0+0.5,0.0 to 10.0,1.2 fc rgb "blue" fs transparent solid 0.15 noborder
set label "Magnet Interior" at 7.75,0.6 center

set label "Magnet Edge" at 7.0,0.25 center
set arrow from 7.0,0.2 to z0+0.5,0.0

set label "Enge Function\nOrigin" at 4.0,0.5 center
set arrow from 4.0,0.4 to z0,0.0

# Plot.
plot f(x) with l lt rgb "red"
