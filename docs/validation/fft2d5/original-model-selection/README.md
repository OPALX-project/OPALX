# Explicit original FFT2D5 selection

Use `TYPE=FFT2D5, PIPEMODE=OPEN, FFT2D5MODEL="ORIGINAL"` to identify the
independent transverse 2D slices and local longitudinal closure (g=6.36 for OPEN).
ORIGINAL is the default and only supported value in this standalone change.
Omitting it preserves existing behavior. The value is printed in the solver
summary; specifying it for another solver type is rejected. No field operator,
pusher, frame handling or physics tolerance is changed by this PR.

## Reference benchmark

The attached reference comparison uses the local benchmark build with ORIGINAL
selected (SHA a2c87276106782554324753c889db203629f1b7f645d35e2f1bd30b00e32fa9e),
clean committed658 fd6ba8dc5 (SHA ab1febaafc26f4dc82904808dfb81b6de7c7476d42290fc7cc79a7090e0602f3),
and native PyORBIT. These are reference data, not PR-head tracking results.
The local benchmark also contains earlier numerical changes which are outside
this isolated PR. Curvature and3D operators were disabled in that benchmark.

R=2m,10m,100steps,5m uniform bunch,10mm transverse RMS,10pC orZeroSC,
beta=.001 protons,1,048,576 common particles,512^3 charged grid with full
widths256/320mm. Source SHA4cf23b08e775af54ab1b4d0cd25d090793d36d0735c6c0bd332ddc1af3939c59.

Data: series.csv (606 records) and endpoint comparison in results.csv.
Plots show RMS sizes and normalized local transverse emittances; longitudinal
emittance is the s–energy RMS area. Only the upperright longitudinal emittance
panel is logarithmic. The ZeroSC radial panel shows separate common-time and common-station analytic
references, using physical-volume disk weights and 256-node quadrature.

The local ORIGINAL and committed 658 results agree closely. At 10 pC, final
radial and vertical RMS sizes agree with PyORBIT within 0.6%; radial and vertical
emittances differ by 12% and 50%. The maximum radial-size difference along the
trajectory is 11%, so endpoint agreement does not establish convergence.

## Analytic ZeroSC reference

Let a be the initial radial offset, R the bend radius, phi=s/R, and p0=beta*gamma.
At common time, r=sqrt(R²+2Ra cos(phi)+a²) and q_u=-p0*a*sin(phi)/r.
At common station, r=a*cos(phi)+sqrt(R²-a²*sin(phi)²) and
q_u=-p0*a*sin(phi)/R. In both cases u=r-R. The mean-subtracted covariance
of (u,q_u) gives the normalized RMS emittance. Moments use the radial density
2*sqrt(b²-a²)/(pi*b²)*(1+a/R), with disk radius b=20 mm.

OPALX samples at common time; PyORBIT TEAPOT samples at common lattice station.
Using the matching reference, the 100-step endpoint emittance errors are 0.028%
for OPALX and 0.025% for PyORBIT. PyORBIT's maximum absolute curve error is
0.59% of the analytic peak. The earlier apparent 15% endpoint discrepancy arose
from comparing these different sampling conventions. This check concerns ZeroSC;
charged trajectories require separate care when comparing sampling conventions.

Recreate the figures with `python plot_results.py` (NumPy, pandas, Matplotlib).
