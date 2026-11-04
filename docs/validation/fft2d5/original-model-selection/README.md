# Explicit original FFT2D5 selection

Use `TYPE=FFT2D5, PIPEMODE=OPEN, FFT2D5MODEL="ORIGINAL"` to identify the
independent transverse 2D slices and local longitudinal closure(g=6.36 forOPEN).
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

R=2m,10m,10steps,5m uniform bunch,10mm transverse RMS,10pC orZeroSC,
beta=.001 protons,1,048,576 common particles,512^3 charged grid with full
widths256/320mm. Source SHA4cf23b08e775af54ab1b4d0cd25d090793d36d0735c6c0bd332ddc1af3939c59.

Data: series.csv (all66 endpoint records) and endpoint comparison in results.csv.
Plots show RMS sizes and normalized local transverse emittances; longitudinal
emittance is the s–energy RMS area. Only the upperright longitudinal emittance
panel is logarithmic. The black dotted curve is the analytic uniformB radial
ZeroSC reference, using physical-volume disk weights and256-node quadrature.

The coarse10step Boris orbit phase differs from PyORBIT; agreement in final
RMS sizes does not establish tracking convergence. Charged emittances remain
discrepant. The local original model agrees closely with committed658 for
these data. The failed curvature-model attempt is outside this comparison.

Recreate the figures with `python plot_results.py` (NumPy,pandas,Matplotlib).
