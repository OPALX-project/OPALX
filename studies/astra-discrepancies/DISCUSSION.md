# What the EuXFEL agreement establishes

This discussion accompanies [issue #636](https://github.com/OPALX-project/OPALX/issues/636).
It summarizes completed **frozen-executable studies**, before the current branch
integration and stationary-limit robustness amendment. No new tracking was
performed for this evidence commit. The current branch's algebraically equivalent
robust `hypot` gamma and new ±15 / 10⁻¹²-perturbed stationary tests must not be
read as the source of these previously measured trajectories.

**The original OPALX spatial method has not been shown to be physically wrong,
and the ASTRA-like research operator has not been shown necessary for an
accurate calculation.** Finite-mesh compatibility differs from physical
accuracy, which requires independent error estimates and absolute convergence.
The accepted CM correction addresses a separate kinematic inconsistency.

At the fixed 10.37303974716 ns observation, the validated 1 ps native-scalar
comparison gives the following differences for the common EuXFEL gun-exit
checkpoint, 100,000 particles, −250 pC, one energy bin, 32³ grid and 1 ps step:

| OPALX variant | εnx − ASTRA | εny − ASTRA | εnz − ASTRA |
|---|---:|---:|---:|
| Original | −4.8334% | −4.8825% | −0.4838% |
| Consistent CM, original spatial operator | −4.8439% | −4.8930% | −0.4478% |
| CM, NGP/spectral, phase 0 | −7.9706% | −7.9457% | +0.5355% |
| CM, CIC/corner, phase 0 | +2.7281% | +2.6449% | −1.1908% |
| CM, NGP/corner, phase 0 | +0.4869% | +0.4112% | −0.5088% |
| CM, NGP/corner, phase ½ | +0.0374% | +0.0380% | +0.0285% |

These use full-population, centred mechanical normalized RMS emittance with
native species masses; common-mass differences are about 0.0008 percentage
points. The [STAT audit](PROVENANCE.md#validation-scopes)
and current-case particle checks validate dense scalar sampling and bound
reference-axis effects. Avoiding sparse covariance interpolation changes the
original endpoint transverse gap by only about 0.005 percentage points.

**NGP/corner at phase 0 already meets 1% in all three emittances at all 25
clocks:** maxima are 0.5279%,0.4867%,0.5088% in xyz. Thus the half-cell source
shift is unnecessary for that criterion. This does not mean the original
CIC/spectral calculation meets it, or that either operator is converged. The
[complete native comparison](README.md#results)
preserves all six controls and nine observables.

[Complete plots, units and case legend](README.md) accompany the copied
[native-observable](tables/observables.csv), [historical-reference](tables/historical-references.csv)
and [flow](tables/flow.csv) tables. Some supporting frozen-state audit results
are summarized and fingerprinted in [PROVENANCE.md](PROVENANCE.md), while their
large raw arrays are outside this portable archive.

## Same Poisson physics, different discrete forces

Both calculations use isolated free-space electrostatics in an approximate
bulk rest frame, an integrated rectangular-cell Coulomb kernel, and laboratory
electric/magnetic reconstruction. Neither comparison activates an image plane
or conducting-cathode boundary. The complete numerical chain is nevertheless

$$
\text{particles}\longrightarrow Q_{\rm grid}\longrightarrow\phi_{\rm grid}
\longrightarrow\mathbf E_{\rm grid}\longrightarrow\mathbf F_{\rm particles}.
$$

Matching models should agree within demonstrated numerical uncertainties.
Persistent disagreement after such control would indicate an implementation
or model difference; this benchmark has not reached that condition.

Original OPALX distributes particle charge among eight samples using CIC,
differentiates spectrally, and gathers from collocated field samples. The
research operator assigns charge to one containing cell, uses the same
integrated Green function for potential, differences eight potentials to form
corner fields, and gathers at those corner locations. The complete phase-½ variant additionally shifts the source-lattice phase
by half a cell; the phase-0 combination does not. Each operation changes the effective
finite-grid force. Matching the Poisson kernel does not match that complete
operator.

“Cell-centred versus node-centred” is insufficient: the original mesh origin
is paired with IPPL's indexing convention. No missing half-cell converter
correction was identified. Both arrangements can be located consistently.
Likewise, the debated STANDARD origin override is inactive: these runs select
INTEGRATED, including a finite self-cell potential. The
[Green-function audit](PROVENANCE.md#supporting-frozen-state-studies)
establishes the common analytic kernel type, not pointwise equality of native
ASTRA kernel arrays.

For a smooth Fourier mode, let $a_j=k_jh_j/2$. The corner derivative,
referenced to its actual sampling point, is

$$
D_i=-\frac{2i}{h_i}\sin a_i\prod_{j\ne i}\cos a_j.
$$

Combining principal-mode NGP/CIC assignment factors and their shared CIC
gather gives, with common integrated kernel,

$$
\frac{T_i^{\rm NGP,corner}}{T_i^{\rm CIC,spectral}}
\simeq\prod_{j\ne i}\frac{a_j}{\tan a_j}
=1-\frac1{12}\sum_{j\ne i}k_j^2h_j^2+O((kh)^4).
$$

This demonstrates directional finite-resolution filtering and a shared
long-wavelength limit. It assumes smooth principal modes without sampling
aliases, away from the special Nyquist treatment; it is not a convergence-order proof for noisy adaptive PIC
trajectories. Consistent grid translation cancels in that approximation.
Actual strong phase dependence therefore points to sampling/alias effects
beyond it, rather than proving one centering convention incorrect.

Independent evidence extends beyond the endpoint:
the [reconstructed operator](PROVENANCE.md#supporting-frozen-state-studies)
matches position-matched printed ASTRA probe-field vectors within 0.15%, and
dense-ellipsoid residual momenta within 0.2% at both tested meshes. Few probes
and limited print precision constrain the field check; neither result
establishes continuum accuracy.

## Why transverse and longitudinal emittance respond differently

For each centred laboratory pair q,p, write

$$
a=\operatorname{Cov}(q,p)/\operatorname{Var}(q),\quad r=p-aq,
\qquad\varepsilon_n=\sigma_q\sigma_r/(mc).
$$

The residual includes nonlinear flow, slice mixing, intrinsic spread and
numerical noise; it is not a measured temperature. For changed forces on the
same particles, the exact force-only sensitivity is

$$
\Delta(d\ln\varepsilon_n/dt)
=\frac{\operatorname{Cov}(r,\Delta F)}{\operatorname{Var}(r)}
=\operatorname{Corr}(r,\Delta F_{\rm res})
 \frac{\sigma_{\Delta F_{\rm res}}}{\sigma_r}.
$$

A constant force or a force linear in the same coordinate cannot contribute
to this instantaneous derivative. Force RMS alone consequently does not
predict emittance error.

The [directional replay](PROVENANCE.md#supporting-frozen-state-studies)
quantifies this on identical EuXFEL particles. At 7.973 ns/32³, the complete
operator change has a longitudinal residual-force RMS **1.54 times larger**
than its x counterpart. Yet σpz,res is 5.97 times larger than σpx,res, and its
force-residual correlation is only 0.291 times as large. Together these yield a
relative εz sensitivity **13.3 times smaller**. The rate changes are
+1.355,+1.343,+0.102 %/ns in xyz; at 10.373 ns they reverse to
−0.453,−0.284,−0.0214 %/ns. Good εz agreement is therefore compatible with a
larger transverse discrepancy. These three instantaneous samples cannot be
integrated into a final-trajectory attribution.

Frozen-factor responses change with state and ordering, and their signs can
reverse. They do not provide a unique endpoint allocation; the native controls
below supply actual trajectory contrasts.

## What the beam distributions and other comparisons tell us

The [flow decomposition](tables/flow.csv)
finds nearly equal transverse coherent slopes but approximately 5.5% smaller
residual momentum in original OPALX. The
[tail/slice budget](PROVENANCE.md#supporting-frozen-state-studies)
places about 91% of the residual-variance gap inside the transverse inner 90%,
with about 58% associated with the longitudinal outer 10%. Approximately
63–64% is within-slice residual variance and 36–37% slice-slope mixing;
centroid mixing is negligible. This localizes the discrepancy but does not
identify its force source. Sparse transverse halo removal would discard
physics without explaining the main gap.

The [final flow diagnostic](tables/flow.csv)
shows phase-0 NGP/corner residual-width gaps of +0.335%,+0.232%,−0.514% in xyz.
At phase ½, x/y residuals are about −0.120%, while sizes are +0.1572%: their
compensation yields the approximately +0.038% emittance gaps. Emittance
agreement therefore does not mean identical sizes or residual momenta.

The dense ellipsoid's transverse momentum variance is more than 99.8% coherent
flow. Emittance measures its small remainder, so close RMS expansion does not
validate emittance. EuXFEL's evolving RF/solenoid and slice structure also
prevents treating that control as a universal correction factor.

## Alternatives constrained, and alternatives still open

The [alternative-cause audit](PROVENANCE.md#validation-scopes)
verifies the shared checkpoint, charge, full cohort, emittance convention,
external-map OFF comparison, physical boundaries and active kernel branch.
Reducing the requested timestep fourfold leaves the endpoint transverse gap
near 4.8%; this does not make Boris and collective RK4 identical.

Both meshes adapt. OPALX freshly deposits and solves at the midpoint, with
kernel refresh when pitches change; ASTRA can rescale/reuse between full
rebuilds. The completed
[no-reuse control](tables/adaptation-effects.csv)
changes ASTRA transverse emittance by at most 0.0036%, leaving the original
gap essentially unchanged. Logged Tcheck ratios are prospective rebuild-check
ratios, not direct measurements of applied intermediate corrections.

The [rotation audit](PROVENANCE.md#supporting-frozen-state-studies)
uses the inferred sequential ASTRA xz/yz algorithm with a fixed CM closure and
consistent inverse field-vector transforms. At the final 32³ state it changes
force RMS by roughly 0.03–0.09%, but the effect increases at 64³ and is sensitive
to extremal NGP boundary assignments. This is not native ASTRA stage capture
and supplies no endpoint bound.

The [kernel audit](PROVENANCE.md#supporting-frozen-state-studies)
identifies binary32 rounding after native ASTRA's one-dimensional FFT passes.
Six-pass-order emulations on three warm states change force RMS by at most
8.47×10⁻⁷ relative and instantaneous emittance rate by 2.99×10⁻⁶ percentage
points/ns. Independent high-precision kernels agree at tested self/near/far
points; the SI coupling differs by 3.13 ppm. These sensitivities are neither
exact native ASTRA FFT reproduction nor an evolved-emittance bound.

The completed [same-build baseline](PROVENANCE.md#six-native-run-bindings)
reproduces the earlier CM original-operator ε/rms history within 8×10⁻¹¹% over
all 10,374 steps. The tested static-build baseline therefore introduces no
meaningful spatial-moment change. The completed
[PRESTEP control](PROVENANCE.md#six-native-run-bindings)
changes transverse emittance by at most 0.0199% over the all-step history;
endpoint ASTRA gaps remain approximately −4.85%/−4.90%. Moving the update to
PRESTEP does not reproduce collective RK4, but this tested timing alternative
does not close the gap.

The completed [native factorial](tables/factorial-effects.csv)
holds those conditions fixed and establishes these endpoint gap changes in
**percentage points**:

| Axis | Assignment with spectral | Gradient/gather with CIC | Assignment–gradient interaction | Phase with NGP/corner |
|---|---:|---:|---:|---:|
| x | −3.1267 | +7.5720 | +0.8854 | −0.4495 |
| y | −3.0528 | +7.5378 | +0.8191 | −0.3732 |
| z | +0.9833 | −0.7430 | −0.3012 | +0.5373 |

The four columns add exactly to the combined change from CIC/spectral+CM to
NGP/corner/phase ½. Assignment alone worsens transverse agreement; the corner
gradient and matching gather alone overshoot ASTRA. Their combination at
phase 0 already satisfies the target, with an interaction that cannot be
discarded. Good longitudinal agreement is not universal: CIC/corner alone
has εz gap −1.191%. These are conditional evolved-trajectory effects,
including feedback, not unique order-independent causal percentages.

The independent [v3/v4 comparison](PROVENANCE.md#validation-scopes)
checks all 10,374 native rows and 25 fixed clocks. Maximum spatial changes
are below 4.37×10⁻⁵%, without bitwise identity. This constrains build/cadence
confounds for the combined result; it does not establish physical accuracy.

## Conservation, accuracy and the justified repair

The [self-force audit](PROVENANCE.md#supporting-frozen-state-studies)
shows a useful advantage of original CIC/spectral: with an even Green kernel,
the odd electric impulse response and matched scatter/gather weights cancel
own-source and total rest-electrostatic internal forces. The fixed-grid replay
retains this cancellation to roundoff. This is a mathematical discrete property,
not a claim of exact cancellation in every native floating-point execution;
signed-kernel evaluation can introduce the tiny conditioning differences
measured in the kernel audit. The mixed research operator does not retain exact
pair reciprocity. On actual EuXFEL states, its replay has a normalized net rest
force imbalance of 0.024–0.093%; its own-force RMS is only 0.47–3.07% of the
complete operator-force contrast and has much smaller emittance-rate effects.
Neither net force nor self-force establishes the endpoint cause. Laboratory
matter momentum also includes exchange with electromagnetic field momentum.

The minimal accepted
[CM repair](README.md#cm-consistency-correction)
uses $p=P/(mc)$, $\bar p=\langle p\rangle$ and
$G=\langle\sqrt{1+p^2}\rangle$. It retains $\beta=\bar p/G$ but uses
$\Gamma_{\rm CM}=G/\sqrt{G^2-|\bar p|^2}$, consistently for rest-frame mesh
stretching, density normalization and field reconstruction. It repairs a
Lorentz-pair inconsistency without new reductions or replacing the spatial
operator, for the supported aligned single-bin OPEN/no-image path. It remains
an instantaneous electrostatic closure, not a treatment
of internal rest currents, retardation or warm-particle simultaneity.

Finally, [32³→64³ refinement](PROVENANCE.md#supporting-frozen-state-studies)
reduces the original cross-code x gap to about 2.19%, while absolute εx changes
by **15.65% in ASTRA and 13.31% in OPALX**. With fixed 100,000 particles, mean
occupied-cuboid population falls from about 29.6 to 3.36 particles per cell.
Neither method is demonstrably converged. Joint mesh/population studies,
independent padding and timestep checks, and known-field tests are required;
duplicating identical particles does not supply a better physical sample.
The original method might reach 1% agreement without the research operator.
Existing evidence neither proves that outcome nor justifies declaring the
ASTRA-like discretization a necessary production fix.


## Questions for the issue discussion

1. Should the consistent CM repair be accepted independently of any optional
   ASTRA-compatibility spatial operator?
2. What conservation, pair-reciprocity and self-force properties should an
   alternative assignment/gradient/gather combination preserve before becoming
   a supported option or changing a default?
3. Which joint mesh/population, independent padding and timestep convergence
   protocol should establish absolute accuracy for this warm RF/solenoid beam,
   rather than only finite-grid cross-code compatibility?
4. Can native ASTRA kernel/grid and stage field/particle dumps provide a stronger
   common-state validation while keeping the field model and observation clocks
   fixed?

The tables demonstrate the stated sampled 1% comparison result. They do not
resolve those production-model and physical-convergence questions.
