# Frozen-study provenance and validation scope

This portable archive accompanies [OPALX issue #636](https://github.com/OPALX-project/OPALX/issues/636).
It copies numerical tables and figures from the completed EuXFEL investigation.
Logical source names below identify retained study records without publishing
machine-specific paths. Hashes identify evidence; they do not make unbundled
raw data or executables available from this directory.

## Executable and source identities

| Identity | Revision or SHA-256 |
|---|---|
| Original executed OPALX source revision | `7c8a8ba2e6e61fd36c07657c6dabd6fdd3b5d770` |
| Frozen research source base, before its retained research patch | `6b819850e9574f34921b95638b0360a775c4d76a` |
| Frozen research OPALX executable | `ea3cd05f30a2358c47995a3992c27478d475bef12c2305efba0781fd1089807a` |
| ASTRA executable | `d72a57c12bfc46af253c5bf3cd57751cc23b0324d1bb899f602bf08eea6e69d7` |
| Shared gun-exit checkpoint | `d18ecd0f1405af0aa53137f47351212286e5f3943db9a6a5e02f92e2f36dc71a` |
| Inspected pinned IPPL revision | `704932fdaa390f883524cf53588c92386294add4` |
| Frozen research manifest | `2f49f78952868abce0cb1b93db9a3a16f51d9cc2337c4d00988ffdbd39071411` |
| Frozen research patch | `9c4f9868fb88243a23cd7f0ba991b69e3593e95a912c0a5c159e2f6fe0256d18` |

The six controlled warm runs use this one frozen executable, one MPI rank and
four OpenMP threads. The build links OPALX statically; external runtime-library
identities are retained in the original manifest. Its research mode is selected
explicitly by fixed process configuration. The archived source patch included
research numerical operators; those research sources and build products are
not copied into this evidence bundle or asserted to be production defaults.

**The tables describe frozen study executables before the current issue-branch
integration and stationary-limit robustness amendment.** The branch integrates
an algebraically equivalent robust `hypot` evaluation and new ±15 and
±15-with-10⁻¹²-perturbation stationary/counterstreaming tests. The archive makes
no claim that its trajectories were regenerated with that later source.
No new tracking occurred to prepare this commit. Fresh current-branch tests
passed all eight selected suites, with commands and the before-fix stationary
regression control recorded in [VALIDATION.md](VALIDATION.md). That test result
is separate from historical frozen-study validation. Associated manual pages
are in the separate `opalx-manual` repository and are not committed here.

## Validation scopes

- All six warm controls completed with exit code zero, no recorded failure,
  unchanged executed inputs and the full 100000-particle, −250 pC cohort.
  Complete/prepared/input, binary, analysis, scalar/comparison CSV, verifier,
  dump and manifest identities were checked by the original independent audit.
- Each new native control has 104 saved H5 states, so its six-control diagnostic
  audit covers 624 states. The separate final-flow audit covers another stated
  selection: Original 52 + historical CM 52 + five spatial cases ×104 =624;
  it does not include PRESTEP flow. ASTRA's 25 full particle states are checked.
- Primary OPALX emittance/RMS values are all-step 1 ps native scalars with the
  same mechanical definitions verified against complete saved particles.
  Established STAT/H5 relative gates are 1e−8; native energy gates are 1e−10.
  Reference-axis projection is bounded rather than silently equated to lab axes.
- Actual warm deposited RHS and rest pitches independently reproduce potential
  and interior E below the unchanged 1e−8 relative L2 gates. This verifies the
  warm Poisson interface and derivative, not deposition from unavailable native
  midpoint particles, gather on those particles or the complete RF/solenoid
  Boris kick. Separate five-mode cold controls have exact-density/field/Boris
  checks. No exact common-state **warm** native Boris replay is claimed.
- The same-build original-operator baseline reproduces historical CM spatial
  scalars within 7.424e−11% over 10374 native rows. PRESTEP timing effects on
  transverse emittance are at most 0.019865% over those rows. The ASTRA no-reuse
  control changes transverse emittance at the 25 clocks by at most 0.003604%.
- The v3/v4 combined-method comparison checks all 10374 native rows and 25 fixed
  clocks; maximum spatial changes are below 4.37e−5%, without bitwise identity.
  Its status verifies provenance/diagnostics, not a numerical accuracy PASS.
- Historical pure validation includes seven focused research helper cases,
  27 analyzer/collector/grid cases, six repeatability cases, seven exact-energy
  cases and eight final-flow cases. The accepted CM study also passed the
  selected eight CTest suites, including the dedicated two-rank/empty-rank CM
  case. None of those counts claims a fresh test of the current branch.

All clocks, particles, species masses, physical maps and numerical acceptance
gates remain fixed. Percentage-point trajectory contrasts retain interactions
and feedback; they are not order-independent causal fractions. Independent
warm electrostatic closure, spatial/population convergence and conservation
questions remain distinct from meeting the sampled 1% comparison criterion.

## Six native run bindings

Each fingerprint below was read from the completed current frozen study.
`conversion report` gives the physical clock/origin metadata supplied with the
unchanged conversion; `analysis` is the checked nine-observable report.
The raw H5/STAT and native-grid files are fingerprinted but are not distributed
here. Portable tables/figures have their own complete `SHA256SUMS` manifest.

### euxfel-baseline

Operator `baseline`, update `MIDPOINT`; 25 comparison clocks, 104 H5 states, full cohort.

| Logical record | SHA-256 |
|---|---|
| completion | `d9c853388b5a2bf845850891233719602e2f9f85f644675a1e6434f1e5839be5` |
| preparation | `57d60dcace68a902efbfe4d44e6379507364c9f969911aa7b911dd004b7ee05c` |
| conversion report | `5848d77c44a9edd8e87e51fad56ae03317fd71ee62d0722919d8e1dd6fee93b8` |
| analysis | `92a543b9d05931e345c8764cc75112fe6aa779a2648e1e1a924685bd5a92b611` |
| comparison CSV | `61445f61f4f5448a3dbd40f24607d8298e338866d451a3a54271802990e080bb` |
| native scalar CSV | `a9f2569670c057fd5ebae784580d6f76154469ca0bcc907a4ad9801376b07d40` |
| warm verification | `78bdc1f760b560ccae6402e3d072151bcedf5690aaad4b15de7dbbb8d4432da5` |
| warm grid dump | `f7a1f838af22829f5597dc9cb2eb35966e71286e5200e12e8e75b8a82a5b4163` |
| native STAT output | `c1c0e6bd7cfa1979f8014e515a8317c57f9aa95bcd77591410a1a7dc7e4cecbb` |
| native H5 output | `72b9591797e6b8b1d44a1b32e8032735e0b3b77ef3e025654a8b971a4af95a94` |

### euxfel-ngp-spectral

Operator `ngp-spectral`, update `MIDPOINT`; 25 comparison clocks, 104 H5 states, full cohort.

| Logical record | SHA-256 |
|---|---|
| completion | `5768b0dddb094ee88019af3998aec6839c32ef778221bf3f08358a8f49ceb37b` |
| preparation | `7462bbefdcb057556ca199b46da297775cfd9499ff633551c1fe76935194fa73` |
| conversion report | `5848d77c44a9edd8e87e51fad56ae03317fd71ee62d0722919d8e1dd6fee93b8` |
| analysis | `acefbc092bd1e16ab3e2ecaeb68b388d4a918fceb25b0cb1be27daec580a99f0` |
| comparison CSV | `2b8bd287119b6eb18e88e55b481cbab624d65fce814528d6ca6cdc01d1bca785` |
| native scalar CSV | `1152cef3eca07d6ba4e5f13032fa171ffba7246a1c1eb540e4dc95398851357d` |
| warm verification | `c8ffa10a8e5cee33cad7f29b5b5b991790ed30257b76428944c79604747749db` |
| warm grid dump | `04e4c1dba46a9250fb151ccd159dfb07c86e91417954c96b513da46650cf0440` |
| native STAT output | `37a94df1224e347f6dd204cf09351efa2f7e2d3197a708c4925c080958d80377` |
| native H5 output | `b38dc284f78df244ed84f06d95a1f2dc4abffc28aa126e7f27a93f0e8f7172db` |

### euxfel-cic-corner

Operator `cic-corner`, update `MIDPOINT`; 25 comparison clocks, 104 H5 states, full cohort.

| Logical record | SHA-256 |
|---|---|
| completion | `525e3668178903e543876081b795f65fe542bfeb7c8f8ad657d8bc04acbbc5b5` |
| preparation | `229cbe2e94fd5839125e10d464c358e2ec1b2bf843a500fa1ed461aada585ff0` |
| conversion report | `5848d77c44a9edd8e87e51fad56ae03317fd71ee62d0722919d8e1dd6fee93b8` |
| analysis | `956cffa9c3d51308bbb3443118d8d1fbefb50677486588255a3a1c55bdf9a05f` |
| comparison CSV | `cbb2c67866db4083cf06087b0eb48d5ae2c83b1f7dd12367204a5297e8b730f9` |
| native scalar CSV | `82535fa37486091381edc9a6950e6322ffb921bc1af68de6525f6c9fdf22c206` |
| warm verification | `2cffa5186648b516959e39d6f3efa06574644f0cc5e54ba3b0a0091778c146ab` |
| warm grid dump | `0d86bb8001735093bf2526d161c05c982311b8a0012fe0da98a25226612063aa` |
| native STAT output | `85942bf4465de75485c02b215fb650591342f2a22ef32524cb5ca4e517032ac2` |
| native H5 output | `9de7bed0a9ed575d7c554240287a8291bb5befa59725ba0359d413a7100de3e2` |

### euxfel-ngp-corner

Operator `ngp-corner`, update `MIDPOINT`; 25 comparison clocks, 104 H5 states, full cohort.

| Logical record | SHA-256 |
|---|---|
| completion | `df52557870c045941cad860f4ed8aa52a0a31f25e4ed779262062ec23f0180c2` |
| preparation | `f618b1f9aea43e68ebe7f4a67d85dcfdde5dcfa9617d5dcb92122732703793d8` |
| conversion report | `5848d77c44a9edd8e87e51fad56ae03317fd71ee62d0722919d8e1dd6fee93b8` |
| analysis | `f31aa5b9a994ffe46c8283773fafe3446d5284532dbf99ca98c0da9528b09a1d` |
| comparison CSV | `51447c6facc2400d6094c4ce9061733df761eb89768d693f604d4d10596f95eb` |
| native scalar CSV | `ed347bde783b90d1355b458b811f1e5db89e9bc014cf3f874b37c3555a041fc8` |
| warm verification | `27288fad92c4678e4704fb04ac209bb68022b70d40d94bfc1ca9e5fbbad61107` |
| warm grid dump | `d8d4c5c6138e70cde292eabfe0b927c4d6a9cbb6a7b9c669adb90cac0dde7827` |
| native STAT output | `bfa1b506cc328655c72a76bde2c3aebcd0627b566bc7e3ac4ac98c5a3fc8f794` |
| native H5 output | `f2179e793e5e2edb8b7d05fd64e0824eeaa1cf792c474317013a2d99ae331d0f` |

### euxfel-astra

Operator `astra`, update `MIDPOINT`; 25 comparison clocks, 104 H5 states, full cohort.

| Logical record | SHA-256 |
|---|---|
| completion | `8a7612928077ab1a7cbc58a062c4bd8c73b59c058774d4565e3f36ad7fd8f8db` |
| preparation | `384e47cf7cdc6a66012c828bc7ceffc5e58748d28674d4648b62a7e6a04f441e` |
| conversion report | `5848d77c44a9edd8e87e51fad56ae03317fd71ee62d0722919d8e1dd6fee93b8` |
| analysis | `feb56e985f590aa85f882607f8a2ebd7c86149ad850ac89cab6c21c63340b450` |
| comparison CSV | `e771761bb7ab23347777f9997721cb3003717fb5f9524178ee9a404f2cb77a7f` |
| native scalar CSV | `c17c3621a4350568e4d69eee62da5bb786b493e16fc8c25ba32e952b29eec352` |
| warm verification | `404058abe3246dc4c097bf81d6196a03a854db48ca9c9d19e3733ff37bce93cc` |
| warm grid dump | `749321108152e0c9707a679c506a25bf048526355d99d82a40e1cc6db4ee4ce0` |
| native STAT output | `6340da46e3fcff3454fb2993b776d186e63a0d8642a75e96c56881cf4346c389` |
| native H5 output | `fdfa3eb6507f7b907f3789c1ce72577da4f7f3aae1dc22c38a82e8de6819956c` |

### euxfel-baseline-prestep

Operator `baseline`, update `PRESTEP`; 25 comparison clocks, 104 H5 states, full cohort.

| Logical record | SHA-256 |
|---|---|
| completion | `b0be94d9cfbdd40f5adb80a17351c4af0ab9c4b1ca2ebf5a785d192350b1783c` |
| preparation | `ff0ec690a45204958a802a08bee11c00dc059af7a80f10e68aaca5149c2dfbfa` |
| conversion report | `5848d77c44a9edd8e87e51fad56ae03317fd71ee62d0722919d8e1dd6fee93b8` |
| analysis | `9df8fd05017c7aff177bdb7c54b9ad6918ca3169175756bdbc69d4fc033960b5` |
| comparison CSV | `2bf2ea35d347933837c0327360263ac19c7ff1202382e03d88bc1534d0459c32` |
| native scalar CSV | `8169352ea2931c74371912c630317cba3b6d79bb28a538cbfbe23de10e782db0` |
| warm verification | `e8a162c6a5c5f4f7c0a5e4c1784551e73a43dcd58ac7a75d0b580bb921df27ce` |
| warm grid dump | `d13484cb14ab1789dfcd0742040ec900eae2f4ff1420e62504993eeff8d422b3` |
| native STAT output | `71b861e188b504f3ca663e20a0c9365a8212ac3fa588b9eb8f01e3bd4ef2dcd6` |
| native H5 output | `0813fdfeb53e1011b3735277fa30643b728be948151b5625a51b5605301a247b` |

## Frozen research source fingerprints

These identify the actual source snapshot retained with the research executable,
not the later current-branch implementation. Filenames below are repository
relative; sources are not duplicated here.

| Frozen source | SHA-256 |
|---|---|
| `src/SpaceCharge/CartesianPIC3D/CartesianDomainUpdater.cpp` | `90d7e5cd636ee73cabdc1a92cd6588ca7aa174f9e12274d8ce73c7e5caa6222a` |
| `src/SpaceCharge/CartesianPIC3D/CartesianPIC3DAlgorithm.cpp` | `3d984961089c44d9abd7024e960cf376a7c303b4c0f29e1a8e6beac776af5d8f` |
| `src/SpaceCharge/CartesianPIC3D/CartesianPIC3DAlgorithm.h` | `51c265deeef39fb0459c8605915c7c64a6cd7f1440a8398c715d55744750770d` |
| `src/SpaceCharge/CartesianPIC3D/ParticleMeshFieldTransfer.cpp` | `fdac60007a94af58ec1797524c5b1f521305e59ed09226a9655902e78cadeb0f` |
| `src/SpaceCharge/CartesianPIC3D/RelativisticFieldComposer.cpp` | `c0d26ea294b779707394dc657fd18fcfdb4e7f078edbc0848f13717da6871749` |
| `src/SpaceCharge/CartesianPIC3D/RelativisticFieldComposer.h` | `1f805aae1d588a87ab716dbc772a42df524168847c69d3897acc62fa69e5ce56` |
| `src/SpaceCharge/CartesianPIC3D/ResearchAstraStencil.h` | `9658515cd933f9f7152f2c564a57cb1c7f0ffbee297ba79ced0e6e3772f1a34d` |
| `src/SpaceCharge/CartesianPIC3D/ResearchOperatorConfig.h` | `77347006b71f1566dd5a2a0ec00014b06df0164dd0d5632a7f30179dcba547e4` |
| `src/SpaceCharge/Poisson/OpenPoissonAdapter.h` | `d9e3957f7840f4b40015f0b2f375a845c729542ba2c51e66c58defbb86dd6a08` |
| `unit_tests/SpaceCharge/CMakeLists.txt` | `fcd9564e6ebec582574ac46468e5476e26cd91a4224718c02e74b4cf675961db` |
| `unit_tests/SpaceCharge/TestCartesianPIC3DAlgorithms.cpp` | `5847db23cc3de02762c0d5e8b47af0e2b67f20d1f8b6fd1197aa80ecd74b0210` |
| `unit_tests/SpaceCharge/TestCartesianPIC3DFields.cpp` | `d1952e004b848f22cc53cc0f5195c4d6fb65cd9be1f8bab8c55bcdca833dff12` |
| `unit_tests/SpaceCharge/TestParticleBinTraversal.cpp` | `e6987e1a9ce9130ea1ec400d4f7393d809cf4fdf3b5ff022060f6872a1dd49b4` |
| `unit_tests/SpaceCharge/TestResearchAstraStencil.cpp` | `6ca2b727f928f27eb98380c7f44c4bcb2be277c3bd4d50039a33772e80d30a2a` |

## Analysis and report identities

These logical fingerprints bind the study's parser/analysis and completed
scientific reports. Some full reports and raw validation records are not
bundled because they contain machine-local links or reference large native
files. The portable README/discussion are rewritten summaries of those sources;
they are hashed as new archive files rather than pretending to be byte-identical
copies of the original reports.

| Logical source | SHA-256 |
|---|---|
| native case analyzer | `f2e1b49a6b962f18c270ec4c60ee440c4ba2052b822832f255fb5345be8b9d61` |
| warm grid verifier | `f570c34df59f37ffb4dc88cca967ca02bfe4b79dac3ff704b30c97e16a6bb087` |
| full collector | `de423e49c2e87afb22e1d759bff9598ae6cb28eae273c41de2db12a8127af6d0` |
| final report writer | `491e1bd68340befb07bb073233e6302825a5de1ca3df6fb0ac3dc01954ef60de` |
| final native scientific report | `289489003ef7f42f2765702038f25aa49c9c160cfe391e081bbdf6e9de962a3e` |
| final scientific explanation | `23b70b5b4c9cd092dae429fddfa66e38f115a56e62a568102371f1c6c1eb2357` |
| alternative-cause audit | `d9ff142f5320edd0b24560a6e399770699925a97fff2ad9b7f0d2d44eca59439` |
| independent completion audit | `1b0e25cc8a9964212c56d3aef2bc2bb4f02c6028e790275faf2af6962264c039` |
| full collector manifest | `0ed5beb30d9349eebdfc70b2b0217848a0065b305006f6c072908c33031134df` |
| machine final summary | `d049320e5239d09d2f56cb7fb73ab3e01cc05f307467e84c151ece133cdb07cb` |
| STAT validation | `a0fc9806e62e8f5612cf6eb29c80af8f824b6034735dc84e0d8eec81014f9f35` |
| energy validation | `97f0dae1540138af50eeb3dbfa2c0d90db58502e0c7715bd1196769185f0d4ff` |
| final-flow script | `09303526c1e7c17f4bcf5116ea4553699aed4ee72f0f55d145526f8cf718b60c` |
| final-flow pure tests | `2e1bda2ce44ea3f7d6f6c9b2147420005fc1cec25086d12892a9840339562e9f` |
| final-flow report | `ab0d432907edd081291d9696f15990f2adc395b47cb0fe195ab950e6433dd0b7` |
| final-flow provenance | `320930f61a10a764e5cd3159e93c5b765039712c2f20e80d22f7360375e0b4ad` |
| kernel and FFT audit | `ded778f764e2b230e90ac2271998b48fe02f5084c572dc7561c4f52523989185` |
| rotation audit | `e181ddd5cb6cf0a91cee60bbf34930437c97e9bee5ab7ae2722c045c062167f2` |
| directional frozen-force audit | `f9394a0d2c435d689bcbf3cde57357123ce382cbbc4d036bbb0f2685325ddda3` |
| self-force audit | `bca69dd5b00045c3c8cd00d3a955f181a7e6764cb772719c273c0bc597ad9d26` |
| tail and slice audit | `4c1fe1590f7460c34bea5b7f8f03c901224453abba883c0fd342942bb2601f81` |
| legacy resolution study | `8296ecc57af8f45b7f6eb8b0a47c15be6fba3f7367343a8ff30554daf61cae2d` |
| adaptation independent audit | `2d42719d9e18df939dcbfb1aee2cd0bb5ae6de8bb5cb8b83e2fb8034be86f1e2` |
| adaptation native scalar report | `8b9ec41dccfe6c6f695cf3c73bf0919c61e72a2402b4d4964b16a8845a039577` |
| v3/v4 repeatability record | `ae890508df2fa3f00562cb321b0e81ef1fb08d77480851c5955cddfece2f21e5` |

## Supporting frozen-state studies

The following supporting results are summarized in the discussion and identified
by the report hashes above. Their complete raw arrays/disassembly are not in
this small portable bundle; they are not substitutes for the copied evolved
trajectory CSVs.

- **Kernel/FFT:** analytic cuboid-integrated `1/r`, finite origin, doubled-grid
  reflection and convolution normalization match mathematically at common
  pitches. Independent 60-digit cuboid/Gaussian-Laplace quadrature agrees after
  binary64 export at 21 tested self/near/far points. Exact native ASTRA kernel
  arrays remain unavailable. The coupling convention differs by about 3.13 ppm.
  Native ASTRA rounds after each one-dimensional FFT pass; six-order controlled
  emulations on three warm states change force RMS by at most 8.47e−7 relative
  and instantaneous relative emittance rate by 2.99e−6 percentage points/ns.
  These are fixed-state sensitivities, not exact native FFT reproduction or an
  evolved-trajectory bound.
- **Directional response:** the same-particle residual-force covariance can
  make the longitudinal relative emittance response far smaller than the
  transverse response even when longitudinal force differences are larger.
  Three sampled instantaneous rates change sign; they are not integrated as
  predicted endpoint differences.
- **Rotations:** inferred sequential ASTRA covariance rotations change final
  32³ force RMS by about 0.03–0.09%; sensitivity grows at 64³ and extremal NGP
  assignments depend on coordinate-path roundoff. No exact native-stage state
  or endpoint rotation bound is claimed.
- **Conservation:** original matched CIC/spectral with an even kernel retains
  own-source and total internal rest-electrostatic cancellation to roundoff in
  fixed-grid replay. The experimental mixed operator does not retain exact
  pair reciprocity. Its sampled normalized net rest-force imbalance is
  0.024–0.093%; own-force RMS is 0.47–3.07% of the full operator contrast. Those
  facts do not alone identify the endpoint cause or laboratory field-momentum
  exchange.
- **Distribution:** approximately 91% of the residual-variance gap is in the
  transverse inner 90%; approximately 58% is associated with the longitudinal
  outer 10%. Within-slice and slice-slope-mixing terms both contribute; removing
  a sparse transverse halo does not explain the discrepancy. These are
  descriptive budgets, not causal isolated force assignments.
- **Legacy refinement:** the earlier, pre-CM 32³→64³ sequence at fixed 100000
  particles changes absolute transverse emittance by about 15.65% in ASTRA and
  13.31% in OPALX. The original cross-code x gap falls to about 2.19%. Those
  values use the retained coarse central-covariance measurement and should not
  be mixed with the new 1 ps STAT primary method. Mean occupied-cuboid-cell
  population falls about 29.6→3.36; padding also changes fractionally, so this
  is not a demonstrated fixed-domain continuum sequence.

## Portable artifact origins

All copied payloads have identical source bytes and exact CSV numeric identity.
The logical source families are `factorial collection`, `final projected flow`
and `adaptation control`; no raw source path is needed to interpret the data.
The filename map is:

| Portable payload | Frozen source family / filename |
|---|---|
| `figures/emittance.png` | factorial collection / `emit.png` |
| `figures/emittance.pdf` | factorial collection / `emit.pdf` |
| `figures/rms.png` | factorial collection / `rms.png` |
| `figures/rms.pdf` | factorial collection / `rms.pdf` |
| `figures/energy.png` | factorial collection / `energy.png` |
| `figures/energy.pdf` | factorial collection / `energy.pdf` |
| `figures/emittance-factor-effects.png` | factorial collection / `emittance-factor-effects.png` |
| `figures/emittance-factor-effects.pdf` | factorial collection / `emittance-factor-effects.pdf` |
| `figures/endpoint-factor-effects.png` | factorial collection / `endpoint-factor-effects.png` |
| `figures/endpoint-factor-effects.pdf` | factorial collection / `endpoint-factor-effects.pdf` |
| `figures/residual-momentum.png` | final projected flow / `residual.png` |
| `figures/residual-momentum.pdf` | final projected flow / `residual.pdf` |
| `figures/coherent-flow.png` | final projected flow / `coherent.png` |
| `figures/coherent-flow.pdf` | final projected flow / `coherent.pdf` |
| `figures/adaptation-grid-geometry.png` | adaptation control / `saved-grid-geometry.png` |
| `figures/adaptation-grid-geometry.pdf` | adaptation control / `saved-grid-geometry.pdf` |
| `figures/adaptation-reuse-effect.png` | adaptation control / `reuse-effect.png` |
| `figures/adaptation-reuse-effect.pdf` | adaptation control / `reuse-effect.pdf` |
| `tables/observables.csv` | factorial collection / `observables.csv` |
| `tables/historical-references.csv` | factorial collection / `historical-reference-observables.csv` |
| `tables/endpoint.csv` | factorial collection / `endpoint.csv` |
| `tables/factorial-effects.csv` | factorial collection / `factorial-effects-percentage-points.csv` |
| `tables/flow.csv` | final projected flow / `flow.csv` |
| `tables/flow-endpoint.csv` | final projected flow / `endpoint.csv` |
| `tables/adaptation-effects.csv` | adaptation native-scalar supplement / `effects.csv` |

## Packaging validation

Copied figure/CSV bytes were matched to their frozen originals; all seven CSVs
were parsed with round-trip floating-point decoding and checked for exact
numeric/data-frame identity. Every local Markdown link is required to resolve
inside this portable bundle, with only the public issue URL external. Text and
binary payloads are scanned for private machine paths. `SHA256SUMS` lists every
bundle file except itself and can be checked with `shasum -a 256 -c SHA256SUMS`.
This packaging checks an archive; it does not rerun simulation or establish
new physical accuracy. Regenerate the manifest after any later validation-note
addition or report edit.
