;;; opal-mode-el -- Major mode for editing OPAL files

;; Author: Oscar Roberto Blanco Garcia
;; email : <oscar.roberto.blanco.garcia@cern.ch>
;; Version: 1.0
;; Created: 17.05.2012
;; Keywords: OPAL major-mode

;; This program is free software; you can redistribute it and/or
;; modify it under the terms of the GNU General Public License as
;; published by the Free Software Foundation; either version 2 of
;; the License, or (at your option) any later version.

;; This program is distributed in the hope that it will be
;; useful, but WITHOUT ANY WARRANTY; without even the implied
;; warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
;; PURPOSE.  See the GNU General Public License for more details.

;; You should have received a copy of the GNU General Public
;; License along with this program; if not, write to the Free
;; Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
;; MA 02111-1307 USA

;;; Commentary:
;;
;; This mode is modified from an example used in a tutorial about Emacs
;; mode creation. The tutorial can be found here:
;; https://www.emacswiki.org/emacs/ModeTutorial

;; Add this to your .emacs file to load and bind it to files with extension
;; .opal

;;; Code:

(defgroup opal nil
 "Major mode to edit OPAL files scripts in emacs"
 :group 'languages
)

(defvar opal-mode-hook nil)

;(defvar opal-mode-map
;  (let ((opal-mode-map (make-keymap)))
;    (define-key opal-mode-map "\C-j" 'newline-and-indent)
;    opal-mode-map)
;  "Keymap for OPAL major mode")

(add-to-list 'auto-mode-alist '("\\.opal\\'" . opal-mode))

; optimiser keywords
(defconst opal-font-lock-keywords-optimise
  (list
  '("\\<\\(DVAR\\|DVARS\\|OBJECTIVE\\|OBJECTIVES\\|CONSTRAINTS\\|OPTIMIZE\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (matchingmet).")


;(regexp-opt '("CONSTRAINT" "ENDMATCH" "LINE" "MATCH" "RUN" "START" "TWISS" "VARY") t)

(defconst opal-font-lock-keywords-simul
  (list
   ; These define the beginning and end of each OPAL entity definition
  '("\\<\\(CONSTRAINT\\|ENDMATCH\\|LINE\\|MATCH\\|RUN\\|START\\|TWISS\\|VARY\\)\\>"
 . font-lock-builtin-face)
 )
 "Highlighting expressions for OPAL mode (simul).")

(defconst opal-font-lock-keywords-programflow
  (list
  '("\\<\\(ELSE\\(?:IF\\)?\\|IF\\|MACRO\\|WHILE\\)\\>"
  . font-lock-keyword-face)
  )
  "Highlighting expressions for OPAL mode (programflow).")

;(regexp-opt '("CALL" "CONST" "ESAVE" "EXIT" "HELP" "OPTION" "PRINT" "QUIT" "REAL" "SAVE" "SELECT" "SHOW" "STOP" "SYSTEM" "TITLE" "VALUE") t)

(defconst opal-font-lock-keywords-controlstm
  (list
  '("\\<\\(C\\(?:ALL\\|ONST\\)\\|E\\(?:SAVE\\|XIT\\)\\|HELP\\|OPTION\\|PRINT\\|QUIT\\|REAL\\|S\\(?:AVE\\|ELECT\\|HOW\\|TOP\\|YSTEM\\)\\|\\(?:TITL\\|VALU\\)E\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (controlstm).")

;(regexp-opt '("CCOLLIMATOR" "CYCLOTRON" "DEGRADER" "DRIFT" "ECOLLIMATOR" "HKICKER" "KICKER" "MARKER" "MATRIX" "MONITOR" "MULTIPOLE" "OCTUPOLE" "PROBE" "QUADRUPOLE" "RBEND" "RCOLLIMATOR" "RFCAVITY" "RINGDEFINITION" "SBEND" "SBEND3D" "SEPTUM" "SEXTUPOLE" "SOLENOID" "SOURCE" "STRIPPER" "TRAVELINGWAVE" "VARIABLE_RF_CAVITY" "VKICKER") t)

(defconst opal-font-lock-keywords-elements
  (list
  '("\\<\\(C\\(?:COLLIMATOR\\|YCLOTRON\\)\\|D\\(?:EGRADER\\|RIFT\\)\\|ECOLLIMATOR\\|HKICKER\\|KICKER\\|M\\(?:A\\(?:RKER\\|TRIX\\)\\|ONITOR\\|ULTIPOLE\\)\\|OCTUPOLE\\|PROBE\\|QUADRUPOLE\\|R\\(?:BEND\\|COLLIMATOR\\|FCAVITY\\|INGDEFINITION\\)\\|S\\(?:BEND\\(?:3D\\)?\\|E\\(?:PTUM\\|XTUPOLE\\)\\|O\\(?:LENOID\\|URCE\\)\\|TRIPPER\\)\\|TRAVELINGWAVE\\|V\\(?:ARIABLE_RF_CAVITY\\|KICKER\\)\\)\\>"
  . font-lock-type-face)
  )
  "Highlighting expressions for OPAL mode (elements).")

;(regexp-opt '("BEAM" "DISTRIBUTION" "FIELDSOLVER" "POLYNOMIAL_TIME_DEPENDENCE") t)

(defconst opal-font-lock-keywords-beamspec
  (list
  '("\\<\\(BEAM\\|DISTRIBUTION\\|FIELDSOLVER\\|POLYNOMIAL_TIME_DEPENDENCE\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (beamspec).")

;(regexp-opt '("LMDIF" "MIGRAD" "SIMPLEX") t)

(defconst opal-font-lock-keywords-matchingmet
  (list
  '("\\<\\(LMDIF\\|MIGRAD\\|SIMPLEX\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (matchingmet).")

;(regexp-opt '("ENDTRACK" "SURVEY" "TRACK") t)

(defconst opal-font-lock-keywords-orbit_corr
  (list
  '("\\<\\(ENDTRACK\\|SURVEY\\|TRACK\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (orbit_corr).")

;(regexp-opt '("EPRINT") t)

(defconst opal-font-lock-keywords-plot
  (list
  '("\\<\\(EPRINT\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (plot).")

;(regexp-opt '("CYCLE" "ENDEDIT" "FLATTEN" "INSTALL" "MOVE" "REFLECT" "REMOVE" "SEQEDIT") t)

(defconst opal-font-lock-keywords-seqediting
  (list
  '("\\<\\(CYCLE\\|ENDEDIT\\|FLATTEN\\|INSTALL\\|MOVE\\|RE\\(?:FLECT\\|MOVE\\)\\|SEQEDIT\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (seqediting).")

;(regexp-opt '("A" "ADD" "AMPLITUDE_MODEL" "ANGLE" "APERTURE" "APVETO" "AT" "AUTOPHASE" "B" "BBOXINCR" "BCFFTT" "BCFFTX" "BCFFTY" "BCURRENT" "BEAM_PHIINIT" "BEAM_PRINIT" "BEAM_RINIT" "BFREQ" "BOUNDPDESTROYFQ" "BY" "CALLS" "CATHTEMP" "CENTRE" "CHARGE" "CLASS" "CMD" "CLEAR" "COLUMN" "CONDUCT" "CONST_LENGTH" "CORRX" "CORRY" "CORRZ" "CSRDUMP" "CUTOFFLONG" "CUTOFFPX" "CUTOFFPY" "CUTOFFPZ" "CUTOFFR" "CUTOFFX" "CUTOFFY" "CYHARMON" "CZERO" "DESIGNENERGY" "DK1" "DK1S" "DK2" "DK2S" "DKN" "DKNR" "DKS" "DKSR" "DLAG" "DPHI" "DPSI" "DS" "DT" "DTHETA" "DUMP" "DVOLT" "DX" "DY" "DZ" "E1" "E2" "EBDUMP" "ECHO" "EKIN" "ELASER" "ELEMEDGE" "ELEMENT" "EMISSIONMODEL" "EMISSIONSTEPS" "EMITTED" "ENABLEHDF5" "END_NORMAL_X" "END_NORMAL_Y" "END_POSITION_X" "END_POSITION_Y" "ENDSEQUENCE" "ESCALE" "ET" "EVERYSTEP" "EX" "EY" "FE" "FGEOM" "FILE" "FINT" "FMAPFN" "FNAME" "FORM" "FREQ" "FREQUENCY_MODEL" "FROM" "FSTYPE" "FTOSCAMPLITUDE" "FTOSCPERIODS" "FULL" "GAMMA" "GAP" "GAPWIDTH" "GEOMETRY" "GREENSF" "H1" "H2" "HAPERT" "HARMON" "HARMONIC_NUMBER" "HGAP" "HKICK" "IDEALIZED" "IMAGENAME" "INFO" "INPUTMOUNITS" "INTENSITYCUT" "INTERPL" "IS_CLOSED" "ITSOLVER" "K0" "K0S" "K1" "K1S" "K2" "K2S" "K3" "K3S" "KEYWORD" "KICK" "KN" "KS" "L" "LAG" "LASERPROFFN" "LATTICE_PHIINIT" "LATTICE_RINIT" "LATTICE_THETAINIT" "LENGTH" "LEVEL" "LOGBENDTRAJECTORY" "LOWER" "MASS" "MAXITERS" "MAXR" "MAXSTEPS" "MAXZ" "MBTC" "METHOD" "MINR" "MINZ" "MODE" "MREX" "MREY" "MSCALX" "MSCALY" "MT" "MX" "MY" "NBIN" "NFREQ" "NLEFT" "NO" "NPART" "NPOINTS" "NRIGHT" "NUMCELLS" "OFFSETPX" "OFFSETPY" "OFFSETPZ" "OFFSETT" "OFFSETX" "OFFSETY" "OFFSETZ" "OPCHARGE" "OPMASS" "OPYIELD" "ORDER" "ORIENTATION" "ORIGIN" "OUTFN" "P0" "P1" "P2" "P3" "P4" "PARFFTT" "PARFFTX" "PARFFTY" "PARTICLE" "PATTERN" "PC" "PHI" "PHI0" "PHIINIT" "PLANE" "POLYORDER" "PRECMODE" "PRINIT" "PSDUMPFRAME" "PSDUMPFREQ" "PSI" "PTC" "PYMULT" "PZINIT" "PZMULT" "R51" "R52" "R61" "R62" "RADIUS" "RANDOM" "RANGE" "REFER" "REFPOS" "REPARTFREQ" "RESET" "RFMAPFN" "RFPHI" "RINIT" "RMAX" "RMIN" "ROTATION" "ROW" "S" "SEED" "SELECTED" "SEQUENCE" "SIGMA" "SIGMAPX" "SIGMAPY" "SIGMAPZ" "SIGMAR" "SIGMAT" "SIGMAX" "SIGMAY" "SIGMAZ" "SIGX" "SIGY" "SLPTC" "SPLIT" "STATDUMPFREQ" "STEP" "STEPSPERTURN" "STOP" "SUPERPOSE" "PARTICLEMATTERINTERACTION" "SCALABLE" "SYMMETRY" "T0" "TABLE" "TAU" "TCR1" "TCR2" "TELL" "TFALL" "THETA" "THIN" "THRESHOLD" "TIME" "TIMEINTEGRATOR" "TMULT" "TO" "TOL" "TOLERANCE" "TPULSEFWHM" "TRACE" "TRISE" "TURNS" "TYPE" "UPPER" "VERIFY" "VERSION" "VKICK" "VMAX" "VMIN" "VOLT" "W" "WAKEF" "WARP" "WARN" "WEIGHT" "WIDTH" "WRITETOFILE" "X" "XEND" "XMA" "XMULT" "XSIZE" "XSTART" "Y" "YEND" "YMA" "YMULT" "YSIZE" "YSTART" "Z" "Z0" "ZEND" "ZINIT" "ZSTART" "ZSTOP") t)


(defconst opal-font-lock-keywords-parameters
  (list
  '("\\<\\(A\\(?:DD\\|MPLITUDE_MODEL\\|NGLE\\|P\\(?:ERTURE\\|VETO\\)\\|T\\|UTOPHASE\\)\\|B\\(?:BOXINCR\\|C\\(?:FFT[TXY]\\|URRENT\\)\\|EAM_\\(?:\\(?:P\\(?:HI\\|R\\)\\|R\\)INIT\\)\\|FREQ\\|OUNDPDESTROYFQ\\|Y\\)\\|C\\(?:A\\(?:LLS\\|THTEMP\\)\\|ENTRE\\|HARGE\\|L\\(?:ASS\\|EAR\\)\\|MD\\|O\\(?:LUMN\\|N\\(?:DUCT\\|ST_LENGTH\\)\\|RR[XYZ]\\)\\|SRDUMP\\|UTOFF\\(?:LONG\\|P[XYZ]\\|[RXY]\\)\\|YHARMON\\|ZERO\\)\\|D\\(?:ESIGNENERGY\\|K\\(?:1S\\|2S\\|[NS]R\\|[12NS]\\)\\|LAG\\|P\\(?:[HS]I\\)\\|THETA\\|UMP\\|VOLT\\|[STXYZ]\\)\\|E\\(?:BDUMP\\|CHO\\|KIN\\|L\\(?:ASER\\|EME\\(?:DGE\\|NT\\)\\)\\|MI\\(?:SSION\\(?:MODEL\\|STEPS\\)\\|TTED\\)\\|N\\(?:ABLEHDF5\\|D\\(?:SEQUENCE\\|_\\(?:NORMAL_[XY]\\|POSITION_[XY]\\)\\)\\)\\|SCALE\\|VERYSTEP\\|[12TXY]\\)\\|F\\(?:E\\|GEOM\\|I\\(?:LE\\|NT\\)\\|MAPFN\\|NAME\\|ORM\\|R\\(?:EQ\\(?:UENCY_MODEL\\)?\\|OM\\)\\|STYPE\\|TOSC\\(?:AMPLITUDE\\|PERIODS\\)\\|ULL\\)\\|G\\(?:A\\(?:MMA\\|P\\(?:WIDTH\\)?\\)\\|EOMETRY\\|REENSF\\)\\|H\\(?:A\\(?:PERT\\|RMON\\(?:IC_NUMBER\\)?\\)\\|GAP\\|KICK\\|[12]\\)\\|I\\(?:DEALIZED\\|MAGENAME\\|N\\(?:FO\\|PUTMOUNITS\\|TE\\(?:NSITYCUT\\|RPL\\)\\)\\|S_CLOSED\\|TSOLVER\\)\\|K\\(?:0S\\|1S\\|2S\\|3S\\|EYWORD\\|ICK\\|[0-3NS]\\)\\|L\\(?:A\\(?:G\\|SERPROFFN\\|TTICE_\\(?:\\(?:PHI\\|R\\|THETA\\)INIT\\)\\)\\|E\\(?:NGTH\\|VEL\\)\\|O\\(?:GBENDTRAJECTORY\\|WER\\)\\)\\|M\\(?:A\\(?:SS\\|X\\(?:ITERS\\|STEPS\\|[RZ]\\)\\)\\|BTC\\|ETHOD\\|IN[RZ]\\|ODE\\|RE[XY]\\|SCAL[XY]\\|[TXY]\\)\\|N\\(?:BIN\\|FREQ\\|LEFT\\|O\\|P\\(?:ART\\|OINTS\\)\\|RIGHT\\|UMCELLS\\)\\|O\\(?:FFSET\\(?:P[XYZ]\\|[TXYZ]\\)\\|P\\(?:CHARGE\\|MASS\\|YIELD\\)\\|R\\(?:DER\\|I\\(?:\\(?:ENTATIO\\|GI\\)N\\)\\)\\|UTFN\\)\\|P\\(?:A\\(?:R\\(?:FFT[TXY]\\|TICLE\\(?:MATTERINTERACTION\\)?\\)\\|TTERN\\)\\|HI\\(?:0\\|INIT\\)?\\|LANE\\|OLYORDER\\|R\\(?:ECMODE\\|INIT\\)\\|S\\(?:DUMPFR\\(?:AME\\|EQ\\)\\|I\\)\\|TC\\|\\(?:YMUL\\|Z\\(?:INI\\|MUL\\)\\)T\\|[0-4C]\\)\\|R\\(?:5[12]\\|6[12]\\|A\\(?:DIUS\\|N\\(?:DOM\\|GE\\)\\)\\|E\\(?:F\\(?:ER\\|POS\\)\\|PARTFREQ\\|SET\\)\\|F\\(?:MAPFN\\|PHI\\)\\|INIT\\|M\\(?:AX\\|IN\\)\\|O\\(?:TATION\\|W\\)\\)\\|S\\(?:CALABLE\\|E\\(?:ED\\|LECTED\\|QUENCE\\)\\|IG\\(?:MA\\(?:P[XYZ]\\|[RTXYZ]\\)?\\|[XY]\\)\\|LPTC\\|PLIT\\|T\\(?:ATDUMPFREQ\\|EP\\(?:SPERTURN\\)?\\|OP\\)\\|UPERPOSE\\|YMMETRY\\)\\|T\\(?:A\\(?:BLE\\|U\\)\\|CR[12]\\|ELL\\|FALL\\|H\\(?:ETA\\|IN\\|RESHOLD\\)\\|IME\\(?:INTEGRATOR\\)?\\|MULT\\|OL\\(?:ERANCE\\)?\\|PULSEFWHM\\|R\\(?:\\(?:AC\\|IS\\)E\\)\\|URNS\\|YPE\\|[0O]\\)\\|UPPER\\|V\\(?:ER\\(?:IFY\\|SION\\)\\|KICK\\|M\\(?:AX\\|IN\\)\\|OLT\\)\\|W\\(?:A\\(?:KEF\\|R[NP]\\)\\|EIGHT\\|IDTH\\|RITETOFILE\\)\\|X\\(?:END\\|M\\(?:A\\|ULT\\)\\|S\\(?:IZE\\|TART\\)\\)\\|Y\\(?:END\\|M\\(?:A\\|ULT\\)\\|S\\(?:IZE\\|TART\\)\\)\\|Z\\(?:0\\|END\\|INIT\\|ST\\(?:ART\\|OP\\)\\)\\|[ABLSW-Z]\\)\\>"
  . font-lock-doc-face)
  )
  "Highlighting expressions for OPAL mode (parameters).")

;(regexp-opt '("EALIGN" "EFCOMP" "ERROR") t)

(defconst opal-font-lock-keywords-errordef
  (list
  '("\\<\\(E\\(?:ALIGN\\|FCOMP\\|RROR\\)\\)\\>"
  . font-lock-warning-face)
  )
  "Highlighting expressions for OPAL mode (errordef).")

;(regexp-opt '("AMR" "ANTIPROTON" "ASTRA" "CARBON" "CENTRE" "CLIGHT" "CMASS" "COLLIM" "DEGRAD" "DEUTERON" "DMASS" "E" "ELECTRON" "EMASS" "ENTRY" "EXIT" "FALSE" "FFT" "FFTPERIODIC" "HMINUS" "HMMASS" "INTEGRATED" "LONG-SHORT-RANGE" "MMASS" "MUON" "NONE" "NONEQUIL" "OPEN" "P3M" "PERIODIC" "PI" "PMASS" "POSITRON" "PROTON" "RADDEG" "SAAMG" "STANDARD" "TRANSV-SHORT-RANGE" "TRUE" "TWOPI" "UMASS" "URANIUM" "XEMASS" "XENON" "1D-CSR") t)

(defconst opal-font-lock-keywords-constants
  (list
  '("\\<\\(1D-CSR\\|A\\(?:MR\\|NTIPROTON\\|STRA\\)\\|C\\(?:ARBON\\|ENTRE\\|LIGHT\\|MASS\\|OLLIM\\)\\|D\\(?:E\\(?:GRAD\\|UTERON\\)\\|MASS\\)\\|E\\(?:LECTRON\\|MASS\\|NTRY\\|XIT\\)?\\|F\\(?:ALSE\\|FT\\(?:PERIODIC\\)?\\)\\|HM\\(?:\\(?:INU\\|MAS\\)S\\)\\|INTEGRATED\\|LONG-SHORT-RANGE\\|M\\(?:MASS\\|UON\\)\\|NONE\\(?:QUIL\\)?\\|OPEN\\|P\\(?:3M\\|ERIODIC\\|I\\|MASS\\|\\(?:OSITR\\|ROT\\)ON\\)\\|RADDEG\\|S\\(?:AAMG\\|TANDARD\\)\\|T\\(?:R\\(?:\\(?:ANSV-SHORT-RANG\\|U\\)E\\)\\|WOPI\\)\\|U\\(?:MASS\\|RANIUM\\)\\|XE\\(?:MASS\\|NON\\)\\)\\>"
  . font-lock-constant-face)
  )
  "Highlighting expressions for OPAL mode (constants).")

;(regexp-opt '("TITLE") t)

(defconst opal-font-lock-keywords-stringatt
  (list
  '("\\<\\(TITLE\\)\\>"
  . font-lock-builtin-face)
  )
  "Highlighting expressions for OPAL mode (stringatt).")

;(regexp-opt '("ABS" "ACOS" "ASIN" "ATAN" "COS" "COSH" "EXP" "GAUSS" "LOG" "LOG10" "RANF" "SIN" "SINH" "SQRT" "TAN" "TANH" "TGAUSS") t)

(defconst opal-font-lock-keywords-functions
  (list
  '("\\<\\(A\\(?:BS\\|COS\\|\\(?:SI\\|TA\\)N\\)\\|COSH?\\|EXP\\|GAUSS\\|LOG\\(?:10\\)?\\|RANF\\|S\\(?:INH?\\|QRT\\)\\|T\\(?:ANH?\\|GAUSS\\)\\)\\>"
  . font-lock-function-name-face)
  )
  "Highlighting expressions for OPAL mode (functions).")

;(regexp-opt '("ALFA[1-3][1-3]" "ALFA[1-3][1-3]P" "ALFX" "ALFY" "APER_[1-4]" "BETA[1-3][1-3]" "BETA[1-3][1-3]P" "BETX" "BETY" "DDPX" "DDPY" "DDQ[1-2]" "DDX" "DDY" "DELTAP" "DISP[1-3]" "DISP1P[1-3]" "DISP2P[1-3]" "DISP3P[1-3]" "DMUX" "DMUY" "DPX" "DPY" "DQ1" "DQ2" "DX" "DY" "EIGN[1-6][1-6]" "ENERGY" "GAMMA[1-3][1-3]" "GAMMA[1-3][1-3]P" "KICK" "K[1-6]" "MU[1-3]" "MX" "MY" "N1" "N1X_M" "N1Y_M" "NAME" "ON_AP" "ON_ELEM" "PHIT" "PHIX" "PHIY" "PTN" "PXN" "PYN" "PT" "PX" "PY" "Q1" "Q2" "RE" "RE[1-6][1-6]" "RM[1-6][1-6]" "RTOL" "R[1-6][1-6]" "THETA" "TM[1-6][1-6][1-6]" "TN" "T[1-6][1-6][1-6]" "WT" "WX" "WY" "XN" "XTOL" "YN" "YTOL" "R" "S" "T" "X" "Y") t)

(defconst opal-font-lock-keywords-variables_opal
  (list
  '("\\<\\(A\\(?:LF\\(?:A[1-3][1-3]P?\\|[XY]\\)\\|PER_[1-4]\\)\\|BET\\(?:A[1-3][1-3]P?\\|[XY]\\)\\|D\\(?:D\\(?:P[XY]\\|Q[12]\\|[XY]\\)\\|ELTAP\\|ISP[1-4]\\(?:P[1-3]\\)?\\|MU[XY]\\|P[XY]\\|Q[12]\\|[XY]\\)\\|E\\(?:IGN[1-6][1-6]\\|NERGY\\)\\|GAMA[1-3][1-3]P?\\|K\\(?:\\(?:ICK\\)?[1-6]\\)\\|MU\\(?:[1-3]\\|[XY]\\)\\|N\\(?:1\\(?:[XY]_M\\)?\\|AME\\)\\|ON_\\(?:AP\\|ELEM\\)\\|P\\(?:HI[TXY]?\\|[TXY]N\\|[TXY]\\)\\|Q[12]\\|R\\(?:E\\(?:[1-6][1-6]\\)?\\|M[1-6][1-6]\\|TOL\\|[1-6][1-6]\\)\\|T\\(?:HETA\\|M[1-6][1-6][1-6]\\|N\\|[1-6][1-6][1-6]\\)\\|W[TXY]\\|X\\(?:N\\|TOL\\)\\|Y\\(?:N\\|TOL\\)\\|[RSTXY]\\)\\>"
;  . font-lock-negation-char-face)
  . font-lock-variable-name-face)
  )
  "Highlighting expressions for OPAL mode (variables_opal).")

(defconst opal-font-lock-special_operators
  (list
   '("\\(->\\|:=\\)"
  . font-lock-warning-face)
  )
  "Highlighting expressions for OPAL mode (variables_opal).")

(defconst opal-font-lock-special_constants
  (list
   '("\\(#[es]\\)"
  . font-lock-constant-face)
  )
  "Highlighting expressions for OPAL mode (variables_opal).")


(defconst opal-font-lock-keywords-3
  (append
     opal-font-lock-keywords-optimise
     opal-font-lock-special_operators
     opal-font-lock-special_constants
     opal-font-lock-keywords-programflow
     opal-font-lock-keywords-simul
     opal-font-lock-keywords-controlstm
     opal-font-lock-keywords-elements
     opal-font-lock-keywords-beamspec
     opal-font-lock-keywords-matchingmet
     opal-font-lock-keywords-orbit_corr
     opal-font-lock-keywords-plot
     opal-font-lock-keywords-seqediting
     opal-font-lock-keywords-parameters
     opal-font-lock-keywords-errordef
     opal-font-lock-keywords-constants
     opal-font-lock-keywords-stringatt
     opal-font-lock-keywords-functions
     opal-font-lock-keywords-variables_opal
     opal-font-lock-special_operators
     opal-font-lock-special_constants
  )
 "Balls-out highlighting in OPAL mode.")

(defvar opal-font-lock-keywords opal-font-lock-keywords-3
  "Default highlighting expressions for OPAL mode.")

(defvar opal-mode-syntax-table
  (let ((opal-mode-syntax-table (make-syntax-table c-mode-syntax-table)))

    ; This is added so entity names with unde rscores can be more easily parsed
	(modify-syntax-entry ?_ "w" opal-mode-syntax-table)
	(modify-syntax-entry ?. "w" opal-mode-syntax-table)

	;  Comment styles are same as C++
	(modify-syntax-entry ?/ ". 124 b" opal-mode-syntax-table)
	(modify-syntax-entry ?* ". 23" opal-mode-syntax-table)
	(modify-syntax-entry ?\n "> b" opal-mode-syntax-table)
	(modify-syntax-entry ?! "< b" opal-mode-syntax-table)
	(modify-syntax-entry ?' "|" opal-mode-syntax-table)
	opal-mode-syntax-table)
  "Syntax table for opal-mode")

;;; ### autoload
(defun opal-mode ()
  "Major mode for editing OPAL script files"
  (interactive)
  (kill-all-local-variables)
  (setq mode-name "OPAL")
  (setq major-mode 'opal-mode)
  (setq comment-start "//")
;  (use-local-map opal-mode-map)
  (set-syntax-table opal-mode-syntax-table)
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '(opal-font-lock-keywords nil t))
;; Set up search
  (add-hook 'opal-mode-hook
     (lambda ()  (setq case-fold-search t)))
  (run-hooks 'opal-mode-hook)
)
(provide 'opal-mode)

;;; opal-mode.el ends here