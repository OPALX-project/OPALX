This archive contains a C++ implementation of the WFG Toolkit for constructing test problems, as reported in:
  [WFG2005a] Simon Huband, Luigi Barone, Lyndon While, and Phil Hingston. A Scalable Multi-objective Test Problem Toolkit. Evolutionary Multi-Criterion Optimization: Third International Conference, EMO 2005. Proceedings, volume 3410 of Lecture Notes in Computer Science, pages 280-294. Springer-Verlag, March 2005.
and later updated in:
  Simon Huband, Phil Hingston, Luigi Barone, and Lyndon While. A Review of Multiobjective Test Problems and a Scalable Test Problem Toolkit. To appear in: IEEE Transactions on Evolutionary Computation.

The implementation includes the suggested benchmark problems, WFG1-WFG9, as well as the I1-I5 problems.

The source code includes a driver, 'main.cpp', that demonstrates how the toolkit can be employed. The driver simply prints out random Pareto optimal fitness vectors for a problem as specified on the command line. It can also be employed to print out values corresponding to the plots shown in Figure 2 of the [WFG2005a] paper.

Note:
  * The I1-I5 functions were not intended to be used as benchmark functions, although they're pretty functional as is, and I5 represents a very difficult kind of problem.
  * The bias in WFG1 is pretty significant, so much so that the limited precision of double's on most systems may cause problems. This is worth bearing in mind when considering the results, and particularly when comparing results between algorithms that were implemented on different systems.

Recommended values for k and l:
  * Whilst it is largely a matter of choice, letting l=20 should be sufficient for the number of distance related parameters.
  * Whilst suggesting a value for k, the number of position related parameters, is more difficult, we suggest the following:
    - If M == 2, then let l == 4. Otherwise, if M >= 3, let l == 2*(M-1). This ensures there is always at least 4 position related parameters, and never the absolute minimum number of position related parameters. For problems exploring a large number of objectives, it may be desirable to just use the minimum number of position related parameters.

This canonical form of this archive is available at http://www.wfg.csse.uwa.edu.au/.

The directory structure of the archive is as follows:

./	          - contains the driver 'main.cpp'
./Toolkit/    - contains the actual toolkit source code

--
For any queries, email s.huband@ecu.edu.au, or for full contact details, visit http://www.wfg.csse.uwa.edu.au/.

Copyright 2006, Walking Fish Group.