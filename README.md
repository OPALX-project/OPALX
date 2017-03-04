This is OPAL >= 1.9.9. Be aware of the fact that the documentation is not
anymore a seperate project. The documentation is in this repository in the 
directory 'doc' (https://gitlab.psi.ch/OPAL/src/tree/master/doc). On machines 
with LaTeX the User's Reference Manual can be built in the same directory as the 
code by configuring cmake using -DENABLE_DOCUMENTATION=1 to enable the compilation.
The manual can then be built by calling 'make doc'. If you want 
to compile the manual every time you compile the code then set 
-DEXCLUDE_DOC_FROM_ALL=0 .