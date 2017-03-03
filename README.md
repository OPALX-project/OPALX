This is OPAL >= 1.9.9. Be ware of the fact that the documentation is not
anymore a seperate project. The documentation is in this repository in the 
directory 'doc' (https://gitlab.psi.ch/OPAL/src/tree/master/doc). On machines 
with LaTex the User's Reference Manual can be built in the same directory as the 
code by configuring cmake using -DENABLE_DOCUMENTATION=1 to enable the compilation
and -DEXCLUDE_DOC_FROM_ALL=1 to exclude it from compiling when make is called 
without target. The User's Reference Manual can then be built by calling 'make doc'.