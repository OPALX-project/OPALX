AC_INIT([OPAL],[1.1.8],[opal@lists.psi.ch])

AC_PROG_CC([mpicc])
AC_PROG_CXX([mpicxx])

AC_CONFIG_AUX_DIR([config])
#disable f77 tests
m4_defun([_LT_AC_LANG_F77_CONFIG], [:])
# library creation
AC_PROG_LIBTOOL
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([./config.h:./config.in])
AM_INIT_AUTOMAKE([opal],[1.1.9])

AC_LANG(C++)
AC_PROG_INSTALL

AC_PATH_PROG(svnversioncommand, svnversion)
if test "X$svnversioncommand" = "X" || test `$svnversioncommand -n '.'` = "exported"; then
AC_DEFINE_UNQUOTED(SVN_VERSION, ["`cat saved_svn_version.txt`"], [repository svn version])
else
AC_DEFINE_UNQUOTED(SVN_VERSION, ["`svnversion -n`"], [repository svn version])
`svnversion -n > saved_svn_version.txt`
fi

CFLAGS="${CFLAGS=}"
CPPFLAGS="${CPPFLAGS=}"

IPPLDEFS="-fpermissive -w -DIPPL_MPI\ 
 -DPARALLEL_IO\
 -DMPICH_SKIP_MPICXX\
 -DIPPL_DEBUG\ 
 -DIPPL_GCC\
 -DIPPL_DONT_POOL\
 -DIPPL_USE_XDIV_RNG\
 -DIPPL_LINUX\ 
 -DIPPL_NO_STRINGSTREAM\
 -DPETE_BITWISE_COPY\
 -DIPPL_HAS_TEMPLATED_COMPLEX\
 -DIPPL_USE_STANDARD_HEADERS\
 -DIPPL_USE_PARTIAL_SPECIALIZATION\
 -DIPPL_STDSTL\
 -DIPPL_LONGLONG\
 -Drestrict=__restrict__ -DNOCTAssert"


# this macro is used to get the arguments supplied
# to the configure script (./configure --enable-debug)
# Check if we have enable debug support.
debug_default="no"
AC_MSG_CHECKING(whether to enable debugging)
AC_ARG_ENABLE(debug, [  --enable-debug=[no/yes] turn on debugging
                     [default=$debug_default]],, enable_debug=$debug_default)

# Yes, shell scripts can be used
if test "x$enable_debug" = "xyes"; then
        CFLAGS="$CFLAGS -g -DDEBUG"
        CPPFLAGS="$CPPFLAGS -g -DDEBUG"
        AC_MSG_RESULT(yes)
else
        CFLAGS="$CFLAGS -O3 "
        CPPFLAGS="$CPPFLAGS -O3 "
        AC_MSG_RESULT(no)
fi


AC_ARG_WITH(ippl-includedir,
        AC_HELP_STRING(
                [--with-ippl-includedir=dir],
                [ippl include files in dir]
        ),
        [
                test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -I$withval -I$HDF5HOME/include ${IPPLDEFS}"
        ]
)
AC_ARG_WITH(ippl-libdir,
        AC_HELP_STRING(
                [--with-ippl-libdir=dir],
                [ippl library libippl.a in dir]
        ),
        [
                test X$withvalue != Xno && LDFLAGS="$LDFLAGS -L/usr/lib64 -L$withval -L$HDF5HOME/lib"
                HAVE_IPPL_=true
        ]

)
AC_ARG_WITH(classic,
	AC_HELP_STRING(
		[--with-classic=dir],
		[classic files in dir]
	),
	[
		if test X$withvalue != Xno; then
			CPPFLAGS="$CPPFLAGS -I$withval"
			LDFLAGS="$LDFLAGS -L$withval"
			HAVE_CLASSIC_=true
		fi
	]
)
AC_ARG_WITH(classic-includedir,
	AC_HELP_STRING(
		[--with-classic-includedir=dir],
		[classic include files in dir]
	),
	[
		test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -I$withval"
	]
)
AC_ARG_WITH(classic-libdir,
	AC_HELP_STRING(
		[--with-classic-libdir=dir],
		[classic library libCLASSIC.a in dir]
	),
	[
		test X$withvalue != Xno && LDFLAGS="$LDFLAGS -L$withval"
        HAVE_CLASSIC_=true
	]
)

AC_ARG_WITH(h5hut,
	AC_HELP_STRING(
		[--with-h5hut=dir],
		[h5hut prefix]
	),
	[
		test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -I$withval/include"
		test X$withvalue != Xno && LDFLAGS="$LDFLAGS -L$withval/lib"
	]
)

AC_ARG_WITH(h5hut-includedir,
	AC_HELP_STRING(
		[--with-h5hut-includedir=dir],
		[h5hut include files in dir]
	),
	[
		test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -I$withval"
	]
)
AC_ARG_WITH(h5part-libdir,
	AC_HELP_STRING(
		[--with-h5hut-libdir=dir],
		[h5hut library in dir]
	),
	[
		test X$withvalue != Xno && LDFLAGS="$LDFLAGS -L$withval"
	]
)
AC_ARG_WITH(gsl,
	AC_HELP_STRING(
		[--with-gsl=dir],
		[gsl prefix]
	),
	[
		test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -I$withval/include"
		test X$withvalue != Xno && LDFLAGS="$LDFLAGS -L$withval/lib"
	]
)
AC_ARG_WITH(gsl-includedir,
	AC_HELP_STRING(
		[--with-gsl-includedir=dir],
		[gsl include files in dir]
	),
	[
		test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -I$withval"
	]
)
AC_ARG_WITH(gsl-libdir,
	AC_HELP_STRING(
		[--with-gsl-libdir=dir],
		[gsl library libgsl.a in dir]
	),
	[
		test X$withvalue != Xno && LDFLAGS="$LDFLAGS -L$withval"
	]
)
AC_ARG_WITH(hdf5,
    AC_HELP_STRING(
        [--with-hdf5=dir],
        [hdf5 include files in dir]
        ),
        [
            test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -I$withval/include"
            test X$withvalue != Xno && LDFLAGS="$LDFLAGS -L$withval/lib"
        ]
)
AC_ARG_WITH(hdf5-includedir,
    AC_HELP_STRING(
        [--with-hdf5-includedir=dir],
        [hdf5 include files in dir]
        ),
        [
            test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -I$withval"
        ]
)
AC_ARG_WITH(hdf5-libdir,
    AC_HELP_STRING(
        [--with-hdf5-libdir=dir],
        [hdf5 library in dir]
        ),
        [
            test X$withvalue != Xno && LDFLAGS="$LDFLAGS -L$withval"
            HAVE_HDF5_=true
        ]
)
AC_ARG_WITH(wake-test,
	AC_HELP_STRING(
		[--with-wake-test],
		[enable test as done on the BS thesis of S. Pauli]
	),
	[
		test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -DnoENABLE_WAKE_TESTS -DEnoNABLE_WAKE_DEBUG"
	]
)
AC_ARG_WITH(trilinos-includedir,
	AC_HELP_STRING(
		[--with-trilinos-includedir=dir],
		[trilinos include files in dir]
	),
	[
		test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -I$withval"
	]
)
AC_ARG_WITH(trilinos-libdir,
	AC_HELP_STRING(
		[--with-trilinos-libdir=dir],
		[trilinos library ml in dir]
	),
	[
		test X$withvalue != Xno && LDFLAGS="$LDFLAGS -L$withval"
        HAVE_TRILINOS_=true
	]
)
AC_ARG_WITH(libdir,
    AC_HELP_STRING(
        [--with-libdir=dir],
        [additional library in dir]
    ),
    [
        test X$withvalue != Xno && LDFLAGS="$LDFLAGS $withval"
    ]
)
AC_ARG_WITH(blas,
    AC_HELP_STRING(
        [--with-blas=lib],
        [blas library or path (i.e --with-blas=/usr/path/lib/libblas.a]
    ),
    [
        test X$withvalue != Xno && BLASLIB="$withval"
        HAVE_BLAS_=true
    ]
)
AC_ARG_WITH(lapack,
    AC_HELP_STRING(
        [--with-lapack=lib],
        [lapack library or path (i.e --with-blas=/usr/path/lib/liblapack.a]
    ),
    [
        test X$withvalue != Xno && LAPACKLIB="$withval"
        HAVE_LAPACK_=true
    ]
)
AC_ARG_ENABLE([ml-solver],
    AC_HELP_STRING([--enable-ml-solver],
        [enable ml solver use [default=no]]),
        [HAVE_MLS_=false,
         test X$withvalue != Xno && HAVE_MLS_=true],
        [HAVE_MLS_=false]
)
################################################################################
## CHECK LIBS
AC_LANG_PUSH([C++])
    if test x$HAVE_LIBS_ = xtrue ; then
     echo check libs
    fi
AC_LANG_POP
## END LIBS
################################################################################

################################################################################
## CHECK CLASSIC AAA
AC_LANG_PUSH([C++])
    if test x$HAVE_CLASSIC_ = xtrue ; then
        AC_CHECK_LIB(CLASSIC, [main],
                         [LIBS="$LIBS -lCLASSIC"],
                         [AC_MSG_ERROR([OPAL needs CLASSIC])])
    fi
AC_LANG_POP
## END CLASSIC
################################################################################

################################################################################
## CHECK IPPL
AC_LANG_PUSH([C++])
    if test x$HAVE_IPPL_ = xtrue ; then
        AC_CHECK_LIB(ippl, [main],
                         [LIBS="$LIBS -lippl -lm"],
                         [AC_MSG_ERROR([OPAL needs IPPL])])
    fi
AC_LANG_POP
## END IPPL
################################################################################

################################################################################
## CHECK_HDF5
AC_LANG_PUSH([C++])
        AC_CHECK_LIB(hdf5, [main],
                         [LIBS="$LIBS -lhdf5 -lm -lz"],
                         [AC_MSG_ERROR([OPAL needs HDF5])])
AC_LANG_POP
## END HDF5
################################################################################

################################################################################
## CHECK_H5PART
AC_LANG_PUSH([C++])
AC_CHECK_LIB(H5hut, [h5_open_file],
                         [LIBS="$LIBS -lH5hutC -lH5hut -lhdf5"],
                         [AC_MSG_ERROR([OPAL needs H5hut])])
AC_LANG_POP
## END H5PART
################################################################################

################################################################################
## CHECK GSL LIB
AC_LANG_PUSH([C++])
	HAVE_GSL_=true
        AC_CHECK_HEADERS(
            [gsl/gsl_rng.h], 
            [], 
            [HAVE_GSL_=false, 
             AC_MSG_ERROR(gsl_rng.h not found!)], 
            [])
        
        AC_CHECK_HEADERS(
            [gsl/gsl_histogram2d.h], 
            [], 
            [HAVE_GSL_=false, 
             AC_MSG_ERROR(gsl_histogram2d.h not found!)], 
            [])

        #gslcblas has to be in path so that gsl test passes
        AC_CHECK_LIB(gslcblas, [main],
                         [GSLCBLASLIB=-lgslcblas,
                          LIBS="$LIBS -lgslcblas"],
                         [HAVE_GSL_=false, 
                          AC_MSG_ERROR([OPAL needs GSL])])
        AC_SUBST([GSLCBLASLIB])
        AC_CHECK_LIB(gsl, [main],
                         [GSLIB=-lgsl],
                         [HAVE_GSL_=false, 
                          AC_MSG_ERROR([OPAL needs GSL])])
        AC_SUBST([GSLIB])

        if test x$HAVE_GSL_ = xtrue ; then
            LIBS="$LIBS $GSLIB"
        else
            AC_MSG_ERROR([Cannot enable GSL due to missing dependencies!])
        fi
AC_LANG_POP
## END GSL
################################################################################

################################################################################
## CHECK header files and libraries for ML-Solver
AC_LANG_PUSH([C++])
    if test x$HAVE_MLS_ = xtrue ; then

        #if test x$HAVE_BLAS_ = xtrue ; then
        #    AC_CHECK_LIB($BLASLIB, [main],
        #                 [],
        #                 [HAVE_MLS_=false, 
        #                  AC_MSG_ERROR([ML solver depends on BLAS])])
        #fi
       
        #if test x$HAVE_LAPACK_ = xtrue ; then
        #    AC_CHECK_LIB($LAPACKLIB, [main],
        #                 [],
        #                 [HAVE_MLS_=false, 
        #                  AC_MSG_ERROR([ML solver depends on LAPACK])])
        #fi

        # Check Trilinos Headers 
        AC_CHECK_HEADERS(
            [ml_include.h], 
            [], 
            [HAVE_MLS_=false, 
             AC_MSG_ERROR(ml_include.h not found!)], 
            [])
        AC_CHECK_HEADERS(
            [ml_MatrixFreePreconditioner.h], 
            [], 
            [HAVE_MLS_=false, 
             AC_MSG_ERROR(ml_MatrixFreePreconditioner.h not found!)], 
            [])
        AC_CHECK_HEADERS(
            [ml_MultiLevelPreconditioner.h], 
            [], 
            [HAVE_MLS_=false, 
             AC_MSG_ERROR(ml_MultiLevelPreconditioner.h not found!)],
            [])
        AC_CHECK_HEADERS(
            [ml_MultiLevelOperator.h], 
            [], 
            [HAVE_MLS_=false, 
             AC_MSG_ERROR(ml_MultiLevelOperator.h not found!)], 
            [])
        AC_CHECK_HEADERS(
            [ml_epetra_utils.h], 
            [], 
            [HAVE_MLS_=false, 
             AC_MSG_ERROR(ml_epetra_utils.h not found!)], 
            [])
        AC_CHECK_HEADERS(
            [Epetra_Map.h], 
            [], 
            [HAVE_MLS_=false, 
             AC_MSG_ERROR(Epetra_Map.h not found!)], 
            [])
        AC_CHECK_HEADERS(
            [Epetra_Vector.h], 
            [], 
            [HAVE_MLS_=false, 
             AC_MSG_ERROR(Epetra_Vector.h not found!)], 
            [])
        AC_CHECK_HEADERS(
            [Epetra_LinearProblem.h], 
            [], 
            [HAVE_MLS_=false, 
             AC_MSG_ERROR(Epetra_LinearProblem.h not found!)], 
            [])
        AC_CHECK_HEADERS(
            [Epetra_Operator.h], 
            [], 
            [HAVE_MLS_=false, 
             AC_MSG_ERROR(Epetra_Operator.h not found!)], 
            [])
        AC_CHECK_HEADERS(
            [AztecOO.h], 
            [], 
            [HAVE_MLS_=false, 
             AC_MSG_ERROR(AztecOO.h not found!)], 
            [])
        AC_CHECK_HEADERS(
            [Teuchos_CommandLineProcessor.hpp], 
            [], 
            [HAVE_MLS_=false, 
             AC_MSG_ERROR(Teuchos_CommandLineProcessor.hpp not found!)], 
            [])
        #needed?
        #ifdef HAVE_MPI #include "Epetra_MpiComm.h" #else #include "Epetra_SerialComm.h" #endif
        #AC_CHECK_HEADERS([Epetra_VbrMatrix.h], [], [HAVE_MLS_=false, AC_MSG_ERROR(Epetra_VbrMatrix.h not found!)], [])
        #AC_CHECK_HEADERS([Epetra_Time.h], [], [HAVE_MLS_=false, AC_MSG_ERROR(Epetra_Time.h not found!)], [])
        #AC_CHECK_HEADERS([Galeri_Maps.h], [], [HAVE_MLS_=false, AC_MSG_ERROR(Galeri_Maps.h not found!)], [])
        #AC_CHECK_HEADERS([Galeri_Utils.h], [], [HAVE_MLS_=false, AC_MSG_ERROR(Galeri_Utils.h not found!)], [])
        #AC_CHECK_HEADERS([Galeri_CrsMatrices.h], [], [HAVE_MLS_=false, AC_MSG_ERROR(Galeri_CrsMatrices.h not found!)], [])
        #AC_CHECK_HEADERS([Galeri_VbrMatrices.h], [], [HAVE_MLS_=false, AC_MSG_ERROR(Galeri_VbrMatrices.h not found!)], [])

        # Check Trilinos Libs
        AC_CHECK_LIB(ml, [main],
                         [MLLIB=-lml],
                         [HAVE_MLS_=false, 
                          AC_MSG_ERROR([ML solver needs Trilinos compiled with --enable-ml])])
        AC_SUBST([MLLIB])
        AC_CHECK_LIB(galeri, [main],
                         [GALERILIB=-lgaleri],
                         [HAVE_MLS_=false, 
                          AC_MSG_ERROR([ML solver needs Trilinos compiled with --enable-galeri])])
        AC_SUBST([GALERILIB])
        
        AC_CHECK_LIB(amesos, [main],
                         [AMESOSLIB=-lamesos],
                         [HAVE_MLS_=false, 
                          AC_MSG_ERROR([ML solver needs Trilinos compiled with --enable-amesos])])
        AC_SUBST([AMESOSLIB])
        
        AC_CHECK_LIB(ifpack, [main],
                         [IFPACKLIB=-lifpack],
                         [HAVE_MLS_=false, 
                          AC_MSG_ERROR([ML solver needs Trilinos compiled with --enable-ifpack])])
        AC_SUBST([IFPACKLIB])
        
        AC_CHECK_LIB(epetraext, [main],
                         [EPETRAEXTLIB=-lepetraext],
                         [HAVE_MLS_=false, 
                          AC_MSG_ERROR([ML solver needs Trilinos compiled with --enable-epetraext])])
        AC_SUBST([EPETRAEXTLIB])
        
        AC_CHECK_LIB(epetra, [main],
                         [EPETRALIB=-lepetra],
                         [HAVE_MLS_=false, 
                          AC_MSG_ERROR([ML solver needs Trilinos compiled with --enable-epetra])])
        AC_SUBST([EPETRALIB])
        
        AC_CHECK_LIB(teuchos, [main],
                         [TEUCHOSLIB=-lteuchos],
                         [HAVE_MLS_=false, 
                          AC_MSG_ERROR([ML solver needs Trilinos compiled with --enable-teuchos])])
        AC_SUBST([TEUCHOSLIB])
        
        AC_CHECK_LIB(triutils, [main],
                         [TRIUTILSLIB=-ltriutils],
                         [HAVE_MLS_=false, 
                          AC_MSG_ERROR([ML solver needs Trilinos compiled with --enable-triutils])])
        AC_SUBST([TRIUTILSLIB])
        
        AC_CHECK_LIB(aztecoo, [main],
                         [AZTECOOLIB=-laztecoo],
                         [HAVE_MLS_=false, 
                          AC_MSG_ERROR([ML solver needs Trilinos compiled with --enable-aztecoo])])
        AC_SUBST([AZTECOOLIB])
		
        # only enable ML solver if all checks are ok
        if test x$HAVE_MLS_ = xtrue ; then
            CPPFLAGS="$CPPFLAGS -DnoHAVE_ML_SOLVER"
            LIBS="$LIBS $MLLIB $GALERILIB $AMESOSLIB $IFPACKLIB $EPETRAEXTLIB $EPETRALIB $TEUCHOSLIB -l$BLASLIB -l$LAPACKLIB -lguide -lpthread $TRIUTILSLIB $AZTECOOLIB -lm $ADDLIBS"
        else
            AC_MSG_ERROR([Cannot enable ML solver due to missing dependencies!])
        fi

        #-L/opt/intel-mkl/mkl-10.0/lib/em64t
        #3rd last: ifcore? -L/opt/intel/intel-10.0/fce-10.0/lib
        #-lguide -lpthread -lgfortran? -lifcore?

    fi
AC_LANG_POP
## END ML SOLVER
################################################################################
AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
