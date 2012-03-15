AC_INIT([MC4],[1.0.0],[andreas.adelmann@psi.ch])
AC_CONFIG_AUX_DIR(config)
#disable f77 tests
m4_defun([_LT_AC_LANG_F77_CONFIG], [:])
# library creation
AC_PROG_LIBTOOL
# standard macros
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([./config.h:./config.in])
AM_INIT_AUTOMAKE([mc4],[1.0.0])
AC_PROG_CC
AC_PROG_CC_STDC
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_INSTALL
AC_LANG(C++)

AC_PATH_PROG(svnversioncommand, svnversion)
if test "X$svnversioncommand" = "X" || test `$svnversioncommand -n '.'` = "exported"; then
AC_DEFINE_UNQUOTED(SVN_VERSION, ["`cat saved_svn_version.txt`"], [repository svn version])
else
AC_DEFINE_UNQUOTED(SVN_VERSION, ["`svnversion -n`"], [repository svn version])
`svnversion -n > saved_svn_version.txt`
fi


IPPLDEFS="-w -DIPPL_MPI\ 
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

CPPFLAGS="$CPPFLAGS -DUSENAMESPACE -DH5_USE_16_API -DDOUBLE_REAL -DMC4HALOFINDER -fpermissive -DFFTW3 -I../initializer -I./HaloFinder/r1 -DDEBUG -DID_32 -DPOSVEL_32 -DGRID_32"

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
                LIBS="-lm $LIBS"
                HAVE_IPPL_=true
        ]
)
AC_ARG_WITH(h5part-includedir,
	AC_HELP_STRING(
		[--with-h5part-includedir=dir],
		[h5part include files in dir]
	),
	[
		test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -I$withval"
	]
)
AC_ARG_WITH(h5part-libdir,
	AC_HELP_STRING(
		[--with-h5part-libdir=dir],
		[h5part library libH5PART.a in dir]
	),
	[
		test X$withvalue != Xno && LDFLAGS="$LDFLAGS -L$withval"
        HAVE_H5PART_=true
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
        HAVE_GSL_=true
	]
)
AC_ARG_WITH(fftw3-includedir,
	AC_HELP_STRING(
		[--with-fftw3-includedir=dir],
		[fftw3 include files in dir]
	),
	[
		test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -I$withval"
	]
)
AC_ARG_WITH(fftw3-libdir,
	AC_HELP_STRING(
		[--with-fftw3-libdir=dir],
		[fftw3 library libfftw3.a in dir]
	),
	[
		test X$withvalue != Xno && LDFLAGS="$LDFLAGS -L$withval"
        HAVE_FFTW3_=true
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
AC_ARG_WITH(libdir,
    AC_HELP_STRING(
        [--with-libdir=dir],
        [additional library in dir]
    ),
    [
        test X$withvalue != Xno && LDFLAGS="$LDFLAGS $withval"
    ]
)
AC_ARG_WITH(libs,
    AC_HELP_STRING(
        [--with-libs=lib],
        [additional library]
    ),
    [
        test X$withvalue != Xno && LDFLAGS="$LDFLAGS $withval"
        #test X$withvalue != Xno && LIBS="$LIBS $withval" && ADDLIBS="$withval"
    ]
)

################################################################################
## CHECK IPPL
AC_LANG_PUSH([C++])
    if test x$HAVE_IPPL_ = xtrue ; then
        AC_CHECK_LIB(ippl, [main],
                         [LIBS="$LIBS -lippl"],
                         [AC_MSG_ERROR([MC4 needs IPPL])])
    fi
AC_LANG_POP
## END IPPL
################################################################################

################################################################################
## CHECK_HDF5
AC_LANG_PUSH([C++])
    if test x$HAVE_HDF5_ = xtrue ; then
        AC_CHECK_LIB(hdf5, [main],
                         [HDF5LIB="-lhdf5"],
                         [AC_MSG_ERROR([MC4 needs HDF5])])
    fi
AC_LANG_POP
## END HDF5
################################################################################

################################################################################
## CHECK_H5PART
AC_LANG_PUSH([C++])
    if test x$HAVE_H5PART_ = xtrue ; then
	if test x$HAVE_HDF5_ = xtrue ; then
          AC_CHECK_LIB(pH5Part, [main],
                         [LIBS="$LIBS -lpH5Part $HDF5LIB"],
                         [AC_MSG_ERROR([MC4 needs H5Part])])
        else
	  AC_MSG_ERROR([MC4 needs HDF5 with H5Part enabled])
        fi
    fi
AC_LANG_POP
## END H5PART
################################################################################

################################################################################
## CHECK_FFTW3
AC_LANG_PUSH([C++])
    if test x$HAVE_FFTW3_ = xtrue ; then
        AC_CHECK_LIB(fftw3, [main],
                         [LIBS="$LIBS -lfftw3"],
                         [AC_MSG_ERROR([MC4 needs FFTW])])
    fi
AC_LANG_POP
## END FFTW3
################################################################################

################################################################################
## CHECK GSL LIB
AC_LANG_PUSH([C++])
    if test x$HAVE_GSL_ = xtrue ; then

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
                          AC_MSG_ERROR([MC4 needs GSL])])
        AC_SUBST([GSLCBLASLIB])
        AC_CHECK_LIB(gsl, [main],
                         [GSLIB=-lgsl],
                         [HAVE_GSL_=false, 
                          AC_MSG_ERROR([MC4 needs GSL])])
        AC_SUBST([GSLIB])

        if test x$HAVE_GSL_ = xtrue ; then
            LIBS="$LIBS $GSLIB"
        else
            AC_MSG_ERROR([Cannot enable GSL due to missing dependencies!])
        fi
    fi
AC_LANG_POP
## END GSL
################################################################################

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
