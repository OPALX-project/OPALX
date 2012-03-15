AC_INIT([classic],[5.1.2],[andreas.adelmann@psi.ch])
AC_CONFIG_AUX_DIR(config)
#disable f77 tests
m4_defun([_LT_AC_LANG_F77_CONFIG], [:])
# library creation
AC_PROG_LIBTOOL
# standard macros
# AC_CONFIG_AUX_DIR(config)
# AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([./config.h:./config.in])
AM_INIT_AUTOMAKE([classic],[5.1.2])
AC_PROG_CC
AC_PROG_CC_STDC
AC_PROG_CXX
AC_PROG_RANLIB
AC_LANG(C++)
AC_CONFIG_FILES([Makefile])

IPPLDEFS="-DIPPL_MPI\
 -DMPICH_SKIP_MPICXX\
 -DIPPL_DEBUG\ 
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
 -Drestrict=__restrict__ -DNOCTAssert -DPARALLEL_IO -w"

OPALSTUFF="$OPAL_ROOT/src"

AC_ARG_WITH(ippl-includedir,
        AC_HELP_STRING(
                [--with-ippl-includedir=dir],
                [ippl include files in dir]
        ),
        [
                test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -I$withval ${IPPLDEFS} -I$OPALSTUFF"
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
AC_ARG_WITH(h5part-includedir,
	AC_HELP_STRING(
		[--with-h5part-includedir=dir],
		[h5part include files in dir]
	),
	[
		test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -DH5_USE_16_API -I$withval"
	]
)
AC_ARG_WITH(h5part-libdir,
	AC_HELP_STRING(
		[--with-h5part-libdir=dir],
		[h5part library libH5PART.a in dir]
	),
	[
		test X$withvalue != Xno && LDFLAGS="$LDFLAGS -L$withval"
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
        ]
)
AC_OUTPUT
