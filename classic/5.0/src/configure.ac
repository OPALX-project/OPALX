AC_INIT([classic],[5.0],[andreas.adelmann@psi.ch])
AC_CONFIG_AUX_DIR(config)
# library creation
AC_PROG_LIBTOOL
# standard macros
# AC_CONFIG_AUX_DIR(config)
# AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([./config.h:./config.in])
AM_INIT_AUTOMAKE([classic],[5.0])
AC_PROG_CC
AC_PROG_CC_STDC
AC_PROG_CXX
AC_PROG_RANLIB
AC_LANG(C++)
AC_CONFIG_FILES([Makefile])

IPPLDEFS="-DIPPL_MPI\
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
 -DWITH_BRICK\
 -DnoCOMP_GNUOLD\
 -DIPPL_STRINGSTREAM\
 -Drestrict=__restrict__ -DNOCTAssert -DPARALLEL_IO -I$H5Part/src"

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
        ]
)
AC_OUTPUT
