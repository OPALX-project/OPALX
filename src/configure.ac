AC_INIT([OPAL],[1.0],[andreas.adelmann@psi.ch])
AC_CONFIG_AUX_DIR(config)
# library creation
AC_PROG_LIBTOOL
# standard macros
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([./config.h:./config.in])
AM_INIT_AUTOMAKE([opal],[1.0])
AC_PROG_CC
AC_PROG_CC_STDC
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_INSTALL
AC_LANG(C++)

IPPLDEFS="-DIPPL_MPI\ 
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
 -DWITH_BRICK\
 -DnoCOMP_GNUOLD\
 -DIPPL_STRINGSTREAM\
 -Drestrict=__restrict__ -DNOCTAssert"

AC_ARG_WITH(ippl-includedir,
        AC_HELP_STRING(
                [--with-ippl-includedir=dir],
                [ippl include files in dir]
        ),
        [
                test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -I$withval -I$HDF5HOME/include -I$HOME/include ${IPPLDEFS}"
        ]
)
AC_ARG_WITH(ippl-libdir,
        AC_HELP_STRING(
                [--with-ippl-libdir=dir],
                [ippl library libippl.a in dir]
        ),
        [
                test X$withvalue != Xno && LDFLAGS="$LDFLAGS -L$withval -L$HOME/lib -L$HDF5HOME/lib"
		LIBS="-lCLASSIC -ldoom -lippl -lfft -lhdf5 -lpH5Part -lz -lfftw3 -lm $LIBS"
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
	]
)
AC_ARG_WITH(doom-includedir,
	AC_HELP_STRING(
		[--with-doom-includedir=dir],
		[doom include files in dir]
	),
	[
		test X$withvalue != Xno && CPPFLAGS="$CPPFLAGS -I$withval"
	]
)
AC_ARG_WITH(doom-libdir,
	AC_HELP_STRING(
		[--with-doom-libdir=dir],
		[doom library libdoom.a in dir]
	),
	[
		test X$withvalue != Xno && LDFLAGS="$LDFLAGS -L$withval"
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
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
