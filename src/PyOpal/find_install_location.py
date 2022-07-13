import distutils.sysconfig
install_location = distutils.sysconfig.get_python_lib(plat_specific=True,standard_lib=False)
print(install_location)