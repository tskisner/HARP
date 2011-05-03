dnl
dnl Configure checks for toy formats
dnl
toy_format=yes
AC_ARG_WITH([toy-format], [AC_HELP_STRING([--without-toy-format], [disable toy formats used for testing])])
if test "$with_toy_format" = "no"; then
   toy_format=no
else
   AC_DEFINE(ENABLE_TOY,1,[Define if we are building the toy file formats])
fi
AM_CONDITIONAL([BUILD_TOY], [test "$toy_format" = "yes"])
