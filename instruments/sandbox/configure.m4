dnl
dnl Configure checks for sandbox formats
dnl
sandbox_format=yes
AC_ARG_WITH([sandbox-format], [AC_HELP_STRING([--without-sandbox-format], [disable sandbox formats used for testing])])
if test "$with_sandbox_format" = "no"; then
   sandbox_format=no
else
   AC_DEFINE(ENABLE_SANDBOX,1,[Define if we are building the sandbox file formats])
fi
AM_CONDITIONAL([BUILD_SANDBOX], [test "$sandbox_format" = "yes"])
