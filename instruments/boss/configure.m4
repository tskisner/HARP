dnl
dnl Configure checks for boss formats
dnl
boss_format=yes
AC_ARG_WITH([boss-format], [AC_HELP_STRING([--without-boss-format], [disable boss formats used for testing])])
if test "$with_boss_format" = "no"; then
   boss_format=no
else
   AC_DEFINE(ENABLE_BOSS,1,[Define if we are building the boss file formats])
fi
AM_CONDITIONAL([BUILD_BOSS], [test "$boss_format" = "yes"])
