# SYNOPSIS
#
# AX_SEARCH_HEADERS([pkgname], [header1 header2], [path1 path2])
#
# DESCRIPTION
#
# Search for headers in all paths until found, resetting autoconf
# cache variables in between modifying CPPFLAGS and calling AC_CHECK_HEADERS.
#
AC_DEFUN([AX_SEARCH_HEADERS], [
for ac_test_location in $3 
do
    dnl Save the current state
    ax_search_headers_save_CPPFLAGS=${CPPFLAGS}

    CPPFLAGS="$CPPFLAGS -I${ac_test_location}"

    AC_MSG_CHECKING([$1 for $2 in ${ac_test_location}])
    AS_ECHO()
    _AS_ECHO_LOG([CPPFLAGS="${CPPFLAGS}"])

    AC_CHECK_HEADERS([$2], [ac_lib_$1=yes], [ac_lib_$1=no])

    dnl We have found the location, leave the loop:
    if test "${ac_lib_$1}" = "yes"
    then
        break;
    fi

    dnl Restore the state to original in case of unsuccessful attempt
    CPPFLAGS=${ax_search_headers_save_CPPFLAGS}
    AX_RESET_HEADERS_CACHE([$2])
done
]) #AX_SEARCH_HEADERS
