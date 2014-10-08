#!/usr/bin/env perl

use strict;

# get git revision and update conf.py

my $gitver = `git rev-parse --short HEAD`;
chomp($gitver);

my $docver = "";
open ( CONF, "./sphinx/source/conf.py" ) || die;
while ( <CONF> ) {
	if ( /^version = .*/ ) {
		my @F = split;
		$docver = $F[2];
	}
}
close ( CONF );

system ( "perl -i -p -e 's/^release =.*/release = \"$docver ($gitver)\"/;' ./sphinx/source/conf.py" );

