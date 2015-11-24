#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];
open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);

    if ($F[5] eq "+") {
        print $F[0] . "\t" . ($F[1] - 4) . "\t" . ($F[2] + 5) . "\t" . join("\t", @F[3 .. 5]) . "\n";
    } else {
        print $F[0] . "\t" . ($F[1] - 5) . "\t" . ($F[2] + 4) . "\t" . join("\t", @F[3 .. 5]) . "\n";
    }
}
close(IN);



