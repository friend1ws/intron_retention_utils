#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];

my %edge2count = ();
open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);
    if ($F[20] eq "1") {
        $edge2count{join("\t", @F[12 .. 17])} = $edge2count{join("\t", @F[12 .. 17])} + 1;
    }
}    
close(IN);

foreach my $edge (sort keys %edge2count) {
    print $edge . "\t" . $edge2count{$edge} . "\n";
}

