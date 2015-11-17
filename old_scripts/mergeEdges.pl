#! /usr/local/bin/perl

use strict;
my $input1 = $ARGV[0];
my $input2 = $ARGV[1];

my %edge2count = ();
open(IN, $input2) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);
    $edge2count{join("\t", @F[0 .. 5])} = $F[6];
}
close(IN);

open(IN, $input1) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);
    my $key = join("\t", @F[0 .. 5]);
    my $count = exists $edge2count{$key} ? $edge2count{$key} : 0;
    print join("\t", @F[0 .. 6]) . "\t" . $count . "\n";
}
close(IN); 
 
