#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];
my $input_ref = $ARGV[1];

my %key2info = ();
open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);

    my $key = join("\t", @F[0 .. 5]);
    $key2info{$key} = $F[6] . "," . $F[7] . "\t" . $F[8] . "," . $F[9];
}
close(IN);

open(IN, $input_ref) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);
    my $key = join("\t", @F[0 .. 5]);

    if (exists $key2info{$key}) {
        print $key . "\t" . $key2info{$key} . "\t" . join("\t", @F[6 .. $#F]) . "\n";
    }
}
close(IN);
