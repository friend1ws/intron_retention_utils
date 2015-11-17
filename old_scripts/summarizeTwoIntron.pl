#! /usr/local/bin/perl

use strict;

my $input_ref1 = $ARGV[0];
my $input_ref2 = $ARGV[1];
my $input_junc = $ARGV[2];

my %key2count_ref1 = ();
open(IN, $input_ref1) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);
    $key2count_ref1{join("\t", @F[0 .. 5])} = $F[6];
}
close(IN);

my %key2count_ref2 = ();
open(IN, $input_ref2) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);
    $key2count_ref2{join("\t", @F[0 .. 5])} = $F[6];
}   
close(IN);


open(IN, $input_junc) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);
    my $key = join("\t", @F[0 .. 5]);

    my $count_ref1 = exists $key2count_ref1{$key} ? $key2count_ref1{$key} : 0;
    my $count_ref2 = exists $key2count_ref2{$key} ? $key2count_ref2{$key} : 0;

    print $key . "\t" . $count_ref1 . "\t" . $count_ref2 . "\t".  $F[6] . "\t" .  $F[7] . "\n";
}
close(IN);


    
