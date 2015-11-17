#! /usr/local/bin/perl

use strict;

my $input1 = $ARGV[0];
my $input2 = $ARGV[1];

my %key2count1 = ();
open(IN, $input1) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);
    my $key = join("\t", @F[0 .. 5]);
    $key2count1{$key} = $F[6] . "\t" . $F[7];
}
close(IN);

my %key2count2 = ();
open(IN, $input2) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);
    my $key = join("\t", @F[0 .. 5]);
    $key2count2{$key} = $F[6] . "\t" . $F[7];
}
close(IN);


my @keys = (keys %key2count1, keys %key2count2);
my %count = ();
@keys = grep {!$count{$_}++} @keys;

foreach my $key (sort @keys) {

    next if (not exists $key2count1{$key});
    next if (not exists $key2count2{$key});

    my $count1 = $key2count1{$key};
    my $count2 = $key2count2{$key};
    my @counts1 = split("\t", $count1);
    my @counts2 = split("\t", $count2);

    next if ($counts1[0] < 5);
    next if ($counts2[0] < 5);
    next if ($counts1[1] / $counts1[0] < 0.05);
    next if ($counts2[1] / $counts2[0] > 0.02);
    next if ($counts1[1] < 3);

    print $key . "\t" . $count1 . "\t" . $count2 . "\n";
}

