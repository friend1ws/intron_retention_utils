#! /usr/local/bin/perl

use strict;

my $input = $ARGV[0];

my %edge2count = ();
my $tempRead = "";
my %tempKeys = ();
open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);
    $F[3] =~ s/\/[12]$//;

    if ($F[3] ne $tempRead) {
        for my $key (keys %tempKeys) {
            $edge2count{$key} = $edge2count{$key} + 1;
        }
        $tempRead = $F[3];
        %tempKeys = ();        
    }

    if ($F[18] eq "10") {
        if ($F[17] eq "+") {
            $F[13] = $F[13] + 4;
            $F[14] = $F[14] - 5;
        } else {
            $F[13] = $F[13] + 5;
            $F[14] = $F[14] - 4;
        }
        $tempKeys{join("\t", @F[12 .. 17])} = 1;
    }
}    
close(IN);

for my $key (keys %tempKeys) {
    $edge2count{$key} = $edge2count{$key} + 1;
}


foreach my $edge (sort keys %edge2count) {
    print $edge . "\t" . $edge2count{$edge} . "\n";
}

