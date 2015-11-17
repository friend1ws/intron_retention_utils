#! /usr/local/bin/perl

use strict;

my %junction2annot = ();
keys(%junction2annot) = 10000000;


print STDERR "Start reading refGene.txt\n";
my $n = 0;
open(IN, "refGene.txt") || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);

    my $chr = $F[2]; 
    my @starts = split(",", $F[9]);
    my @ends = split(",", $F[10]);
    my $strand = $F[3];
    my $exonNum = $F[8];
    my $gene = $F[1];
    my $symbol = $F[12];

    if ($symbol eq "HIST2H4A") {
        $DB::single = 1;
    }

    for (my $i = 0; $i <= $#starts; $i++) {

        my $annot = "";

        my $key = $chr . "\t" . $starts[$i] . "\t" . ($starts[$i] + 1) . "\t" . "-";
        if ($strand eq "+") {
            $annot = $symbol . "(" . $gene . ")" . "." . $i . ".start";
        } else {
            $annot = $symbol . "(" . $gene . ")" . "." . ($exonNum - $i) . ".end";
        }

        if (exists $junction2annot{$key}) {
            $junction2annot{$key} = $junction2annot{$key} . "," . $annot;
        } else {
            $junction2annot{$key} = $annot;
        }
    }

    for (my $i = 0; $i <= $#ends; $i++) {
        
        my $annot = ""; 
        my $key = $chr . "\t" . ($ends[$i] - 1) . "\t" . $ends[$i] . "\t" . "+";
        if ($strand eq "+") {
            $annot = $symbol . "(" . $gene . ")" . "." . $i . ".end";
        } else {
            $annot = $symbol . "(" . $gene . ")" . "." . ($exonNum - $i - 1) . ".start";
        }
        
        if (exists $junction2annot{$key}) {
            $junction2annot{$key} = $junction2annot{$key} . "," . $annot;
        } else {
            $junction2annot{$key} = $annot;
        }
    }

    $n = $n + 1;
    if ($n % 1000 == 0) {
        print STDERR "$n genes completed.\n";
    }

}
close(IN);
print STDERR "Reading refGene.txt completed.\n";

print STDERR "Start reading endGene.txt.\n";
my $n = 0;
open(IN, "ensGene.txt") || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);
    
    my @starts = split(",", $F[9]);
    my @ends = split(",", $F[10]);
    my $chr = $F[2];
    my @starts = split(",", $F[9]);
    my @ends = split(",", $F[10]);
    my $strand = $F[3];
    my $exonNum = $F[8];
    my $gene = $F[1];

    for (my $i = 0; $i <= $#starts; $i++) {

        my $annot = "";
        my $key = $chr . "\t" . $starts[$i] . "\t" . ($starts[$i] + 1) . "\t" . "-";
        if ($strand eq "+") {
            $annot = $gene . "." . $i . ".start";
        } else {
           $annot = $gene . "." . ($exonNum - $i) . ".end";
        }

        if (exists $junction2annot{$key}) {
            $junction2annot{$key} = $junction2annot{$key} . "," . $annot;
        } else {
            $junction2annot{$key} = $annot;
        }
    }

    for (my $i = 0; $i <= $#ends; $i++) {

        my $annot = "";
        my $dir = $strand eq "+" ? "-" : "+";
        my $key = $chr . "\t" . ($ends[$i] - 1) . "\t" . $ends[$i] . "\t" . "+";
        if ($strand eq "+") {
            $annot = $gene . "." . $i . ".end";
        } else {
            $annot = $gene . "." . ($exonNum - $i - 1) . ".start";
        }

        if (exists $junction2annot{$key}) {
            $junction2annot{$key} = $junction2annot{$key} . "," . $annot;
        } else {
            $junction2annot{$key} = $annot;
        }
    }

    $n = $n + 1;
    if ($n % 1000 == 0) {
        print STDERR "$n genes completed.\n";
    }

}
close(IN);   

print STDERR "Reading ensGene.txt completed.\n";

print STDERR "Start reading knownGene.txt.\n";

my $n = 0;
open(IN, "knownGene.txt") || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    my @F = split("\t", $_);

    my @starts = split(",", $F[8]);
    my @ends = split(",", $F[9]);
    my $chr = $F[1];
    my $strand = $F[2];
    my $exonNum = $F[7];
    my $gene = $F[0];

    for (my $i = 0; $i <= $#starts; $i++) {

        my $annot = "";
        my $dir = $strand eq "+" ? "-" : "+";
        my $key = $chr . "\t" . $starts[$i] . "\t" . ($starts[$i] + 1) . "\t" . "-";
        if ($strand eq "+") {
            $annot = $gene . "." . $i . ".start";
        } else {
           $annot = $gene . "." . ($exonNum - $i) . ".end";
        }

        if (exists $junction2annot{$key}) {
            $junction2annot{$key} = $junction2annot{$key} . "," . $annot;
        } else {
            $junction2annot{$key} = $annot;
        }
    }

    for (my $i = 0; $i <= $#ends; $i++) {

        my $annot = "";
        my $dir = $strand eq "+" ? "-" : "+";
        my $key = $chr . "\t" . ($ends[$i] - 1) . "\t" . $ends[$i] . "\t" . "+";
        if ($strand eq "+") {
            $annot = $gene . "." . $i . ".end";
        } else {
            $annot = $gene . "." . ($exonNum - $i - 1) . ".start";
        }

        if (exists $junction2annot{$key}) {
            $junction2annot{$key} = $junction2annot{$key} . "," . $annot;
        } else {
            $junction2annot{$key} = $annot;
        }
    }
          
    $n = $n + 1;
    if ($n % 1000 == 0) {
        print STDERR "$n genes completed.\n";
    }
 
}
close(IN);

print STDERR "Reading knownGene.txt completed.\n";

foreach my $junction (sort keys %junction2annot) {
    my @temp = split("\t", $junction);
    next if (not $temp[0] =~ /^chr[\dXYM]+$/);
    print join("\t", @temp[0 .. 2]) . "\t" . $junction2annot{$junction} . "\t" . "0" . "\t" . $temp[3] . "\n";
}


