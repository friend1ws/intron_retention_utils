#! /usr/local/bin/perl

use strict;

my $n = 0;
open(IN, "refGene.merged.coding.exon.bed") || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    print $_ . "\t" . "coding.exon_" . $n . "\n";
    $n = $n + 1;
}
close(IN);

my $n = 0;
open(IN, "refGene.merged.coding.intron.bed") || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    print $_ . "\t" . "coding.intron_" . $n . "\n";
    $n = $n + 1;
}
close(IN);

my $n = 0;
open(IN, "refGene.merged.coding.5putr.bed") || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    print $_ . "\t" . "coding.5putr_" . $n . "\n";
    $n = $n + 1;
}
close(IN);

my $n = 0;
open(IN, "refGene.merged.coding.3putr.bed") || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    print $_ . "\t" . "coding.3putr_" . $n . "\n";
    $n = $n + 1;
}
close(IN);

my $n = 0;
open(IN, "refGene.merged.noncoding.exon.bed") || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    print $_ . "\t" . "noncoding.exon_" . $n . "\n";
    $n = $n + 1;
}
close(IN);

my $n = 0;
open(IN, "refGene.merged.noncoding.intron.bed") || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    print $_ . "\t" . "noncoding.intron_" . $n . "\n";
    $n = $n + 1;
}
close(IN);

my $n = 0;
open(IN, "refGene.merged.intergene.bed") || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;
    print $_ . "\t" . "intergene_" . $n . "\n";
    $n = $n + 1;
}
close(IN);




