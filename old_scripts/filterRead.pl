#! /usr/local/bin/perl
#$ -S /usr/local/bin/perl
#$ -cwd

use strict;

my $input = $ARGV[0];
my $thres = $ARGV[1];

open(IN, $input) || die "cannot open $!";
while(<IN>) {
    s/[\r\n\"]//g;  
    my @F = split("\t", $_);

    if ($F[0] =~ /^\@/) {
        print join("\t", @F) . "\n";
    }
    next if ($F[4] < $thres);
    next if ($F[6] ne "=");

    $DB::single = 1;
    my @flags = reverse split("", sprintf("%011b", $F[1]));

    next if ($flags[1] == 0);
    
    if ( ($flags[4] == 0 and $flags[5] == 1) and ($F[3] <= $F[7]) ) {
        print join("\t", @F) . "\n";
    }    

    if ( ($flags[4] == 1 and $flags[5] == 0) and ($F[3] >= $F[7]) ) {
        print join("\t", @F) . "\n";
    }       


}
close(IN);


        
