#!/usr/bin/perl
use strict;
use warnings;
my $file=shift; ## a list containing hmm profiles with their optimized threshold
open(F,"$file");
my @data=<F>;
close F;
foreach my $line(@data)
{
        chomp $line;
        my @temp=split '\s+',$line;
        system("sed '11a\ GA    $temp[1]   0.0' $temp[0] > $temp[0].txt");
}
