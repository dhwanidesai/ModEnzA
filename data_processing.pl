#!/usr/bin/perl
use Parallel::ForkManager; ##Install the fork manager (CPAN module)
use HMMmodE;
my $max_processors=40; #Input maximum number of processeors to be used
my $fork= new Parallel::ForkManager($max_processors);
my $file=shift;
open(LIST,$file);
my @file=<LIST>;
close LIST;
foreach (@file)
    {
    $fork->start and next; # do the fork
    chomp $_;
    system "formatdb -p T -i $_ -o T";
    system "blastall -i $_ -p blastp -d $_ -o $_.blast";
    $fork->finish; # do the exit in the child process
    }
$fork->wait_all_children;
foreach (@file)
    {
    chomp $_;
    system "mcxdeblast --score=b --sort=a --bcut=5 $_.blast";
    system "mcxassemble -b $_.blast -q -r max --map";           # transform raw cooccurrence data to mcl matrix format
                                                                                    # -b is the base name of files to be processed and output
                                                                                    # -q Warning options, we turn these on if we expect the raw data to be free of repeats   

    system "mcl $_.blast.sym -scheme 7 -I 1.2 -te 5 -o $_.clus --append-log=yes";        # performing mcl clustering of blast.sym files
                                                                                    # -scheme use a preset resource scheme,high schemes result/
                                                                                    # /in more n more expensive computations that may possibly be more accurate,defualt is 4.
    system "rm $_.blast *.psq *.psd *.phr *.pin *.psi *.err *.map *.raw *.hdr";
    my $clusfile=$_.".clus";
    my $tabfile=$_.".blast.tab";
    my $clusterObj = HMMmodE->new(-file => $_, -clusType => 'MCL', -clusFile => $clusfile, -tabFile => $tabfile);
    my $ref_array=$clusterObj->nsCluster();
    my @ns=@$ref_array;
    $ns_hash{$_}=\@ns;
        foreach my $clu_no(@ns)
            {
                my $nameref=$clusterObj->{'sample_outfile'}{$clu_no};
            }
    }
END { warn time - $^T, " seconds elapsed\n" }
