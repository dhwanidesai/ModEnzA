#!/usr/bin/perl 
use lib '/HOME_ann/BII/biiswatis/modenza_archive/enzyme-nov28-2012/BioPerl-1.6.901';
use lib '/HOME_ann/BII/biiswatis/modenza_archive/enzyme-nov28-2012/BioPerl-Run-1.006900/lib/';  
# script 2 process the Swissprot sequence Fasta file, parse the Enzyme.dat flatfile and retrieve sequences associated with each EC number 
use Bio::DB::Fasta;
use Bio::SeqIO;

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
 
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
 
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

$USAGE = q/
USAGE:
perl 1_prepare-enzyme-files.pl       OPTIONS:

				     REQUIRED:		
					       -p, --spf <filename> this is the swissprot fasta file....uniprot_sprot.fasta or some such
                                               -e, --enz <filename> This should be th enzyme.dat file...
                                               
                                               
                                     OPTIONAL:
				               
				               -c <filename> copy the processed Swissprot file with the given filename; the program pre-processes 
						             the Swissprot file to modify the Fasta description lines, so if a filename is   
				                             specified, the copy will be saved with given name; By default, the program saves the 
				                             copy as "modified-swissprot.fasta"; if you don't want to make a copy, use the -x 
				                             option
                                               -x <T|F>      Process the Swissprot fasta file without making a copy; this will modify the Fasta 
                                                             description lines of original Swissprot Fasta file that you downloaded; see the -c 
							     option; if none of -c or -x is used, a copy is saved as "modified-swissprot.fasta" 
							     by default;  if -c and -x are both used, then -x is over-ridden and a copy will be 
							     saved by the filename specified by -c; 
                                               -o, --out <path\/to\/output\/directory> 
				                             The path to the directory (including the directory name) where the enzyme 
				                             files will be created; Enzyme files with less than 3 sequences can be moved to a 
				                             separate directory to start Tier2 ModenzA profiles (see the -s option);
				                             By default it is the present directory (the 
				                             directory from where the script is run)
                                               -s, --split <directory name> 
					                     If a directory name is specified, the program will move the enzyme files containing 
					                     less than 3 sequences to this directory; these sequences can then be used as 
					                     queries for constructing Tier2 ModEnzA profiles
                                        


EXAMPLES:

## Write out all sequences in the current directory.....
## the sequences will be stored in separate files for each EC number in the enzyme database --- e.g. 1.1.1.1.seq etc...

perl 1_prepare-enzyme-files.pl -p uniprot_sprot.fasta -e enzyme.dat 


## Write out files with less than 3 sequences to a separate directory: make sure that the directory exists 
  
perl 1_prepare-enzyme-files.pl --split tier2seq -p uniprot_sprot.fasta -e enzyme.dat 


## write out the tier1 and tier2 enzyme files to the specified directories

perl 1_prepare-enzyme-files.pl -p uniprot_sprot-19-jan-2010-test.fasta -e enzyme-19-jan-2010.dat -o t1directory -s t2directory


AUTHOR:
 
Dhwani Desai
 
/;

use Getopt::Long;

GetOptions (
'p|spf=s' => \$file,
'e|enz=s' => \$idfile,
'c=s' => \$copyname,
'x=s' => \$delete,
's|split=s' => \$t2directory,
'o|out=s' => \$outdir
) or die $USAGE;

die $USAGE if !$file or !$idfile;

if ($outdir)
{
check_directory($outdir);
}
if ($t2directory)
{
check_directory($t2directory);
}


#### modify the Swissprot Fasta file description/headers

open (F,$file);
chomp (@sp=<F>);
close F;

if ($copyname)
{
$file=$copyname;
}
elsif ($delete) 
{
$file=$file;
}
else
{
$file="modified-swissprot.fasta";
}

open (OUT,">$file");

foreach $line(@sp)
{
chomp $line;

	if ($line=~/^>/)
	{
	$line =~ s/^>sp\|//g;
	$line =~ s/\|/ /g;
	print OUT ">$line\n";
	}
	else
	{
	print OUT "$line\n";
	}
}
close OUT;

my $db = Bio::DB::Fasta->new($file);  # one file or many files.....create an index of the db file
open (LIST,">list_of_all_training_seq");
open (FH, $idfile); #parse the enzyme.dat file to get ids for each EC no.
chomp (@file=<FH>);
for ($line=0;$line<=$#file;$line++)
{
chomp($file[$line]);
	if ($file[$line] =~ /^ID/)
        {
        $file[$line] =~ s/ID//;
        $file[$line] =~ s/\s+//;
        
        $enzyme=$file[$line];
        print " Enzyme:$enzyme\n"; #get the EC no
        
        $outfile="$enzyme.seq";    #the filename to write out the sequences

        	while ($file[$line] !~ /\/\//)
        	{
        	$line++;
        	next if $file[$line] !~ /^DR/;
			if ($file[$line] =~ /^DR/)
			{
			$file[$line]=~s/^DR//;
			$file[$line] =~ s/\s+//;
			#print "$file[$line]\n";
			@seq = (@seq,split(/;/,$file[$line])); #Getting the acc. ids from the DR field
			}
        	
        	}
	#print "@seq\n";
		if (scalar @seq < 3 && scalar @seq >= 1 && $t2directory)
		{
		$out = Bio::SeqIO->new('-file' => ">$t2directory/$outfile",'-format' => 'Fasta'); # open the output file for writing
		}
		elsif (scalar @seq < 3 && scalar @seq >=1 && !$t2directory)
		{
		#print "no separate diretory for TierII given; using current dir\n";
		$out = Bio::SeqIO->new('-file' => ">$outfile",'-format' => 'Fasta'); # open the output file for writing
		}
		elsif (scalar @seq >= 3 && $outdir)
		{
		$out = Bio::SeqIO->new('-file' => ">$outdir/$outfile",'-format' => 'Fasta'); # open the output file for writing
		}
		elsif (scalar @seq >= 3 && !$outdir)
		{
		#print "no separate diretory for TierI given; using current dir\n";
		$out = Bio::SeqIO->new('-file' => ">$outfile",'-format' => 'Fasta'); # open the output file for writing
		print LIST "$outfile\n";
		}

		
		
		foreach $id (@seq)
		{
		@temp=split(/,/,$id);
		$accno=$temp[0];
		$accno=~s/\s+//g;
		$temp[1]=~s/\s+//g;
		#print "$accno\n";
		$id=$accno."|".$temp[1];

		$desc=$db->header($accno); # Access the seq description using Bio::DB::Fasta obj
		#$desc=~s/\s+//;
		
		#print "description $desc\n";
			if ($desc =~ /ragment/) # Check for fragments
			{
			print "$id is a fragment\n";
			$seqobj=();
			}
			else
			{
			$seqobj = $db->get_Seq_by_id($accno); # get a PrimarySeq obj
			}
			
		next if(!$seqobj);		
		$out->write_seq($seqobj); # write out the sequence to $outfile

		}
	
	@seq=();      
        }
	

}



#### subroutines #####

sub check_directory
{
my $directory=shift;
# checking if the test directory already exists
        if (-d $directory)
        { # the -d bit means check if there is a directory called ...
        print "$directory already exists, sequence files will be saved in $directory\n";
        return 1;
        }
        else
        {
        die "$directory does not exist.....make sure $directory exists and the path is correct.....\n\n";
        #mkdir ($directory, 0777) || die "sorry system is unable to create output directory $directory";
        }


}
END { warn time - $^T, " seconds elapsed\n" }
