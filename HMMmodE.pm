#----------------------------------------------------------------------------
# PACKAGE : HMMmodE
# PURPOSE : To encapsulate code for parsing MCL cluster output and generating HMMer profiles of each cluster
# AUTHOR  : Dhwani Desai (dhwanidesai@gmail.com) and Soumyadeep Nandi (soumyadeep_nandi@yahoo.com) under the guidence of Dr. Andrew M.Lynn
# REVISION: Revised to work with HMMER 3 by Swati Sinha (swati.scisjnu@gmail.com) under the guidence of Dr. Andrew M.Lynn
# VERSION: HMMmodE 3
# CHange Log - modified generate samples.....Andrew 6/7/2014		
#----------------------------------------------------------------------------
BEGIN {$ENV{MUSCLEDIR} = '/HOME_ann/sw/usr/bin/'; }
use lib '/HOME_ann/BII/biiswatis/modenza_archive/enzyme-nov28-2012/BioPerl-1.6.901/';
use lib '/HOME_ann/BII/biiswatis/modenza_archive/enzyme-nov28-2012/BioPerl-Run-1.006900/lib/';
package HMMmodE;
use Bio::DB::Fasta;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::SeqIO;
use Bio::SearchIO;
use Data::Dumper;
use POSIX;
use Statistics::Lite qw(:all);
use Hmm3;

## POD Documentation:

=head2 new

Title   : new
Usage   : $clusterObj = HMMmodE->new(-file => $seqfilename, -clusType => 'MCL', -clusFile => '/path/to/clusterfile', -tabFile => '/path/to/MCL .tab file')
Function: Creates a new "MCL cluster file" object
Returns : An HMMmodE stream initialised with the appropriate format
Args    : Named parameters:
             -file => $filename 		 #Compulsory argument - file containing all the sequences of a given family 
             -clusType => 'MCL' 		 #Compulsory argument - For methods other than MCL, -clusType should be 'o' and the clusters should be formatted as described later;
	     -clusFile => '/path/to/clusterfile' #Compulsory argument - The file containing the clusters	
	     -tabFile => '/path/to/MCL .tab file'#Optional - if -clusType is 'MCL', specify the path to the .tab file 
=cut
sub new {
	my ($class, %param) = @_;
	my $self = bless {}, $class;	
	my $file = $param{-file}||die "Sequence file not passed;quitting...\n";
	$self->{'file'} = $file;
	my $db = new Bio::DB::Fasta($file);  # one file or many files
	my $methref=_get_method_ref($db,'get_Seq_by_id');
	my $clus_type=$param{-clusType};
	if (!$clus_type){die "Clustering method not specified; Set -clusType => 'MCL' or 'o'\n";}
	my $clus_file=$param{-clusFile}||die "Cluster file not passed;quitting...\n";	
	$self->{'clusfile'} = $clus_file;
	my $tab_file=$param{-tabFile};
	
# Parsing the MCL output file file using some temporary variables

my @node_line=();
my $nodes;
my @cluster_line=();
my $clusters;
my @file=();
my $mclheader;
my $line;
my @tab=();
my %tab_hash=();
my @id_nos=();
my $clus_no;
my $temp_clus;
my $temp_ids;
my @clus=();
my %clus_hash=();
my %seq_hash=();
my %aln_hash=();
my $flag=0;
my %sample_out_hash=();
if ($clus_type eq 'MCL')
{
	if (-s $tab_file)
	{
	open (TAB,$tab_file);
		while (<TAB>)
		{
		chomp $_;
	#	print "$_\n";
		next if (/^#/);
		@tab = split (/\s+/, $_);
		$tab_hash{$tab[0]}=$tab[1];
		}
	close (TAB);
	open(FH,"$clus_file")|| die "can not open $clus_file";
	print $clus_file;
	chomp (@file=<FH>);
	if ($file[0]=~/\#/){shift @file;};
	$mclheader=$file[0];
	if ($mclheader !~ /\(mclheader/){die "$self->{'clusfile'} not a valid MCL cluster file\n";}
		for ($line=0;$line<=$#file;$line++)
		{
		chomp $file[$line];					
			if ($file[$line]=~/^number of nodes/)
			{
			@node_line=split (" ", $file[$line]);
			$nodes=$node_line[$#node_line];
			$nodes=~s/[<>]//g;
			$line++;
			}
			if ($file[$line]=~/^number of clusters/)
			{
			@cluster_line=split (" ", $file[$line]);
			$clusters=$cluster_line[$#cluster_line];
			$clusters=~s/[<>]//g;
			$line++;
			}
				
			next if $file[$line]=~/^(100){1,20}/;	
 			
 			if ($file[$line] =~ /^[0-9]+/)
 			{
 			@id_nos=split (" ", $file[$line]);
 			$clus_no=shift (@id_nos);
 				if ($file[$line] !~ /\$$/)
 				{
 				$line++;
 				}
 				else
 				{
 				$flag=1;
 				pop @id_nos;
 				}
 				
 				
				while ($file[$line] !~ /\$$/)
				{
				push (@id_nos,split (" ",$file[$line]));
				$line++;
				}
				
				if (!$flag)
				{
				push (@id_nos,split (" ",$file[$line]));
				pop @id_nos;
				}		
 			}
 			
 			else
 			{
 			next;
 			}
 			
 			my @ids=();
 			@ids=map {$tab_hash{$_}} @id_nos;
 			$clus_hash{($clus_no+1)}=[@ids];
 			$self->{'cluster_ids'} = \%clus_hash;	
			my @seqobjs=();
			@seqobjs=map {&$methref("$_")} @ids;			
			$seq_hash{($clus_no+1)}=[@seqobjs];
			$self->{'cluster_seq'} = \%seq_hash;
			my $temp_clus=$clus_no+1;	
			my $outfile=$file;
			$outfile=~s/\.seq//;
			$outfile = "$outfile"."_$temp_clus";
			$sample_out_hash{$temp_clus}=$outfile;
			
		}
	$self->{'sample_outfile'}=\%sample_out_hash;
	close FH;
	}
	
	else 
	{
	die "Cluster type is MCL, but the .tab file is missing; Can't proceed; Try -clusType => 'o'\n";
	}
}

else
{
my @ids=();
my $outfile=$file;
	if (-s $clus_file)
	{
	open(FH,"$clus_file"); 
			
		while ($line=<FH>)
		{
		
			if ($. == 1 && $line !~/Cluster/)
			{
			die "$clus_file does not match the prescribed format for -clusType => 'o'\n ";
			}
			
			if ($line =~ /Cluster/)
			{
			my $outfile=$file;
			$outfile=~s/\.seq//;
			$line=~s/[.,-]//g;
			@id_nos=split(" ",$line);
			$clus_no=$id_nos[1];
				if ($temp_clus)
				{
				$outfile = "$outfile"."_$temp_clus";
				#print "out $outfile\n"; 
				$sample_out_hash{$temp_clus}=$outfile;
				$clus_hash{$temp_clus}=[@ids];
				my @seqobjs=();
				@seqobjs=map {&$methref("$_")} @ids;
				$seq_hash{$temp_clus}=[@seqobjs];
				$temp_ids=scalar (@ids);
				$nodes += $temp_ids;
				}
			@ids=();
			$temp_clus=$clus_no;
			}
			else
			{
			$line =~ s/\s+//;
				if ($line ne "")
				{	
				push (@ids,$line);
				}
			}
		}
	$outfile = "$outfile"."_$temp_clus";
	$sample_out_hash{$temp_clus}=$outfile;	
	$clus_hash{$temp_clus}=[@ids];
	$temp_ids=scalar (@ids);
	$nodes += $temp_ids;
	@clus = keys %clus_hash;
	$clusters=scalar @clus;
	$self->{'cluster_ids'} = \%clus_hash;
	
	my @seqobjs=();
	@seqobjs=map {&$methref("$_")} @ids;
	$seq_hash{$temp_clus}=[@seqobjs];
	$self->{'cluster_seq'} = \%seq_hash;
					
				
	}
	
	else
	{
	die "Cluster type is 'o', but the cluster file is missing; Can't proceed; \n";
	}
}

# End of parsing the MCL file

$self->{'nodes'} = $nodes;
$self->{'clusters'} = $clusters;


return $self;
}


# =head2 destroy
# 
#  Title   : destroy
#  Usage   : $HMMmodEObj->destroy()
#  Function: Release the memory
# 
# =cut
# 
# 
# sub destroy 
# {
# my $self = shift;
# printf ("Cleaning object \"$self\" from memory at %s\n", scalar localtime);
# }


=head2 cleanUp

 Title   : cleanUp
 Usage   : $HMMmodEObj->cleanUp();
 Function: Delete the temp aln and hmm files in the train_samples directory;

=cut

sub cleanUp
{
my $self=shift;
my $path_train=$self->{'path_train'};
my $path_test=$self->{'path_test'};
system("rm $path_test/*");
system ("rm $path_train/*.aln");
system ("rm $path_train/*prof.hmm");
system ("rm $path_train/*.index");
system ("rm $path_train/*.seq");
system ("rm $path_train/*_?.??.hmm");
system ("rm $path_train/*_?.?.hmm");
system ("rm $path_train/neg.seq.*");
system ("rm $path_train/tot_fp.hmm");
}


# =head2 clusInfo

#  Title   : clusInfo
#  Usage   : $arrayref=$HMMmodEObj->clusInfo(Cluster number 1 upto $clusters);
#  Function: Fetch the sequence ids for a cluster number
#  Returns : reference to an array of ids; when called without arguments- returns a hash of cluster-number and no. of ids in each cluster
#  Args    : Cluster number or id (1 up to the total number of clusters)
 
# =cut
  
sub clusInfo 
{
my ($self,$clus_no) = @_;
	if($clus_no)
	{
		return $self->{'cluster_ids'}{$clus_no};
	}
	else
	{
		my $arrayref;
		my @ids=();
		my $id_no;
		my %info_hash;
		my $clusters=$self->{'clusters'};
		for (my $a=1; $a<=$clusters;$a++)
		{
			$arrayref=$self->{'cluster_ids'}{$a};
			@ids=@$arrayref;
			$id_no=scalar @ids;
			$info_hash{$a}=$id_no;
		}
	return \%info_hash;
	}
}

# =head2 writeSeq
 
#  Title   : writeSeq
#  Usage   : $HMMmodEObj->writeSeq(Cluster number 1 upto $clusters);
#  Function: Fetch the ids for cluster number, fetch sequences for the ids and write out the sequence file using BioPerl sequence objects 
#  Returns : ref_to_array of sequence objects
#  Args    : Cluster number or id (1 up to the total number of clusters) and the output file for writing the sequences to; If called without the outfile, it just returns a ref_to_array of sequence objects....
 
# =cut

 
sub writeSeq 
{
my ($self,$clus_no,$outfile,$combine_flag) = @_;
my @seqobjs=();
my $ref = $self->{'cluster_seq'}{$clus_no};
my @seqobjs = @$ref;

	if ($outfile)
	{
		if ($combine_flag eq "c")
		{
		my $out = Bio::SeqIO->new('-file' => ">>$outfile", '-format' => 'Fasta');
		$out->write_seq(@seqobjs);
		}
	
		else
		{
		my $out = Bio::SeqIO->new('-file' => ">$outfile", '-format' => 'Fasta');
		$out->write_seq(@seqobjs);
		}
	}
	else
	{
	die "Outfile not provided; can't write sequence to file\n";
	}

}

# =head2 mAlign
 
#  Title   : mAlign
#  Usage   : $HMMmodEObj->mAlign("alignment method", \@params, $clus_no-cluster number or id) #Currently works only for Muscle or Clustalw using Bio::Tools::Run::Alignment::Muscle or Bio::Tools::Run::Alignment::ClustalW from BioPerl
#  Function: Use Muscle to multiple-align a group of unaligned sequences
#  Returns : A Bio::SimpleAlign object
#  Args    : 	"alignment method" -> string indicating the alignment program (currently only Muscle or ClustalW)
#		\@params -> parameters for running Muscle (passed as a reference)
#		Cluster number or id for which alignment needs to be done

# =cut

sub mAlign 
{

my ($self,$method,$ref_params,$clus_no,$seq1,$seq2,$prof_flag) = @_;

my @params=@$ref_params;

	if (defined $seq1 && defined $seq2 && defined $prof_flag)
	{
	my $out=$self->{'sample_outfile'}{$clus_no};
	my $path_train=$self->{'path_train'};
		
		if ($method =~ /muscle/)
		{	
		system ("/HOME_ann/sw/usr/bin/muscle -profile -in1 $seq1 -in2 $seq2 -quiet -out $path_train/$out.both.aln");
		return "$path_train/$out.both.aln";
		}
	
		if ($method =~ /clustalw/)
		{
		my $factory = Bio::Tools::Run::Alignment::Clustalw -> new  (@params);
		my $aln=$factory->profile_align($seq1,$seq2);
		return $aln;
		}
	}
	else
	{	
		if ($method =~ /muscle/)
		{
		my $factory = new Bio::Tools::Run::Alignment::Muscle (@params);
		my $aln=$factory->align($seq1);
		return $aln;
		}
	
		if ($method =~ /clustalw/)
		{
		my $factory = Bio::Tools::Run::Alignment::Clustalw -> new  (@params);
		my $aln=$factory->align($seq1);
		return $aln;
		}
	
	}
	

}

# =head2 nsCluster
 
#  Title   : nsCluster
#  Usage   : $HMMmodEObj->nsCluster();
#  Function: Generate distinct clusters with non-single (ns) members
#  Returns : 
#  Args    : None
 
# =cut
 
sub nsCluster 
{
my ($self) = @_;
my $file=$self->{'file'};
my $ref_hash=$self->{'cluster_ids'};
my %clus_hash=%$ref_hash;
my @ns=();

foreach (sort keys %clus_hash)
{
my $idref=$self->{'cluster_seq'}{$_};
my @tempid=@$idref;
my $outfile=$self->{'sample_outfile'}{$_};	
	if (scalar @tempid >= 3)
	{
	$self->writeSeq($_,"$outfile".".seq");
	push (@ns,$_);
	}
	else
	{		
	$self->writeSeq($_,"$file".".rest","c");
	}
}
print "@ns >>>\n";
$self->{'ns_clusters'}=\@ns;
return \@ns;
}

# =head2 generateSamples
 
#  Title   : generateSamples
#  Usage   : $HMMmodEObj->generateSamples
#  Function: Generate n test and train samples of a cluster 'cluster number'; 
#            This is exhaustive sampling meaning each sequence appears at least once as part of a test set;
#            Default of 10 samples is good enough for most cases
#		For further details see ref     
#		Srivastava PK*, Desai DK*, Nandi S and Lynn AM BMC Bioinformatics 2007, 8:104 (27 March 2007)

#  Returns : Writes out the sequence files for the n test and train samples to 
#            '/path/to/the/test|train/samples/directory';
#  	     When called with a cluster number returns a hash of sample-number keys with ref_to_array of test or train sequence objects
#  Args    : 	=>Cluster number, 
#		=>'/path/to/the/test|train/samples/directory' 
#		=> n (optional) ; default is 10;
 
# =cut


 
sub generateSamples
{
my ($self,$cluster_no,$ref_fp) = @_;

my $path_test=$self->{'path_test'};
my $path_train=$self->{'path_train'};
my %hash_index=();
my $samples_tp;
my $samples_fp;
my %sample_test=();
my %sample_train=();
my %sample_fp=();
#Added by Andrew - 6/7/2014
my %sample_fptrain=();
my $outfile=$self->{'sample_outfile'}{$cluster_no};	


my $seqref=$self->{'cluster_seq'}{$cluster_no};
my @seq=@$seqref;
my $tot_seq=scalar @seq;

my @fp_seq=@$ref_fp;
my $tot_fp=scalar @fp_seq;
#print "$tot_fp ....fps in samples\n";
my $size_test_tp;
my $size_test_fp;
my $n;
my $flag;
my $exhaustive_samples_tp;
my $exhaustive_samples_fp;
my $left_over_tp;
my $left_over_fp;

	if ($tot_seq >= 10 && $tot_fp >= 10)
	{
	$n=10;
	$size_test_tp=ceil($tot_seq/$n);
	$size_test_fp=ceil($tot_fp/$n);
	$exhaustive_samples_tp=int($tot_seq/$size_test_tp);
	$exhaustive_samples_fp=int($tot_fp/$size_test_fp);
	$flag="a";
	}
	
	elsif ($tot_seq < 10 && $tot_fp < 10)
	{
		if ($tot_fp<$tot_seq)
		{
		$exhaustive_samples_fp=$n=$tot_fp;
		if ($n){$size_test_tp=ceil ($tot_seq/$n);}else {$n=0;}
		if ($size_test_tp){$exhaustive_samples_tp=int($tot_seq/$size_test_tp);}else {$size_test_tp=0}	
		$size_test_fp=1;
		}
		else
		{
		$exhaustive_samples_tp=$n=$tot_seq;
		if ($n){$size_test_fp=ceil ($tot_fp/$n);}else {$n=0;}
		if ($size_test_fp){$exhaustive_samples_fp=int($tot_fp/$size_test_fp);}else {$size_test_fp=0;}
		$size_test_tp=1;	
		}
	$flag="b";
	}
	
	elsif ($tot_seq>=10 && $tot_fp <10)
	{
	$n=$exhaustive_samples_fp=$tot_fp;
	if ($n) {$size_test_tp=ceil($tot_seq/$n);}else {$n=0;}
	$size_test_fp=1;
	if ($size_test_tp){$exhaustive_samples_tp=int ($tot_seq/$size_test_tp);}else {$size_test_tp=0;}
	$flag="c";
	}
	
	elsif ($tot_seq < 10 && $tot_fp>=10)
	{
	$n=$exhaustive_samples_tp=$tot_seq;
	$size_test_tp=1;
	if ($n) {$size_test_fp=ceil($tot_fp/$n);}else {$n=0;}
	if ($size_test_fp){$exhaustive_samples_fp=int($tot_fp/$size_test_fp);}else {$size_test_fp=0;}
	$flag="d";
	}
	
if ($size_test_tp){$left_over_tp=$tot_seq % $size_test_tp;}else {$left_over_tp=0;}
if ($size_test_fp){$left_over_fp=$tot_fp % $size_test_fp; }else {$left_over_fp=0;}
$samples_tp+=$exhaustive_samples_tp;
$samples_fp+=$exhaustive_samples_fp;
print "$n sample sets\t$tot_seq total sequences\t$size_test_tp => tp sample size\t$size_test_fp => fp sample size\t$exhaustive_samples_tp tp ex-samples possible\t$exhaustive_samples_fp fp ex-sampels possible\t$left_over_tp tp sequence(s) will be left out\t$left_over_fp fp sequence(s) will be left out \n";

my $counter=0;
	for (my $set=1;$set<=$exhaustive_samples_tp;$set++)
	{
	my %hash_index=();
	my @test=();
	my @train=();
		for (my $test_seq=$counter;$test_seq<$counter+$size_test_tp;$test_seq++)
		{
		#print "test begins with $counter in set $set\n";
		push (@test,$seq[$test_seq]);
		$hash_index{$seq[$test_seq]}=1;
		}
	$sample_test{$set}=\@test;	
	my $out_test = Bio::SeqIO->new('-file' => ">$path_test/$outfile.$set.test.seq", '-format' => 'Fasta');
	$out_test->write_seq(@test);
		for (my $train=0;$train<$tot_seq;$train++)
		{
			if ($hash_index{$seq[$train]} != 1)
			{
			#print $seq[$train],"\n";
			push (@train,$seq[$train]);
			}
		}
 		
	$sample_train{$set}=\@train;	
	my $out_train = Bio::SeqIO->new('-file' => ">$path_train/$outfile.$set.train.seq", '-format' => 'Fasta');
	$out_train->write_seq(@train);
	$counter += $size_test_tp;	
	}

my $counter=0;
	for (my $set=1;$set<=$exhaustive_samples_fp;$set++)
	{
	my %hash_index=();
	my @test=();
	my @train=();
		for (my $test_seq=$counter;$test_seq<$counter+$size_test_fp;$test_seq++)
		{
		#print "test begins with $counter in set $set\n";
		push (@test,$fp_seq[$test_seq]);
# Modified by Andrew 6/7/2014 - to change earlier $seq to $fp_seq. This was probably not used in the earlier script...?
		$hash_index{$fp_seq[$test_seq]}=1;
		}
 	
	$sample_fp{$set}=\@test;	
	my $out_test = Bio::SeqIO->new('-file' => ">$path_test/$outfile.$set.fp.seq", '-format' => 'Fasta');
	$out_test->write_seq(@test);
# Added by Andrew 6/7/2014 ---start
# added output file $outfile.$set.fptrain.seq and variable $sample_fptrain{$set}
		for (my $train=0;$train<$tot_fp;$train++)
		{
			if ($hash_index{$fp_seq[$train]} != 1)
			{
			#print $seq[$train],"\n";
			push (@train,$fp_seq[$train]);
			}
		}
 		
	$sample_fptrain{$set}=\@train;	
	my $out_train = Bio::SeqIO->new('-file' => ">$path_train/$outfile.$set.fptrain.seq", '-format' => 'Fasta');
	$out_train->write_seq(@train);
#Added by Andrew 6/7/2014 - end
	$counter += $size_test_fp;	
	}

	if ($left_over_tp)
	{
	my ($ref_test,$ref_train,$set_fill)=$self->_filler_samples($cluster_no,$left_over_tp,$size_test_tp,\@seq,$samples_tp,"test",1);
	$samples_tp=$set_fill;
	my @test=@$ref_test;
	my @train=@$ref_train;
	$sample_train{$samples_tp}=\@train;
	$sample_test{$samples_tp}=\@test;
	}
 	print "current tp set number is $samples_tp\n";
 	  
	if ($left_over_fp)
	{
	my ($ref_test,$set_fill)=$self->_filler_samples($cluster_no,$left_over_fp,$size_test_fp,\@fp_seq,$samples_fp,"fp");
	$samples_fp=$set_fill;
	my @test=@$ref_test;
	$sample_fp{$samples_fp}=\@test;
	}
 	print "current fp set number is $samples_fp\n";
 	  
 	
	if ($samples_tp<$n)
	{
	#my $diff=$n-$samples_tp;
#	print "$diff more sample(s) to go....\n";
	
		while ($samples_tp != $n)
		{
		$samples_tp++;
		my ($ref_test,$ref_train)=$self->splitSample($cluster_no,$size_test_tp);
		my @test=@$ref_test;
		my @train=@$ref_train;
		
#		print "split @test\n";
		
		$sample_test{$samples_tp}=\@test;
		$sample_train{$samples_tp}=\@train;
		
		my $out_test = Bio::SeqIO->new('-file' => ">$path_test/$outfile.$samples_tp.test.seq", '-format' => 'Fasta');
		$out_test->write_seq(@test);
		
		my $out_train = Bio::SeqIO->new('-file' => ">$path_train/$outfile.$samples_tp.train.seq", '-format' => 'Fasta');
		$out_train->write_seq(@train);
		
		}
	
	}
	
	if ($samples_fp<$n)
	{
	my $diff=$n-$samples_fp;
	
		while ($samples_fp != $n)
		{
		$samples_fp++;
		my @rand_array=_gen_rand_no($size_test_fp,($#fp_seq+1));
		my @uniq=_remove_dupes(($#fp_seq+1),@rand_array);
	
		my @rand_seq= map {$fp_seq[$_]} @uniq;
		$sample_fp{$samples_fp}=\@rand_seq;
		my $out_fp = Bio::SeqIO->new('-file' => ">$path_test/$outfile.$samples_fp.fp.seq", '-format' => 'Fasta');
		$out_fp->write_seq(@rand_seq);
		}
	
	}

$self->{'test_samples'}{$cluster_no}=\%sample_test;
$self->{'train_samples'}{$cluster_no}=\%sample_train;
$self->{'fp_samples'}{$cluster_no}=\%sample_fp;
#Added by Andrew 6/7/2014
$self->{'fptrain_samples'}{$cluster_no}=\%sample_fptrain;	
}


sub _filler_samples
{
my ($self,$cluster_no,$left_over_seq,$sample_size,$seqref,$current_set,$string,$flag)=@_;

#print "In filler...$left_over_seq, $string, $flag, $current_set, @$seqref\n";
my $path_train=$self->{'path_train'};
my $path_test=$self->{'path_test'};
my $outfile=$self->{'sample_outfile'}{$cluster_no};

my %hash_index=();
my @test=();
my @train=();
my $set=$current_set+1;
my @seq=@$seqref;
my $tot_seq=scalar @seq;
my @left_out=splice(@seq,-$left_over_seq);
push (@test,@left_out);
	
	
my $remaining=$sample_size-$left_over_seq;
my @rand_array=_gen_rand_no($remaining,($#seq+1));
my @uniq=_remove_dupes(($#seq+1),@rand_array);
my @rand_seq= map {$seq[$_]} @uniq;
push (@test,@rand_seq);
my $out_test = Bio::SeqIO->new('-file' => ">$path_test/$outfile.$	set.$string.seq", '-format' => 'Fasta');
$out_test->write_seq(@test);
	
	if ($flag)
	{
	@hash_index{@test}=();
		for (my $train=0;$train<$tot_seq;$train++)
		{
		next if (exists $hash_index{$seq[$train]});
			if ($seq[$train])
			{
			push (@train,$seq[$train]);
			}
		}
	my $out_train = Bio::SeqIO->new('-file' => ">$path_train/$outfile.$set.train.seq", '-format' => 'Fasta');
	$out_train->write_seq(@train);
	return (\@test,\@train,$set);
	}
	else
	{
	return (\@test,$set);
	}
}


# =head2 loadSamples
 
#  Title   : loadSamples
#  Usage   : $HMMmodEObj->loadSamples ($cluster_no,$n,$path_test,$path_train)
#  Function: Load test and train samples into the $HMMmodEObj object
#  Returns : 
#  Args    : 	$cluster_no => Cluster number, 
#		$path_test  =>'/path/to/the/test_samples/directory' 
#	        $path_train =>'/path/to/the/train_samples/directory' 
#		$n => n (optional) ; default is 10;
 
# =cut


sub loadSamples
{
my ($self,$cluster_no,$n,$path_test,$path_train) = @_;
my %sample_test=();
my %sample_train=();

	if (!$self->{'path_test'} && !$self->{'path_train'})
	{
	$self->{'path_test'}=$path_test;
	$self->{'path_train'}=$path_train;
	}
my $tot_samples=$n?$n:10;
my $outfile=$self->{'sample_outfile'}{$cluster_no};
print "out ", $outfile,"\n";	
	for (my $samples=1;$samples<=$tot_samples;$samples++)
	{
	my @test=();
	my @train=();
	
	my $in_test = Bio::SeqIO->new('-file' => "$path_test/$outfile.$samples.test.seq", '-format' => 'Fasta');
		
		while (my $testseq = $in_test->next_seq)
		{
		push (@test,$testseq);
		}
	
	my $in_train = Bio::SeqIO->new('-file' => "$path_train/$outfile.$samples.train.seq", '-format' => 'Fasta');
		
		while (my $seq = $in_train->next_seq)
		{
		push (@train,$seq);
		}
	
	$sample_test{$samples}=\@test;
	$sample_train{$samples}=\@train;	
	}
	
$self->{'test_samples'}{$cluster_no}=\%sample_test;
$self->{'train_samples'}{$cluster_no}=\%sample_train;
}

sub loadPath
{
my ($self,$cluster_no,$path_test,$path_train) = @_;
	
	if (!$path_test && !$path_train)
	{
	print "Path to store test and train samples not provided; will attempt to create /tmp/test_samples and /tmp/train_samples\n";
	_check_directory("/tmp/test_samples");
	_check_directory("/tmp/train_samples");
	$self->{'path_test'}="/tmp/test_samples";
	$self->{'path_train'}="/tmp/train_samples";
	}
	else
	{
	_check_directory("$path_test");
	_check_directory("$path_train");
	$self->{'path_test'}=$path_test;
	$self->{'path_train'}=$path_train;
	}
}

sub path
{
my ($self,$cluster_no,$path_train) = @_;

        if (!$path_train)
        {
        print "Path to store tp profiles; will attempt to create /tmp/train_samples\n";
        _check_directory("/tmp/train_samples");
        $self->{'path_train'}="/tmp/train_samples";
        }
        else
        {
        _check_directory("$path_train");
        $self->{'path_train'}=$path_train;
        }
}

# =head2 _get_method_ref
 
#  Title   : _get_method_ref($obj,$method)
#  Usage   : Only to be called internally
#  Function: Get a reference to a method of an object
#  Returns : a reference to a sub-routine called by an object
#  Args    : 	$obj -> any object of a class
#		$method -> a method called by $obj
 
# =cut

sub _get_method_ref
{
my ($self1,$methodname)=@_;
my $methref = sub {return $self1->$methodname(@_);};
return $methref;
}


# =head2 _gen_rand_no
 
#  Title   : _gen_rand_no($size,$limit)
#  Usage   : Only to be called internally
#  Function: Generating an array of random no.s
#  Returns : an array of random numbers
#  Args    : 	$size -> length of the array
#		$limit -> upper limit of random integers
 
# =cut


sub _gen_rand_no
{
my ($array_len,$limit)=@_;

my $len=();
my @array=();
	while ($len < $array_len)
	{
	my $rand_no=int(rand($limit));
	push (@array,$rand_no);
	$len = scalar(@array);
	}
return @array;
}


# =head2 _remove_dupes
 
#  Title   : _remove_dupes($limit,@rand_array)
#  Usage   : Only to be called internally
#  Function: Removing duplicate entries from an array of random no.s and adding new ones to make up the size
#  Returns : an array of unique random numbers
#  Args    : 	 \@rand_array -> array of rnadom numbers
#		$limit -> upper limit of random integers
 
# =cut



sub _remove_dupes
{

my $u_lim=shift;
my @old_array=@_;
my %seen=();
my $old_len=scalar (@old_array);
my @uniq = grep(!$seen{$_}++,@old_array);
my $new_len=scalar (@uniq);
my $diff=$old_len-$new_len;
        if ($diff > 0)
        {
        my @additional=_gen_rand_no($diff,$u_lim);
        splice (@uniq,($#uniq+1),0,@additional);
        _remove_dupes($u_lim,@uniq);
        }
        else
        {
        return @uniq;
        }

}


# =head2 sub _check_directory
# 
#  Title   : _check_directory($path)
#  Usage   : Only to be called internally
#  Function: Check if /path/to/sample/test|train/directory exists; If the directory names are not provided then tries to create 2 directories test_samples and train_samples in /tmp to store the test and train sample sequence files
#  Returns : Error msg if /path/to/sample/test|train/directory cannot be found in the current directory OR it fails to #	     create the directories in /tmp
#  Args    : /path/to/sample/test|train/directory
# 
# =cut
# 
# 

sub _check_directory
{
my $directory=shift;
# checking if the test directory already exists
	if (-d $directory)
	{ # the -d bit means check if there is a directory called ...
		print "$directory already exists, samples will be saved in $directory\n";
	}
	else
	{
	print "$directory does not exist.....creating $directory.....\n\n";
	mkdir ($directory, 0777) || die "sorry system is unable to create output directory $directory";
	}

}
 

# =head2 splitSample
 
#  Title   : splitSample($sample_size,$samples)
#  Usage   : Only to be called internally
#  Function: Generating split samples 
#  Returns : references to arrays of test and train sequences
#  Args    : 	$sample_size -> number of test sequences
#		$cluster_no -> integer indicating the cluster number
 
# =cut

sub splitSample
{

my ($self,$cluster_no,$size)=@_;
my %hash_index=();
my (@test,@train)=();
my $seqref=$self->{'cluster_seq'}{$cluster_no};
my @seq = @$seqref;


my $tot_seq=scalar @seq;

my @rand_array=_gen_rand_no($size,$tot_seq);
my @uniq=_remove_dupes($tot_seq,@rand_array);

my @rand_seq=map {$seq[$_]} @uniq;


@hash_index{@rand_seq}=();
	for (my $train=0;$train<$tot_seq;$train++)
	{
	next if (exists $hash_index{$seq[$train]});
	push (@train,$seq[$train]);
	}

return (\@rand_seq,\@train);

}
# =head2 crossValidate
 
#  Title   : crossValidate
#  Usage   : $HMMmodEObj->crossValidate(cluster id)
#  Function: Train models of HMMER 3 from the train samples and test them on test samples, 
#	     parse Hmm hits and generate a cumulative histogram of Sensitivity, specificity and MCC for each set
#	     Calculate the optimised threshold for the given cluster (Mode of the avg MCC distribution) 	    
#  Returns : scalar- optimised cut-off score for a given cluster
#  Args    : Cluster id
 
# =cut
 
 
sub crossValidate
{
my ($self,$cluster_no) = @_;

my %cumu_mcc;
my %seen;
my %cumu_mode_mcc=();
my %seen_mode=();
my $seqfile=$self->{'file'};
my $ref_test_hash=$self->{'test_samples'}{$cluster_no};
my %test_hash=%$ref_test_hash;

my $ref_train_hash=$self->{'train_samples'}{$cluster_no};
my %train_hash=%$ref_train_hash;

my $ref_fptrain_hash=$self->{'fptrain_samples'}{$cluster_no};
my %fptrain_hash=%$ref_fptrain_hash;


my $ref_fp_hash=$self->{'fp_samples'}{$cluster_no};
my %fp_hash=%$ref_fp_hash;

my $path_test=$self->{'path_test'};
my $path_train=$self->{'path_train'};
my $outfile=$self->{'sample_outfile'}{$cluster_no};
my @samples_test=sort {$a<=>$b} keys %test_hash;
my @samples_fptrain=sort {$a<=>$b} keys %fptrain_hash;
my @samples=();
print scalar @samples_test;
print "\n";
print scalar @samples_fptrain;
print "\n";
if (scalar @samples_test >   scalar @samples_fptrain)
{
	@samples=@samples_fptrain;
}
else
{
	 @samples=@samples_test;
}
my $mode_counter=0;
	foreach my $s_no(@samples)
	{
		my $ref_train_array=$train_hash{$s_no};
		my $ref_test_array=$test_hash{$s_no};
		my $ref_fp_seq=$fp_hash{$s_no};
# Get the test sequences into a array
		my @tmp_test_ids=@$ref_test_array;
        	my @test_ids=map{$_->display_id()}@tmp_test_ids;
		foreach(@test_ids)
		{
#		print "test $_\n";
		}
# get the fp sequences into a array
		my @tmp_fp_ids=@$ref_fp_seq;
        	my @fps_ids=map{$_->display_id()}@tmp_fp_ids;
        	foreach(@fps_ids)
        	{
 #               print "fp $_\n";
        	}
#initialise %tp_test_hash to contain @tp_ids as keys and values set to 0
                my %tp_test_hash;
                foreach(@test_ids)
                {
                        $tp_test_hash{$_}=0;
                }
#              while ( my ($key, $value) = each(%tp_test_hash) ) {
#              print "inside TP CV3 $key => $value\n";
#              }
#initialise %fp_test_hash to contain @fp_ids as keys and values set to 0
                my %fp_test_hash;
                foreach(@fps_ids)
                {
                        $fp_test_hash{$_}=0;
                }
		my %fp_train_hash;
                foreach(@fps_ids)
                {
                        $fp_train_hash{$_}=0;
                }
#               while ( my ($key, $value) = each(%fp_test_hash) ) {
#               print "inside FP CV3 $key => $value\n";
#               }
		my @params = (OUT => "$path_train/$outfile.$s_no.tp.aln", MAXITERS => 2,"QUIET" => "QUIET" ,MAXMB => 3800);
		my $ref_params=\@params;
		my $train_aln=$self->mAlign("/HOME_ann/sw/usr/bin/muscle",$ref_params,$cluster_no,$ref_train_array);
        	system ("/HOME_ann/sw/usr/hmmer3/hmmer-3.0-linux-intel-x86_64/binaries/hmmbuild --informat afa $path_train/$outfile.$s_no.hmm $path_train/$outfile.$s_no.tp.aln");
		my $ref_fp_test_hash=$self->hmmsearchParse_CV("$path_test/$outfile.$s_no.fp.seq","$path_train/$outfile.$s_no.hmm","neg.hmms");
		my %new_fp_train_hash=(%fp_train_hash,%$ref_fp_test_hash);
		my @fp_train_score=values %new_fp_train_hash;
		foreach(@fp_train_score)
		{
			print "fp train score is $_\n";
		}
		my $ref_tp_test_hash=$self->hmmsearchParse_CV("$path_test/$outfile.$s_no.test.seq","$path_train/$outfile.$s_no.hmm", "pos.hmms",-9999);
		my %new_tp_test_hash=(%tp_test_hash,%$ref_tp_test_hash);
		my @tp_score=values %new_tp_test_hash;
		my %new_fp_test_hash=(%fp_test_hash,%$ref_fp_test_hash);
		my @fp_score=values %new_fp_test_hash;
		my $ref_mcc=$self->writeHisto(\@tp_score,\@fp_score,"$path_train/$outfile.$s_no.hmmer.histo");
		%cumu_mcc=%$ref_mcc;
		foreach (keys %cumu_mcc)
		{
		$seen{$_}=$seen{$_}+$cumu_mcc{$_};
		}
		if (defined $fp_train_score[0])
		{
			my @fp=();
	#		@fp=@$ref_fp_seq;
			my $ref_fptrain_seq=$fptrain_hash{$s_no};
			@fp=@$ref_fptrain_seq;
			my @params = (OUT => "$path_train/$outfile.$s_no.fp.aln", MAXITERS => 2,"QUIET" => "QUIET", MAXMB => 3800);	
			my $fp_aln=$self->mAlign("/HOME_ann/sw/usr/bin/muscle",\@params,$cluster_no,\@fp);
			my $aligned=$self->mAlign("/HOME_ann/sw/usr/bin/muscle",\@params,$cluster_no,"$path_train/$outfile.$s_no.tp.aln","$path_train/$outfile.$s_no.fp.aln","p");
			my @tp=@$ref_train_array;
 			my @tp_ids=map{$_->display_id()}@tp;
			my @fp_ids=map{$_->display_id()}@fp;
			my $db = new Bio::DB::Fasta($aligned);  # one file or many files
			my $tp_out=Bio::SeqIO->new('-file' => ">$path_train/$outfile.$s_no.tp.prof.aln",'-format' => 'Fasta');
			foreach my $id (@tp_ids)
			{
				my $seqobj = $db->get_Seq_by_id($id); # get a PrimarySeq obj
				if($seqobj)
				{
				$tp_out->write_seq($seqobj);
				}
			}
			my $fp_out=Bio::SeqIO->new('-file' => ">$path_train/$outfile.$s_no.fp.prof.aln",'-format' => 'Fasta');
			foreach my $id (@fp_ids)
			{
			my $seqobj = $db->get_Seq_by_id($id); # get a PrimarySeq obj
				if($seqobj)
				{
				$fp_out->write_seq($seqobj);
				}
			}
 		system ("/HOME_ann/sw/usr/hmmer3/hmmer-3.0-linux-intel-x86_64/binaries/hmmbuild --informat afa $path_train/$outfile.$s_no.tp.prof.hmm $path_train/$outfile.$s_no.tp.prof.aln");
 	 	system ("/HOME_ann/sw/usr/hmmer3/hmmer-3.0-linux-intel-x86_64/binaries/hmmbuild --informat afa $path_train/$outfile.$s_no.fp.prof.hmm $path_train/$outfile.$s_no.fp.prof.aln");
		my $hmmmodefile=$self->modE($cluster_no,$s_no,"$path_train/$outfile.$s_no.tp.prof.hmm","$path_train/$outfile.$s_no.fp.prof.hmm");
		my $ref_tp_test_hash=$self->hmmsearchParse_CV("$path_test/$outfile.$s_no.test.seq","$hmmmodefile","mode.pos.hmms",-9999);
# %tp_test_hash should already contain IDs as keys and values initialised to 0
#		my %tp_test_hash=%$ref_tp_test_hash;
		my %new_tp_test_hash=(%tp_test_hash,%$ref_tp_test_hash);
		#while ( my ($key, $value) = each(%new_tp_test_hash) ) {
                #print "inside TP test hash $key => $value\n";
                #}
		my @tp_score=values %new_tp_test_hash;
		my $ref_fp_test_hash=$self->hmmsearchParse_CV("$path_test/$outfile.$s_no.fp.seq","$hmmmodefile","mode.neg.hmms",-9999);
# similarly %fp_test_hash should contain FP IDs as keys and values initialised to 0
#		my %fp_test_hash=%$ref_fp_test_hash;
		my %new_fp_test_hash=(%fp_test_hash,%$ref_fp_test_hash);
#print fp_test_hash
		#while ( my ($key, $value) = each(%new_fp_test_hash) ) {
                #print "inside FP test hash $key => $value\n";
                #}
		my @fp_score=values %new_fp_test_hash;
		my $ref_mode_mcc=$self->writeHisto(\@tp_score,\@fp_score,"$path_train/$outfile.$s_no.mode.histo");
		%cumu_mode_mcc=%$ref_mode_mcc;
		
			foreach (keys %cumu_mode_mcc)
			{
			$seen_mode{$_}=$seen_mode{$_}+$cumu_mode_mcc{$_};
			}
		$mode_counter++;
		}	
	}
my $samples=scalar @samples;	
my $hmmer_cut=$self->_thresholdScore(\%seen,$cluster_no,"hmmer",$samples);
	if ($mode_counter)
	{
	my $mode_cut=$self->_thresholdScore(\%seen_mode,$cluster_no,"mode",$mode_counter);
	return ($hmmer_cut,$mode_cut);
	}
	else
	{
	return $hmmer_cut;
	}
}
 
 
# =head2 modE
 
#  Title   : modE
#  Usage   : $HMMmodEObj->modE("path/to/true_positive/profile_aligned/hmm","path/to/false_positive/profile_aligned/hmm")
#  Function: Given the profile_aligned TP and FP hmm files, identifies the columns which are similarly conserved in both using a relative entropy filter
#	     Modifies the emission probabilities of the common columns of TP profile with those from the same columns in the FP profile 
#	     Writes out the modified profile hmm file of HMMER 3 format
#  Returns : $path/to/the/modified/hmmfile
#  Args    : "path/to/true_positive/profile_aligned/hmm"
#            "path/to/false_positive/profile_aligned/hmm"

# =cut 
 
 
sub modE
{
	my ($self,$cluster_no,$sample_id,$file_tp,$file_fp)=@_;
	my $outfile=$self->{'sample_outfile'}{$cluster_no};
	$outfile=$outfile.".$sample_id".".ModE";
	my $path_train=$self->{'path_train'};
	my $fp_hmm = new Hmm3(-file=>"$file_fp");
	my $tp_hmm = new Hmm3(-file=>"$file_tp");
	my @maps_tp=();
	my @maps_fp=();
	my $sum=0;
	while ( my ($k0, $v0) = each %{$tp_hmm->matchEmissionMap})
	{
		while ( my ($k1, $v1) = each %{$fp_hmm->matchEmissionMap})
		{
			if($v0 == $v1)
			{
				push(@maps_tp,$k0);
				push(@maps_fp,$k1);
			}
		}
	}
	my $arr_map_tp_size=@maps_tp;
	for(my $i=0;$i<=$#maps_tp;$i++)
	{   
	my $REPositiveByNegative=0;
		for(my $j=0;$j<=scalar @{$fp_hmm->{'BASE'}} - 1;$j++)
		{
			$REPositiveByNegative+= ${${$tp_hmm->matchEmissionProbablity}{$maps_tp[$i]}}{${$tp_hmm->{'BASE'}}[$j]} * log (${${$tp_hmm->matchEmissionProbablity}{$maps_tp[$i]}}{${$tp_hmm->{'BASE'}}[$j]}/${${$fp_hmm->matchEmissionProbablity}{$maps_fp[$i]}}{${$fp_hmm->{'BASE'}}[$j]});
		}
	my $REPositiveByNull=0;
		for(my $k=0;$k<=scalar @{$fp_hmm->{'BASE'}} - 1;$k++)
		{
			$REPositiveByNull += ${${$tp_hmm->matchEmissionProbablity}{$maps_tp[$i]}}{${$tp_hmm->{'BASE'}}[$k]} * log ( ${${$tp_hmm->matchEmissionProbablity}{$maps_tp[$i]}}{${$fp_hmm->{'BASE'}}[$k]} / ${$fp_hmm->nuleProbablity}{${$fp_hmm->{'BASE'}}[$k]} );
		}
		if($REPositiveByNegative > $REPositiveByNull)
		{
			my %newNull;
			my %newP;
			my $denomenator=0;
			for(my $k=0;$k<=scalar @{$fp_hmm->{'BASE'}} - 1;$k++)
			{
$newNull{ ${$fp_hmm->{'BASE'}}[$k] } = (${$fp_hmm->nuleProbablity}{${$fp_hmm->{'BASE'}}[$k]} > ${${$fp_hmm->matchEmissionProbablity}{$maps_fp[$i]}}{${$fp_hmm->{'BASE'}}[$k]} ) ? ${$fp_hmm->nuleProbablity}{${$fp_hmm->{'BASE'}}[$k]} : ${${$fp_hmm->matchEmissionProbablity}{$maps_fp[$i]}}{${$fp_hmm->{'BASE'}}[$k]};
				if (${$fp_hmm->nuleProbablity}{${$fp_hmm->{'BASE'}}[$k]} > ${${$fp_hmm->matchEmissionProbablity}{$maps_fp[$i]}}{${$fp_hmm->{'BASE'}}[$k]} )
				{
					$newP{ ${$tp_hmm->{'BASE'}}[$k] }=${${$tp_hmm->matchEmissionProbablity}{$maps_tp[$i]}}{${$tp_hmm->{'BASE'}}[$k]};
				}
				else
				{
			$newP{ ${$tp_hmm->{'BASE'}}[$k] }=((${${$tp_hmm->matchEmissionProbablity}{$maps_tp[$i]}}{${$tp_hmm->{'BASE'}}[$k]}*${$fp_hmm->nuleProbablity}{${$fp_hmm->{'BASE'}}[$k]})/$newNull{ ${$fp_hmm->{'BASE'}}[$k] });
				}
			$denomenator +=  $newP{ ${$tp_hmm->{'BASE'}}[$k] };
			}
			while (my ($k, $v) = each %newP)
                        {
                        	$newP{$k}=$v/$denomenator;
                       		${${$tp_hmm->matchEmission}{$maps_tp[$i]}}{$k} = -log($newP{$k});
                	}
		}	
	} 
$fp_hmm->destroy;
$tp_hmm->writeProfile("$path_train/$outfile");
$tp_hmm->destroy;
return ("$path_train/$outfile");
} 
 
 
 
 
# =head2 _thresholdScore
 
#  Title   : _thresholdScore
#  Usage   : Only to be called internally-$self->_thresholdScore(ref_to_hash of cumulative mcc for all "n" sets, $cluster_no, $string, $samples)
#  Function: Identifies the cut-off score from the average of the cumulative mcc distribution for all "n" sets 
#	     Cut-off is the median of scores corresponding to the max value of the mcc distribution averaged over "n" #            sets
#  Returns : scalar => $cutoff_score
#  Args    : ref_to_hash of cumulative mcc for all "n" sets
#            $cluster_no => cluster number
#	     $string => "hmmer" OR "ModE" depending on the profile 
#	     $samples => total number of samples
# =cut
 
 
sub _thresholdScore
{
my ($self,$ref_hash,$cluster_no,$string,$samples)=@_;
my %avg;
my %mcc=%$ref_hash;
my $path_train=$self->{'path_train'};
my $outfile=$self->{'sample_outfile'}{$cluster_no};
$outfile=$outfile.".$string".".avg";
open (OUT ,">$path_train/$outfile");

	foreach (sort {$a<=>$b} keys %mcc)
	{
	my $avg=$mcc{$_}/$samples;
	print OUT "$_\t$avg\n";
	$avg{$avg}=$avg{$avg}.",".$_;
	}

my @avg_values=keys %avg;
my $max = max(@avg_values);
my @scores=split (",",$avg{$max});
my $median=median(@scores);
#print "max $max cutoff $median \n";
return $median;
} 
 

# =head2 hmmsearchParse3
 
#  Title   : hmmsearchParse3
#  Usage   : $HMMmodEObj->hmmsearchParse(ref_to_array of seqobjs, hmm model, $cutoff)
#  Function: runs hmmsearch on the array of seqobjs 
#	     parse the results to get hit ids and scores >= $cutoff
#  Returns : ref_to_hash of ids and scores id=>score
#  Args    : ref_to_array of seqobjs => the sequences to be scored with the model
#            hmm model => path/to/the/hmm/file  
#            $cutoff => the cut-off score; default is 0

# =cut


sub hmmsearchParse
{
my ($self,$seq,$hmm,$outfile,$cutoff)=@_;
my %score_hash=();
my $cutoffscore=$cutoff?$cutoff:0;
#my $cutoffscore=-99;
# --max    : Turn all heuristic filters off (less speed, more power)
#system("hmmsearch --cpu=1 --max $hmm $seq > $outfile");## Added the -max option for during the BMC research notes paper
system("/HOME_ann/sw/usr/hmmer3/hmmer-3.0-linux-intel-x86_64/binaries/hmmsearch $hmm $seq > $outfile");	
my $search = new Bio::SearchIO(-format => 'hmmer3', -file   => "$outfile");
		while (my $result = $search->next_result)
		{
			while(my $hit = $result->next_hit)
			{
				if ($hit->score() > $cutoffscore)
				{
				my $id=$hit->name();
				#print $id,"--\t",$hit->score(),"\n";
				$score_hash{$id}=$hit->score();
				}
			}
		}
return \%score_hash;
}



sub hmmsearchParse_CV
{
my ($self,$seq,$hmm,$outfile,$cutoff)=@_;
my %score_hash=();
#my $cutoffscore=$cutoff?$cutoff:0;
my $cutoffscore=-9999;
# --max    : Turn all heuristic filters off (less speed, more power)
system("/HOME_ann/sw/usr/hmmer3/hmmer-3.0-linux-intel-x86_64/binaries/hmmsearch --max $hmm $seq > $outfile");## Added the -max option for during the BMC research notes paper
#system("hmmsearch $hmm $seq > $outfile");       
my $search = new Bio::SearchIO(-format => 'hmmer3', -file   => "$outfile");
                while (my $result = $search->next_result)
                {
                        while(my $hit = $result->next_hit)
                        {
                          if ($hit->score() > $cutoffscore)
                                {
                               my $id=$hit->name();
                               print $id,"--\t",$hit->score(),"\n";
                               $score_hash{$id}=$hit->score();
                               }
                         }
                 }
                               return \%score_hash;
}                               

# =head2 _catSeq
 
#  Title   : _catSeq
#  Usage   : Only to be called internally 
#  Function: Cat sequences of samples of different groups with the same "sample_number" except for the group 
# 	     ;
#	     	For example: _catSeq(1,1,"test") cats  
#  Returns : An ref_to_array of seqobjs from the catted samples
#  Args    : cluster_id
#	     sample id,
#            $string => either "test" or "train" to cat test or train samples
 
# =cut

sub _catSeq
{
my ($self,$cluster_id,$sample_id,$string)=@_;
my %hash_index=();
my @cat=();
my @other_seqs=();
my $internal_flag;
my $seqin;

my $ref_own_seq=$sample_id && $string?$self->{$string."_samples"}{$cluster_id}{$sample_id}:$self->{'cluster_seq'}{$cluster_id};
my @own_seq=@$ref_own_seq;
	if ($string)
	{
	$internal_flag=1;
	}
my @own_ids=map{$_->display_id()}@own_seq;
@hash_index{@own_ids}=();
	if ($internal_flag)
	{
	my $path=$string=~/train/?$self->{'path_train'}:$self->{'path_test'};
	$seqin = Bio::SeqIO->new(-file   => "/bin/cat $path/*.$sample_id.$string.seq |", -format => 'Fasta');
	}
	else
	{
	$seqin = Bio::SeqIO->new(-file   => "/bin/cat ./*_*.seq |", -format => 'Fasta');
	}
	
	
	while (my $inseq = $seqin->next_seq)
	{
	push (@cat,$inseq);
	}
	for (my $s=0;$s<=$#cat;$s++)
	{
	my $seqobj=$cat[$s];
	my $id=$seqobj->display_id();
	next if (exists $hash_index{$id});
	push (@other_seqs,$cat[$s]);
	}
#print scalar @other_seqs,"  number of other seqs\n";	
return \@other_seqs;
}


# =head2 writeHisto
 
#  Title   : writeHisto
#  Usage   : $HMMmodEObj->writeHisto(ref_to_array of tp scores, ref_to_array of fp scores,$histo_out)
#  Function: Generates a cumulative distribution of sensitivity, specifictificity and MCC
#	     Stores the distributions in a text file; this file can be used to visualize the plots 
#  Returns : 
#  Args    : ref_to_array of tp scores => array of score ranges of the tp sequences
#            ref_to_array of tp scores => array of score ranges of the tp sequences
#            $histo_out => /path/to/histo/output/file

# =cut


sub writeHisto

{
my ($self,$ref_tp_score,$ref_fp_score,$fileout)=@_;

my $TP;
my $FP;
my $TN;
my $FN;
my $MIN;
my $MAX;
my %mcc=();

open (FH, ">$fileout");

my @data_tp=@$ref_tp_score;
my @data_fp=@$ref_fp_score;
my $size_tp=@data_tp;
my $size_fp=@data_fp;

#$size of $tphmms is TP+FN
#FN = $sizetphmms - $sumtp
#$sumtp of $tphmms is TP
#$sumfp of $fphmms is FP
# TN = $sizefphmms -$sumfp
# Calculate Sens = TP/(TP+FN) Spec = TP/TP+FP
#
# hash{range}={TP}=value
#            ={FP}=value
#
my %hash_tp = statshash(@data_tp);
my $start_range_tp = $hash_tp{'min'};
my $end_range_tp = $hash_tp{'max'};


my %hash_fp = statshash(@data_fp);
my $start_range_fp = $hash_fp{'min'};
my $end_range_fp = $hash_fp{'max'};

my $increment = 1;

	if($start_range_tp < $start_range_fp) 
	{
	$MIN = floor($start_range_tp)
	}
	else
	{
	$MIN = floor($start_range_fp)
	}
	
	if($end_range_tp > $end_range_fp) 
	{
	$MAX = floor($end_range_tp)
	}
	else
	{
	$MAX = floor($end_range_fp)
	}

my %bin;
my @range=();
	while($MIN <= $MAX)
	{
	push(@range,$MIN);
	$bin{$MIN}{'TP'}=0;
	$bin{$MIN}{'FP'}=0;
	$bin{$MIN}{'SENS'}=0;
	$bin{$MIN}{'SPEC'}=0;
	$bin{$MIN}{'MCC'}=0;
	$MIN += $increment;
	}

	foreach my $el(@data_tp)
	{
	chomp $el;
	my $element = floor($el);
		if($element =~ /^#/){ next}
	
	my $i=0;
		foreach my $range(@range)
		{
			if($element < ($range))
			{
			$bin{$range}{'TP'}++;
			}
		}
	}

	foreach my $el(@data_fp)
	{
	chomp $el;
	my $element = floor($el);
		if($element =~ /^#/){ next}
	
	my $i=0;
		foreach my $range(@range)
		{
			if($element < ($range))
			{
			$bin{$range}{'FP'}++;
			}
		}
	}

	foreach my $range (@range)
	{
	$TP=( $size_tp - $bin{$range}{'TP'} );
	$TN= $bin{$range}{'FP'};
	$FP=( $size_fp - $bin{$range}{'FP'} );
	$FN= $bin{$range}{'TP'};
		if ($TP+$FN)
		{
		$bin{$range}{'SENS'}=$TP/($TP+$FN);
		}
		if ($TP+$FP)
		{
		$bin{$range}{'SPEC'}=$TP/($TP+$FP);
		}

		if (($TP+$FN)*($TP+$FP)*($TN+$FP)*($TN+$FN))
		{
		$bin{$range}{'MCC'} = (($TP*$TN) - ($FP * $FN))/sqrt( ($TP+$FN)*($TP+$FP)*($TN+$FP)*($TN+$FN) );
		}
	$mcc{$range}=$bin{$range}{'MCC'};
	print FH $range,"\t",$bin{$range}{'SENS'},"\t",$bin{$range}{'SPEC'},"\t",$bin{$range}{'MCC'},"\t",$TP,"\t",$TN,"\t",$FP,"\t",$FN,"\n";
	}
return \%mcc;
}
1;
