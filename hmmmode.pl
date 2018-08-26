#!/usr/bin/perl
use strict;
use HMMmodE;				
use Hmm3;				
use Bio::Tools::Run::Alignment::Muscle;
use Bio::Tools::Run::Hmmer; 
use Bio::AlignIO;
use Bio::SeqIO;
use Data::Dumper;
use Bio::DB::Fasta;
use Statistics::Lite qw(:all);
my $seqfilename;
my $clusterObj;	
my %hash_info=();
my @fp_seq=();
my @own_seq;
my $file;
my $path_train;	
my $clus_no;
my $fpseqfile;
my @clus=();	
my @fp_final=();
my @fpscores=();
my @tpscores=();
#print "Please enter the file having list of pre-classified training sequences\n";
#my $list=<STDIN>;
my $list=shift;
open (LIST,$list);
## CHECK WHETHER all.fast exists,IF NOT PRINT A ERROR MESSAGE
my $filename = 'all.fasta';
if (-e $filename)
{
	 print "all.fasta file Exists!\n";
}
else
{
	print "ERROR::There is no all.fasta file, CREATE it and RE-RUN the code\n";
	exit;
}
## REMOVE DUPLICATES (IF ANY)
my $index_file='all.fasta.index';
if(!-e $index_file)
{
remove_dupes("all.fasta");
}
## CREATE A INDEX FASTA FILE 
my $db = new Bio::DB::Fasta("all.fasta");
## FILE TO WRITE OUT THE OPTIMIZED THRESHOLDS
open (OUT,">$list.list_cut.hmm"); 
open (CUTOFF,">$list.list_TC_NC");		
## LOOPS THROUGH EACH FILE
while ($seqfilename=<LIST>) 
{
	chomp $seqfilename;			
	my $clusfile=$seqfilename.".clus";		
	my $tabfile=$seqfilename.".blast.tab";
	$clusterObj = HMMmodE->new(-file => $seqfilename, -clusType => 'MCL', -clusFile => $clusfile, -tabFile => $tabfile);	
	my $ref=$clusterObj->clusInfo();					 
	%hash_info=%$ref;								
	print Dumper (\%hash_info);
	foreach $clus_no (sort keys %hash_info)			 					
	{
		print "Cluster number = $clus_no\n";
		my $outfile=$clusterObj->{'sample_outfile'}{$clus_no};
		print "name of out file is $outfile\n";
		print "running for cluster $clus_no ...\n";
		print $hash_info{$clus_no}, " sequences ----\n";
		my $ref_own=$clusterObj->{'cluster_seq'}{$clus_no};
		@own_seq=@$ref_own;				
		print "Sequence ids are @own_seq ======\n";
		if (scalar @own_seq >=3)			 
		{
			my %hash_index=();
			my @fp_ids=();
			my %fp_train_hash=();
			my @other_seqs=();
			my @own_ids=();
			$clusterObj->loadPath($clus_no,"./$list.test_samples","./$list.train_samples");	
			$path_train=$clusterObj->{'path_train'};
			@own_ids=map{$_->display_id()}@own_seq;				
			@hash_index{@own_ids}=();					
			my @params = (OUT => "$path_train/$outfile.aln",MAXITERS => 2,"QUIET" => "QUIET"); 
			my $train_aln=$clusterObj->mAlign("muscle",\@params,$clus_no,$ref_own);				
			
			#### GENERATE TP HMM PROFILE
			print "Generating TP HMM profile for $outfile\n";
			system ("hmmbuild --informat afa $path_train/$outfile.hmm $path_train/$outfile.aln");
			
			### SEARCH AND SCORE ALL THE TRAINING SEQUENCES AGAINST ALL SEQUENCES
			my $ref_fp_train_hash=$clusterObj->hmmsearchParse("all.fasta","$path_train/$outfile.hmm","$outfile.neg");
			%fp_train_hash=%$ref_fp_train_hash;
			
 			## GET THE FALSE POSITIVE SEQUENCES into an array @fp_ids
 			foreach (sort keys %fp_train_hash)		
 			{
 				$_=~s/\s+//g;
 				next if exists $hash_index{$_};
 				if ($_)
 				{
 					push (@fp_ids,$_);
					push (@fpscores,$fp_train_hash{$_}); 
 				}
 			}
			foreach (sort keys %fp_train_hash)
                        {
                                $_=~s/\s+//g;
                                if (exists $hash_index{$_})
                                {
                                        push (@tpscores,$fp_train_hash{$_});
                                }
                        }

			print "Number of False positive ids are","\t", scalar @fp_ids,"\n";			
			my $max_fp_NC= max(@fpscores);
			my $min_tp_TC=min(@tpscores);
			print CUTOFF "$seqfilename\tNC" ,"\t", $max_fp_NC,"\n"; 
			print CUTOFF "$seqfilename\tTC" ,"\t", $min_tp_TC,"\n";
			
			## IF FP>3 PROCEED FOR HMM-ModE
			if (defined $fp_ids[2])				 
			{
				@fp_seq=();
				foreach my $id (@fp_ids)
				{
					my $seqobj = $db->get_Seq_by_id($id); 	 		
					if ($seqobj)
					{
						push (@fp_seq,$seqobj);
					}
				}
			my $out_fp = Bio::SeqIO->new('-file' => ">$path_train/neg.seq", '-format' => 'Fasta');
			$out_fp->write_seq(@fp_seq);
			print "FP ids are","\n", "@fp_seq\n";
			$file="$path_train/neg.seq";
			
			## IF FPs MORE THAN 200,GET REPRESENTATIVVE 200 FPs BY RANDOM SAMPLING
			if (scalar @fp_ids > 200)		
			{
				#system "~/ncbi-blast-2.2.29+/bin/makeblastdb -in $file -dbtype prot"; 
				#system "~/ncbi-blast-2.2.29+/bin/blastp -in $file -db $file -out $file.blast";
				system "formatdb -p T -i $file -o T"; 
				system "blastall  -i $file -p blastp -d $file -a 24 -o $file.blast";
				system "mcxdeblast --score=b --sort=a --bcut=5 $file.blast";
				system "mcxassemble -b $file.blast -q -r max --map";
				system "mcl $file.blast.sym -scheme 7 -I 1.2 -te 20 -o $file.clus --append-log=yes";
				system "rm $file.blast $path_train/*.psq $path_train/*.psd $path_train/*.phr $path_train/*.pin $path_train/*.psi $path_train/*.err $path_train/*.map $path_train/*.raw $path_train/*.hdr";	
				$clusfile="$file.clus";
				$tabfile ="$file.blast.tab";		
				my $fpclusterObj=HMMmodE->new(-file => $file, -clusType => 'MCL', -clusFile => $clusfile, -tabFile => $tabfile);	
				my $ref=$fpclusterObj->clusInfo();			
				%hash_info=%$ref;														
				print "fp object....\n";			
				print Dumper (\%hash_info);
				@clus=keys %hash_info;
				my $tot_fp_seq=scalar @fp_ids;
				@fp_final=();
				
				## RANDOM SAMPLING 				
				foreach my $clusno(@clus)							{
					my @temp=();
					my $percent_cluster=$hash_info{$clusno}/$tot_fp_seq;
					my $seq_no=int ($percent_cluster*200);
					my $ref_fp_clus_seq=$fpclusterObj->{'cluster_seq'}{$clusno};
					my @fp_clus_seq=@$ref_fp_clus_seq;
					my @rand_array=_gen_rand_no($seq_no,($#fp_clus_seq+1));
					my @uniq=_remove_dupes(($#fp_clus_seq+1),@rand_array);
					@temp=map {$fp_clus_seq[$_]}@uniq;
					push (@fp_final,@temp);
				}
				my $out_fp_samp = Bio::SeqIO->new('-file' => ">$path_train/neg.samp.seq", '-format' => 'Fasta');
				$out_fp_samp->write_seq(@fp_final);
				$fpseqfile="$path_train/neg.samp.seq";
				print scalar @fp_final, "fps found;\nRunning threshold optimisation and ModE...\n";
			}
			else
			{
				@fp_final=@fp_seq;
				$fpseqfile="$path_train/neg.seq";	
				print scalar @fp_final, "fps found;\nRunning threshold optimisation and ModE...\n";
			}
			
			## 10 fold Cross validation
			$clusterObj->generateSamples($clus_no,\@fp_final);	
			#exit;
			my @cutoff=$clusterObj->crossValidate($clus_no);
			#exit;
				
 			#print "Array of cutoff scores" ,"\t", "@cutoff\n";
 			print OUT "$path_train/$outfile.hmm\t$cutoff[0]\n";
	
			# here begins the modification of the HMM profile....HMM-ModE 	
			
			# first..make an alignment of the 'False Positive' sequences
			system ("muscle -in  $fpseqfile -maxiters 2 -maxmb 4200 -quiet -out $path_train/tot_fp.aln");	
			
			# next ..build an HMM from the alignment
			system ("hmmbuild --informat afa $path_train/tot_fp.hmm $path_train/tot_fp.aln");
			
			# ALign the 'False Positive' alignment with the 'True Positive' one 			
			my $aligned=$clusterObj->mAlign("muscle",\@params,$clus_no,"$path_train/$outfile.aln","$path_train/tot_fp.aln","p"); 
			
			# create a Bio::DB object from the combined alignment file to enable separation of the 'profile aligned' 'True Positive' & fp sequences
			my $db = new Bio::DB::Fasta($aligned);  
			
			# SeqIO object for 'True Positive' alignment
			my $tp_out=Bio::SeqIO->new('-file' => ">$path_train/$outfile.tp.prof.aln",'-format' => 'Fasta');
			foreach my $id (@own_ids)
			{
				my $seqobj = $db->get_Seq_by_id($id); # get a PrimarySeq obj
				if($seqobj)
				{
					$tp_out->write_seq($seqobj);  # Fetch the 'True Positive' ids out from the profile-alignment and write it out as FASTA
				}
			}
			
			# SeqIO object for 'False Positive' alignment
			my $fp_out=Bio::SeqIO->new('-file' => ">$path_train/$outfile.fp.prof.aln",'-format' => 'Fasta');
			foreach my $id (@fp_ids)
			{
				my $seqobj = $db->get_Seq_by_id($id); # get a PrimarySeq obj
				if($seqobj)
				{
					$fp_out->write_seq($seqobj);  # Fetch the 'False Positive' ids out from the profile-alignment and write it out as FASTA
				}
			}
			
			# Build the HMM from the profile-aligned 'True positive' alignment
			system ("hmmbuild --informat afa $path_train/$outfile.tp.prof.hmm $path_train/$outfile.tp.prof.aln");
			
			# Build the HMM from the profile-aligned 'False positive' alignment
			system ("hmmbuild --informat afa $path_train/$outfile.fp.prof.hmm $path_train/$outfile.fp.prof.aln");
			
			# Modify the 'True positive' profile emission probabilites -ModE
			my $hmmmodefile=$clusterObj->modE($clus_no,"-","$path_train/$outfile.tp.prof.hmm","$path_train/$outfile.fp.prof.hmm");
			
			# Writes out the HMM-ModE threshold
			print OUT "$hmmmodefile\t$cutoff[1]\n";	
			}
			
			else
			{
				print OUT "$path_train/$outfile.hmm\t0\n";
				print "No fps...use Hmmer profile with default cut-off\n";
			}
			$clusterObj->cleanUp();
		}	
	}
}		
sub remove_dupes
{
	print "Removing duplicate sequences if any\n";
	my $file_name=shift;
	my %hash=();
	my $key="";
	my $old_key="";
	my $seq="";
	open(IN,"$file_name");
	while (my $line=<IN>)
	{
		chomp $line;
        	if ($line =~ /\>/)
        	{
        		$key=$line;
        		$key=~s/\>//;
        		$hash{$old_key}=$seq;
        		$seq ="";
        		$old_key=$key;
        	}
        	else
        	{
        		$seq .= $line;
        	}
	$hash{$old_key}=$seq;
	}
	my @keys=sort keys %hash;
	shift @keys;
	close IN;
	open (OUT,">$file_name");
	foreach $key(@keys)
	{
        	print OUT ">$key\n$hash{$key}\n";
	}
	close (OUT);
}
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
END { warn time - $^T, " seconds elapsed\n" }		
