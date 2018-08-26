#----------------------------------------------------------------------------
# PACKAGE : hmmer3
# PURPOSE : To encapsulate code for parsing HMMER 3 profile.
# AUTHOR  : Swati Sinha under the guidence of Dr. Andrew  Lynn
# ANKNOLEDGEMENT: To be filled
#		
#----------------------------------------------------------------------------
package Hmm3;
use strict;

## POD Documentation:

=head1 NAME

hmmer - HMMER profile parser.

=head1 SYNOPSIS

=head2 Create file object using method "new":

    use hmmer;

    $HmmerObj = hmmer->new( -file => "input profile");

=head2 Parsing the file with file object

    Example- Extraction of Emission probability

    while(my ($position, $hashEmissionProbablityPerAlphabet)= each %{$HmmerObj->matchEmission})
    {
       print "$position  \n";
       while(my ($Alphabet, $probablity)=each %{$hashEmissionProbablityPerAlphabet})
       {
          print "$Alphabet -----   $probablity\n";
       }
    }


=head1 DESCRIPTION

    To be filled in

=head1 CONSTRUCTORS

=head2 hmmer->new()

    $HmmerObj = hmmer->new( -file => "input profile");
The new() class method constructs a new hmmer object.  The
returned object can be used to retrieve or print hmmer objects. 
At present new() accepts the following parameters:

=over 4

=item -file
A file path to be opened for reading.  The usual Perl
conventions apply:

   'file'       # open file for reading

=back

=head1 OBJECT METHODS

See below for more detailed summaries.  The main methods are:

=head2 $BtoDtransitionScore = $hmmerObj-E<gt>BtoDtransitionScore()

Fetch the B-E<gt>D Transition scores from hmmer profile.


=head2 $matchEmission = $hmmerObj-E<gt>matchEmission()

Fetch the match emission from hmmer profile.


=head2 $matchEmissionMap = $hmmerObj-E<gt>matchEmissionMap()

Fetch the match emission map from hmmer profile.


=head2 $stateTransition = $hmmerObj-E<gt>stateTransition()

Fetch the state transitions from hmmer profile.

=head2 $matchEmissionProbablity = $hmmerObj-E<gt>matchEmissionProbablity()

Calculate the match emission probablity from hmmer profile.

=head2 $stateTransitions = $hmmerObj-E<gt>stateTransitions()

Fetch the state transitions from hmmer profile.

=head2 $insertEmission = $hmmerObj-E<gt>insertEmission()

Fetch the insert emission from hmmer profile.

=head2 $writeProfile = $hmmerObj-E<gt>writeProfile($hmmerObj,"filename")

Writes the profile into another file as hmmer profile.

=head2 $hmmerObj-E<gt>destroy()

Cleans up the object from the memory.

=head1 PAKAGE VARIABLES

=head2 Variables as array

    BASE = array of the symbol alphabets in the hmmer profile

=head2 Variable as string

    HMMER3/b = File format version
    NAME     = Model name
    ACC	     = Accession number
    DESC     = Description line
    LENG     = Model length
    ALPH     = Symbol alphabet
    RF       = Reference annotation flag
    CS       = Consensus structure annotation flag
    MAP      = Map annotation flag
    DATE     = Creation date
    COM      = Command line log
    NSEQ     = Sequence number
    EFFN     = Effective Sequence Number
    CKSUM    = Training alignment checksum
    GA       = Pfam gathering thresholds GA1 and GA2
    TC       = Pfam trusted cutoffs TC1 and TC2
    NC       = Pfam noise cutoffs NC1 and NC2
    STATS    = Statistical parameters needed for E-vlaue calculations
    HMM      = String of the symbol alphabets in the hmmer profile
    COMPO    = This is the first line in the main model section: these are the model’s overall average match state emission probabilities, which are
	       used as a background residue composition in the “filter null” model.

See Hmmer user guide "http://hmmer.wustl.edu/ " for more detailed summaries of each variables.


=head1 EXAMPLE

    use hmmer;

    $HmmerObj = hmmer->new( -file => "input profile");

    while(my ($position, $hashEmissionProbablityPerAlphabet)= each %{$HmmerObj->matchEmission})
    {
       print "$position  \n";
       while(my ($Alphabet, $probablity)=each %{$hashEmissionProbablityPerAlphabet})
       {
          print "$Alphabet -----   $probablity\n";
       }
    }

    foreach $nuleValue (@{$HmmerObj->{'NULE'}}) { printf("%7s",$nukeValue)} print "\n";

    foreach $symbolAlphabet (@{$HmmerObj->{'BASE'}}) { printf("%7s",$symbolAlphabet)} print "\n";



=head2 new

 Title   : new
 Usage   : $hmmerObj = hmmer->new(-file => $filename)
 Function: Creates a new hmmer profile object
 Returns : A hmmer stream initialised with the appropriate format
 Args    : Named parameters:
             -file => $filename

=cut

sub new
{
	my ($class, %param) = @_;
	my $self = bless {}, $class;
	my $file = $param{-file};
	$self->{'file'} = $file;
	my $linePointer=0;
	my $main_bool=0;
	my %items;
	my @items;
	my %matchEmission;
	my @matchEmissionLine=();
	my %matchEmissionMap;
	my @insertEmissionLine=();
	my %insertEmission;
	my @stateTransitionLine=();
	my @stateTransitions=();
	my %stateTransition;
	my $line=0;
	my $compo;
	my $x;
	my @alph=();
	my @compoScore=();
	my @firstlineScore=();
	my @secondlineScore=();
	my %prob;
    	my %prob_nule;
	open(FH,$self->{'file'});
	while(<FH>)
    	{
		if(/^HMM\s+(.*)/)
		{
	    		@alph = split(" ",$1);	    
		}
		if(/^\s+COMPO\s+/)
                {
                @compoScore = split(" ",$_);
		shift @compoScore;
			for(my $i=0;$i<scalar @alph;$i++)
			{
				$prob_nule{$alph[$i]}=exp(-$compoScore[$i]);
			}
		}
                if($_ =~ /m->m/)
		{
	        	$main_bool=1;
	        	$linePointer=0;
	    		foreach my $ele (split (" ",$_))
	    		{
				$ele =~ s/->/to/; 
				push(@stateTransitions,$ele); 
	    		}
	    		next
		}
        	if($_ =~ /\/\//) # Checking the end of hmm
		{
	        	$main_bool=0;
		}
		# Reading the main section.   
		if($main_bool)
		{
	        # Reading each line in main section in hmm 
	    	# the main section starts after the following line
	    	# m->m   m->i   m->d   i->m   i->i   d->m   d->d   
	    	# The indexing is from zero.
	    	#
	    	# The first two lines (after the line starting with COMPO) in the main model section are atypical. They contain information for the core model’s BEGIN node.
	    	# these lines is indexed as "one" and "two"
	    	# First match emission line is indexed as "three"
	    	# First Insert emission line is indexed as "four"
	    	# First state transition line is indexed as "five"
	    	# and so on...
		#
		#  Reading the COMPO line in the main section of hmm ie the first line.
		# 
			if(/^\s+COMPO\s+/)	
			{
				$linePointer=0;
				@compoScore = split(" ",$_);
				$linePointer++;
				next
			}
			# End of reading COMPO line
			#
			# Reading the next two lines which are atypical
			#
			if($linePointer==1)
			{
	   			@firstlineScore = split(" ",$_);
	  			$linePointer++;
				next
			}
			if($linePointer==2)
			{
	   			@secondlineScore = split(" ",$_);
	  			$linePointer++;
				next
			}
			# End of reading first three lines
			#
			#  Reading the match emission line in the main section of hmm
			# 
			if($linePointer==3)
			{
		    		@matchEmissionLine = split(" ",$_);
		    		for(my $j=1;$j<=scalar @alph;$j++)
		    		{	
					$matchEmission{$matchEmissionLine[0]}{$alph[$j-1]} = $matchEmissionLine[$j];
					$line = $matchEmissionLine[0];
                        		# Calculating emission probablity
					$prob{$matchEmissionLine[0]}{$alph[$j-1]} = exp(-$matchEmissionLine[$j]); 
		    		}
                    		$matchEmissionMap{$line} = $matchEmissionLine[21];
                   		$linePointer++;
                    		next
			}
			# End of reading match emission line
			#
			#  Reading the insert emission line in the main section of hmm
			# 
			if($linePointer==4)
			{
		    		@insertEmissionLine = split(" ",$_);		    
		    		for(my $j=1;$j<=scalar @alph;$j++)
		    		{
					$insertEmission{$line}{$alph[$j-1]} = $insertEmissionLine[$j-1];
		    		}
		    		$linePointer++;
                    		next
			}
			# End of reading insert emission line
			#
			#  Reading the state transition line in the main section of hmm
			#
			if($linePointer==5)
			{
		    		@stateTransitionLine = split(" ",$_);
		    		for(my $j=1;$j<=scalar @stateTransitions;$j++)
		    		{
					$stateTransition{$line}{$stateTransitions[$j-1]} = $stateTransitionLine[$j-1];	        
		    		}
		    		$linePointer=3;
                    		next
			}
			# End of reading state transition line
		} 
		# End of if($bool_main)
		if(/^(\S+)\s+(.*)/)
		{
         		push(@items,$_);
    			# if {$1}=='STATS'{$item{$1}
         		$items{$1}=$2;
		}
    	}
	# End of while(<FH>)
	# End of parsing hmmer profile file
	
	$self->{'compoScore'} = \@compoScore;
	$self->{'firstlineScore'} = \@firstlineScore;
	$self->{'secondlineScore'} = \@secondlineScore;
	$self->{'matchEmission'} = \%matchEmission;
	$self->{'matchEmissionMap'} = \%matchEmissionMap;
	$self->{'insertEmission'} = \%insertEmission;
	$self->{'stateTransition'} = \%stateTransition;
	$self->{'matchEmissionProbablity'} = \%prob;
	$self->{'stateTransitions'} = \@stateTransitions;
	$self->{'nuleProbablity'} = \%prob_nule;
	# Assigning the header parameters from hmmer profile file
    	while(my ($k, $v)= each %items)
    	{
        	$self->{$k}= $v;
    	}	
    	$self->{'header'} = \@items;
    	$self->{'BASE'}=\@alph;
    	return $self;
}
# End of new

=head2 destroy

 Title   : destroy
 Usage   : $hmmerObj->destroy()
 Function: Release the memory

=cut

sub destroy 
{
   my $self = shift;
   printf("Cleaning object \"$self\" from memory at %s\n", scalar localtime);
}

=head2 compoScore

 Title   : compoScore
 Usage   : $hmmerObj->Compo
 Function: Fetch the Compo line scores
 Returns : Array with the scores in the COMPO line
 Args    : None

=cut

sub compoScore 
{
	my ($self) = @_;
	return $self->{'compoScore'};
}

=head2 firstlineScore

 Title   : firstlineScore
 Usage   : $hmmerObj->firstlineScore
 Function: Fetch the scores of the firstline
 Returns : Array with the scores in the first line of the model
 Args    : None

=cut

sub firstLine
{
	my ($self)=@_;
	return $self->{'firstlineScore'};
}

=head2 secondlineScore

 Title   : secondlineScore
 Usage   : $hmmerObj->secondlineScore
 Function: Fetch the scores of the secondline
 Returns : Array with the scores in the first line of the model
 Args    : None

=cut

sub secondlineScore
{
	my ($self)=@_;
	return $self->{'secondlineScore'};
}

=head2 nuleProbablity

 Title   : nuleProbablity
 Usage   : $hmmerObj->nuleProbablity
 Function: Convert the nuke score to nule probablity per symbol
           The probability is calculated as 
           probability = (1/20)*(2**(nule_score/INTSCALE)) where  INTSCALE = 1000
 Returns : hash - the key of the hash is "symbol" and the value is nule probability
 Args    : None

=cut

sub nuleProbablity
{
    my ($self) = @_;
    return $self->{'nuleProbablity'};
}

=head2 matchEmission

 Title   : matchEmission
 Usage   : $hmmerObj->matchEmission
 Function: Fetch the match emission score per symbol
 Returns : nested hash - the key of the hash is "node" and the key of the nested hash is symbol and the value is match emission scores
 Args    : None

=cut


sub matchEmission {
    my ($self) = @_;
    return $self->{'matchEmission'};
}

=head2 matchEmissionMap

 Title   : matchEmissionMap
 Usage   : $hmmerObj->matchEmissionMap
 Function: If MAP was yes, then there is one more number on this line, representing the alignment column index for this match state
Fetch the match emission score per symbol
 Returns : nested hash - the key of the hash is "node" and the key of the nested hash is symbol and the value is scores
 Args    : None

=cut

sub matchEmissionMap {
    my ($self) = @_;
    return $self->{'matchEmissionMap'};
}

=head2 insertEmission

 Title   : insertEmission
 Usage   : $hmmerObj->insertEmission
 Function: Fetch the insert emission per symbol from the hmmer profile
 Returns : nested hash - the key of the hash is "node" and the key of the nested hash is symbol and the value is insert emission scores
 Args    : None

=cut

sub insertEmission {
   my ($self) = @_;
   return $self->{'insertEmission'};
}

=head2 stateTransitions

 Title   : stateTransitions
 Usage   : $hmmerObj->stateTransitions
 Function: Fetch the state transitions from the hmmer profile
 Returns : Array - the elements are m->m   m->i   m->d   i->m   i->i   d->m   d->d   b->m   m->e
 Args    : None

=cut


sub stateTransitions {
   my ($self) = @_;
   return $self->{'stateTransitions'};
}

=head2 stateTransition

 Title   : stateTransition
 Usage   : $hmmerObj->stateTransition
 Function: Fetch the state transition scores from the hmmer profile
 Returns : nested hash -- the key of the hash is "node" and each element of the array "$hmmerObj->stateTransitions" is the key of the nested hash and the values corresponding values is the transition score
 Args    : None

=cut

sub stateTransition {
    my ($self) = @_;
    return $self->{'stateTransition'};
}

=head2 matchEmissionProbablity

 Title   : matchEmissionProbablity
 Usage   : $hmmerObj->matchEmissionProbablity
 Function: Convert the match emission score to match emission probablity per symbol
	   The probability is calculated as 
	   probability =antitlog(score)
 Returns : nested hash - the key of the hash is "node" and the key of the nested hash is symbol and the value is match emission probability
 Args    : None

=cut


sub matchEmissionProbablity 
{
    my ($self) = @_;
    return $self->{'matchEmissionProbablity'};
}

=head2 writeProfile

 Title   : writeProfile
 Usage   : $hmmerObj->writeProfile
 Function: writes the profile into given file
 Args    : object_of_the_profile and file_name_to_write_in

=cut

sub writeProfile
{
my ($self, $file) = @_;
open(FH_OUT,">$file");
foreach(@{$self->{'header'}})
{
	unless($_ =~ /\/\//){print FH_OUT $_};
}
	#########################
	print FH_OUT spacer(7);
	my @stateTransition = @{$self->{'stateTransitions'}};
	foreach (@stateTransition)
	{
		$_ =~ s/to/->/;
		printf FH_OUT "%9s",$_;
	}
	print FH_OUT "\n";
	##########################
	print FH_OUT spacer(2);
        my @compoScore = @{$self->{'compoScore'}};
	my $compo=shift(@compoScore);
	print FH_OUT $compo;
	print FH_OUT spacer(1);
	foreach (@compoScore)
	{
		printf FH_OUT "%9s",$_;
	}
	print FH_OUT "\n";
	##########################
	print FH_OUT spacer(8);
	foreach(@{$self->{'firstlineScore'}})
	{
		printf FH_OUT "%9s",$_;
	}
	print FH_OUT "\n";
	###########################
	print FH_OUT spacer(8);
	foreach(@{$self->{'secondlineScore'}})
	{
		printf FH_OUT "%9s",$_;
	}
	print FH_OUT "\n";
	###########################
	my @amino=("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y");
	my @columns= sort { $a <=> $b } keys %{$self->matchEmission};
	foreach my $col (@columns)
	{
		printf FH_OUT "%7s",$col;
		print FH_OUT spacer(1);
		foreach(@amino)
		{
			printf FH_OUT spacer(2);
			printf FH_OUT "%.5f",${$self->matchEmission}{$col}{$_};# round off the score to 5 digits precision
		}
		printf FH_OUT "%7s",${$self->matchEmissionMap}{$col};
		printf FH_OUT "%2s","-" if($self->{'RF'} eq "no");
		printf FH_OUT "%2s","-" if($self->{'CS'} eq "no");
		print FH_OUT "\n";
		print FH_OUT spacer(8);
		foreach(@amino)
		{
			printf FH_OUT "%9s",${$self->insertEmission}{$col}{$_};
		}
		print FH_OUT "\n";
		print FH_OUT spacer(8);             
		foreach (@{$self->{'stateTransitions'}})
		{
		$_ =~ s/->/to/;
		printf FH_OUT "%9s",${$self->stateTransition}{$col}{$_};
		}
		print FH_OUT "\n";
	}
	print FH_OUT "//\n";
	sub spacer()
	{
		my $space = shift;
		for(my $i=0;$i<$space;$i++)
		{
			print FH_OUT " ";
		}
	}
}
1;
