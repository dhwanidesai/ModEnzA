# ModEnzA Protocol for generating accurate HMM profiles for Enzymes (EC nmumbers)

## Installation

The scripts and modules in this protocol require BioPerl to be installed. YOu can get it from https://bioperl.org/

The two Perl modules HMMmodE.pm and Hmm3.pm need to be placed in the folder where the scripts are being run

Example: If you are planning to generate the profiles in a folder called ModEnzA-2018 and if your modules and scripts are in ~/Downloads/ModEnzA

<code>mkdir ModEnzA-2018 </code>

<code>cp ~/Downloads/ModEnzA/*.pm ModEnzA-2018 </code>

The Perl scripts (files with extension .pl) for processing the Enzyme data could be placed in the local bin folder (either /usr/local/bin or /home/username/bin). If the Perl scripts are placed at these locations, make sure they are executable (use chmod).

Alternatively, these scripts can also be copied to the folder where the profiles are being generated.

<code>cp ~/Downloads/ModEnzA/*.pl ModEnzA-2018 </code>

To run the scripts, you need to change directory to where you are going to run the scripts

<code>cd ModEnzA-2018 </code>
  
## Usage
The ModEnzA protocol for generating Enzyme HMM profiles is acccomplished using the following scripts. If the scrips are copied as executables to a folder included in the $PATH environment variable, you can run the scripts without typing "perl" in the beginning

### PREPARE THE SEQUENCE FILES FROM UNIPROT AND ENZYME DATABASE 

( this will create "list_of_all_training_seq" along with all the sequence files)

<code>perl prepare_seq_files_for_modenza.pl </code>

### RUN MCL CLUSTERING 
<code>perl data_processing.pl list_of_all_training_seq </code>

### CREATE PROFILES (this will create "list_of_all_training_seq.list_cut.hmm" along with other files) 
<code>perl hmmmode.pl list_of_all_training_seq </code>

### ADD GA THREHOLD IN THE CREATED FILES
<code>perl add_ga_threshold.pl list_of_all_training_seq.list_cut.hmm </code>

### CONCATENATE ALL THE PROFILES
<code>cd list_of_all_training_seq.train_samples </code>

<code>cat *.ModE.txt > ModEnzA-ModE-Enzyme-2018Jan31.hmm (2018Jan31 enzyme release date, change accordingly) </code>

<code>cat *.hmm.txt > ModEnzA-hmm-Enzyme-2018Jan31.hmm (2018Jan31 enzyme release date, change accordingly) </code>


## Citation
If you find these scripts useful, please cite:

* DK Desai, S Nandi, PK Srivastava, AM Lynn (2011) ModEnzA: accurate identification of metabolic enzymes using function specific profile HMMs with optimised discrimination threshold and modified emission probabilities
Advances in bioinformatics, Volume 2011, Article ID 743782
http://dx.doi.org/10.1155/2011/743782

For the HMMer3 implemenation of the methods please cite:

* Sinha S and Lynn AM (2014) HMM‐ModE: implementation, benchmarking and validation with HMMER3. BMC Res Notes 7, 483
https://doi.org/10.1186/1756-0500-7-483
