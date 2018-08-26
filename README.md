# ModEnzA Protocol

## ! CREATE A FOLDER AND KEEP ALL THE SCRIPTS AND MODULES INSIDE THE DIRECTORY
## ! PREPARE THE SEQUENCE FILES FROM UNIPROT AND ENZYME DATABASE 

( this will create "list_of_all_training_seq" along with all the sequence files)

<code> $ perl prepare_seq_files_for_modenza.pl <\code>

! RUN MCL CLUSTERING 
$  perl data_processing.pl list_of_all_training_seq

! CREATE PROFILES (this will create "list_of_all_training_seq.list_cut.hmm" along with other files) 
$ perl hmmmode.pl list_of_all_training_seq

! ADD GA THREHOLD IN THE CREATED FILES
$ perl add_ga_threshold.pl list_of_all_training_seq.list_cut.hmm

!CONCATENATE ALL THE PROFILES
$ cd list_of_all_training_seq.train_samples
$ cat *.ModE.txt > ModEnzA-ModE-Enzyme-2018Jan31.hmm (2018Jan31 enzyme release date, change accordingly)
$ cat *.hmm.txt > ModEnzA-hmm-Enzyme-2018Jan31.hmm (2018Jan31 enzyme release date, change accordingly)
