# ModEnzA Protocol

## CREATE A FOLDER AND KEEP ALL THE SCRIPTS AND MODULES INSIDE THE DIRECTORY
## PREPARE THE SEQUENCE FILES FROM UNIPROT AND ENZYME DATABASE 

( this will create "list_of_all_training_seq" along with all the sequence files)

<code> $ perl prepare_seq_files_for_modenza.pl </code>

## RUN MCL CLUSTERING 
<code> $  perl data_processing.pl list_of_all_training_seq </code>

## CREATE PROFILES (this will create "list_of_all_training_seq.list_cut.hmm" along with other files) 
<code> $ perl hmmmode.pl list_of_all_training_seq </code>

## ADD GA THREHOLD IN THE CREATED FILES
<code> $ perl add_ga_threshold.pl list_of_all_training_seq.list_cut.hmm </code>

## CONCATENATE ALL THE PROFILES
<code> $ cd list_of_all_training_seq.train_samples </code>
<code> $ cat *.ModE.txt > ModEnzA-ModE-Enzyme-2018Jan31.hmm (2018Jan31 enzyme release date, change accordingly) </code>
<code> $ cat *.hmm.txt > ModEnzA-hmm-Enzyme-2018Jan31.hmm (2018Jan31 enzyme release date, change accordingly) </code>
