# Functional Gene Pipeline Scripts

Intro...

Functional Gene Pipeline Scripts contains a set of python scripts that allows to run one or more 
individual tools offered by RDP FunGene Pipeline (http://fungene.cme.msu.edu). These tools are offered a modular fashion 
allowing researchers to choose the appropriate subset based on their needs.

## Setup

### Required Software
* RDPTools (https://github.com/rdpstaff/RDPTools)
* HMMER3 (http://hmmer.janelia.org protein aligner )	
* Infernal 1.1 (http://infernal.janelia.org 16S rRNA Aligner)	
* USEARCH (http://www.drive5.com/usearch chimera check)
	
### Configuration

The path to the programs and resources are specfied in file `config.ini`. Modify the example `config_skel.ini` and rename it as `config.ini`.

<a name="Usage"></a>
### Usage
The directory `examplefiles` contains the following files:

* sample 454 sequencing data files `1.TCA.454Reads.fna` and `1.TCA.454Reads.qual`. 
* sample nifH nucleotide sequence files `nifH_s1.fa` and `nifH_s2.fa`. 
* an example run descriptor file `16S_rundata.txt`, a tab-delimited file containing a list of sample name, and forward and reverse primers 
for each barcode in each region. Note the header line must kept unchanged. 
* an example gene option file `gene_option.txt`, specifies the parameters used for each program if different from the default. 
The section "[general]" specifies parameters used for all genes of interest. 
The section starts with "[gene]" can be used to override default parameters and those specified in the "[general]" section for that gene. 
* example command option files `framebot_options.txt` and `16s_options.txt`. The command option file should contain the following information in order, one item per line:
`gene_name, basedir, workdir, user_email, status_file, result_tar and mail_file`.
* an example command file `framebot_commands.txt`, `framebot_cluster_commands.txt`. The command file lists the command to be executed in order. 
Each line contains a tab-separated list containing a command name and the parameters used for that command. 
See the complete list of commands in `fgp_wrapper.py`.


#### Pipeline with command and options file

* Here is an example command to run Initial Process, Alignment and Clustering using `fgp_wrapper.py`. 
Modify the file `16S_options.txt` with the correct path to the output directory on your local machine and your email address before run the command:
	
		/path/to/fungene_pipeline/fgp_wrapper.py /path/to/fungene_pipeline/examplefiles/16S_options.txt /path/to/fungene_pipeline/examplefiles/init_cluster_commands.txt 
		/path/to/fungene_pipeline/examplefiles/1.TCA.454Reads.fna,/path/to/fungene_pipeline/examplefiles/1.TCA.454Reads.qual
	
* An example command to run FrameBot using `fgp_wrapper.py`. Modify the file `framebot_options.txt` with the correct path 
to the output directory on your local machine and your email address before run the command:

		/path/to/fungene_pipeline/fgp_wrapper.py /path/to/fungene_pipeline/examplefiles/framebot_options.txt /path/to/fungene_pipeline/examplefiles/framebot_commands.txt  
		/path/to/fungene_pipeline/examplefiles/nifH_s1.fa  /path/to/fungene_pipeline/examplefiles/nifH_s2.fa 

* Another example command that runs FrameBot and Clustering using `fgp_wrapper.py`. 

		/path/tofungene_pipeline/fgp_wrapper.py /path/tofungene_pipeline/examplefiles/framebot_options.txt /path/tofungene_pipeline/examplefiles/framebot_cluster_commands.txt 
		/path/tofungene_pipeline/examplefiles/nifH_s1.fa /path/tofungene_pipeline/examplefiles/nifH_s2.fa

#### With run descriptor file
* As mentioned above, the run descriptor file allows you to specify mapping for multiple barcodes at different regions in an entire sequencing run. 
Here is an example command that first runs Initial Process which splits and filters the samples by barcode, and then for each barcoded-file, runs
the list of commands in `pipeline.py` [described below](#pipeline_py):

    	/path/to/fungene_pipeline/titanium_run_processor.py /path/to/fungene_pipeline/16S_454examplefiles/16S_rundata.txt 
    	/path/to/fungene_pipeline/16S_454examplefiles/gene_option.txt /path/to/fungene_pipeline/16S_454examplefiles/

#### Interactive Python



### What's include in the Package
* `config_skel.ini`

	Configuration file containing path of the program and resource files.
* `resources`

	Contains the alignment models and reference database for each gene
* `fgp_wrapper.py`

	Runs commands listed in a command file for one input sequence file.
* `parseErrorAnalysis.py`

	Generates summary output file for the Defined Community Analysis tool.
<a name="pipeline_py"></a>	
* `pipeline.py`

	Takes one sequence file (and optional quality file), runs the following commands in order, Initial Process, Chimera Check, 
	Defined Community Analysis (mock control sample only), FrameBot (protein reads only), Alignment, Clustering, 
	Diversity Analysis(Shanno and Chao), Rarefaction. 

* `pipeline_core.py`

	Contains a set of basic functions to run each individual commands called by other python scripts
	
* `pipeline_wrapper.sh`

	A shell script wrapper that runs `fgp_wrapper.py`
* `seq_trimmer_model.py`

	Slices a set of alignment sequences.
	
* `titanium_options.py`
	
	Parses a gene option file. Used by other scripts.
* `titanium_run_processor.py`

	Takes a directory of sequence data files and a run descriptor file. Runs `Initial Process` and a listed of commands called by `pipeline.py`.
* `examplefiles`

	Contains the sample data files and example command and option files used by the [Usage example](#Usage).	


