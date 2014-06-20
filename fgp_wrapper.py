#!/usr/bin/env python

import sys
import os
from pipeline_core import SequenceFile
import pipeline_core
import subprocess
from Bio import SeqIO

mail_wrapper_class = "edu/msu/cme/rdp/misc/MailerCmdWrapper"

def check_seq_counts(seq_files):
	total_seqs = 0
	for seq_file in seq_files:
		total_seqs += check_seq_counts_file(seq_file.seq_file)
				
	return total_seqs

def check_seq_counts_file(seq_file):
	total_seqs = 0
	stream = open(seq_file)
	first_char = stream.read(1)
	if first_char == "@":
		format = "fastq"
	else:
		format = "fasta"
	stream.close()

	for seq in SeqIO.parse(open(seq_file), format):
		if seq.id[0] != "#":
			total_seqs += 1
	return total_seqs

def send_mail(email_address, subject, mail_file):
	cmd = ["cafe", mail_wrapper_class, "--to_email", email_address, "--subject_name", subject, "--message_file", mail_file]
	subprocess.call(cmd)

def main(options_file, command_file, infiles):
	seq_files = []
	status_stream = None

	for infile in infiles:
		lexemes = infile.split(",")
		if len(lexemes) == 2:
			seq_files.append(SequenceFile(os.path.abspath(lexemes[0]), os.path.abspath(lexemes[1])))
		else:
			seq_files.append(SequenceFile(os.path.abspath(lexemes[0])))

	matrix_files = seq_files
	cluster_files = seq_files

	lines = [x.strip() for x in open(options_file).read().strip().split("\n")]
	if len(lines) < 7:
		raise IOError("Invalid number of line in options file")
	
	gene_name = lines[0]
	basedir = lines[1]
	workdir = lines[2]
	user_email = lines[3]
	status_file = lines[4]
	result_tar = lines[5]
	mail_file = lines[6]
	
	try:	
		status_stream = open(status_file, "w")
		if os.path.exists(workdir):
			print "Good heavens...the workdir exists"
			sys.exit(1)

		os.mkdir(workdir)

		out = open(os.path.join(basedir, "cmd.txt"), "w")
		out.write(" ".join(sys.argv) + "\n")
		out.close()

		trace = open(os.path.join(basedir, "trace.txt"), "w")

		for line in open(command_file):
			lexemes = line.strip().split()
			
			if lexemes[0] == "init_proc":
				forward_primers = lexemes[1]
				fedit = lexemes[2]
				
				min_length = lexemes[3]
				max_ns = lexemes[4]
				min_qual = lexemes[5]
	
				tag_file = lexemes[6]
	
				process_notag = lexemes[7] == "true"
				keep_primers = lexemes[8] == "true"
	
				if len(lexemes) == 11:
					reverse_primers = lexemes[9]
					redit = lexemes[10]
				else:
					reverse_primers = None
					redit = 0
				
				status_stream.write("Starting initial processing...")
				status_stream.flush()
				seq_files = pipeline_core.init_process(seq_files, tag_file, forward_primers, fedit, reverse_primers, redit, gene_name, min_length, min_qual, max_ns, process_notag, keep_primers, workdir, trace)
				status_stream.write("Done\n")
				status_stream.flush()
			elif lexemes[0] == "check_seq_count":
				min_seqs = int(lexemes[1])
				max_seqs = int(lexemes[2])
				
				status_stream.write("Checking number of sequences...")
				status_stream.flush()
				total_seqs = check_seq_counts(seq_files)
				if (total_seqs < min_seqs) or (total_seqs > max_seqs):
					status_stream.write(" Error occurred: Minimum amount of sequences allowed: %s, Maximum amount of sequences allowed: %s, Total sequences from last step: %s.\n" % (min_seqs, max_seqs, total_seqs))
					status_stream.flush()
					raise Exception()
				else:
					status_stream.write("%s seqs, Okay\n" % total_seqs)
					status_stream.flush()
					
			elif lexemes[0] == "framebot":
				alignment_mode = lexemes[1]
				min_length = lexemes[2]
				min_ident = lexemes[3]
				
				status_stream.write("Running framebot")
				if len(lexemes) == 5:
					status_stream.write(" with custom reference file")
					framebot_ref = lexemes[4]
				else:
					framebot_ref = os.path.join(pipeline_core.resources_dir, gene_name + "/framebot.idx")
					if not os.path.exists(framebot_ref):	
						framebot_ref = os.path.join(pipeline_core.resources_dir, gene_name + "/framebot.fasta")
					else:
						status_stream.write(" using metric index")
				
				status_stream.write("...")
				status_stream.flush()
				seq_files = pipeline_core.framebot(seq_files, framebot_ref, min_ident, min_length, alignment_mode, False, workdir, trace)
				status_stream.write("Done\n")
				status_stream.flush()
			elif lexemes[0] == "align":
				status_stream.write("Running HMMER3 aligner...")
				status_stream.flush()
				seq_files = pipeline_core.align(seq_files, gene_name, workdir, trace)
				
				status_stream.write("Done\n")
				status_stream.flush()
			elif lexemes[0] == "dereplicate":
				derep_mode = lexemes[1]
				mask_seq = None
				prefix = "all_seqs"

				if derep_mode == "unaligned":
					unaligned = True

					if len(lexemes) == 3:
						prefix = lexemes[2]
				elif derep_mode == "aligned":
					unaligned = False

					if len(lexemes) > 2:
						mask_seq = lexemes[2]
					if len(lexemes) > 3:
						prefix = lexemes[3]

				status_stream.write("Dereplicating %s sequences..." % derep_mode)
				status_stream.flush()
				seq_files = pipeline_core.dereplicate(seq_files, unaligned, mask_seq, prefix, workdir, trace)
				status_stream.write("done\n")
				status_stream.flush()
			elif lexemes[0] == "error_analysis":
				status_stream.write("Running control analysis...")
				status_stream.flush()
				
				mock_community_file = lexemes[1]
				seq_cnt = check_seq_counts_file(mock_community_file)
				if seq_cnt > 250:
					status_stream.write("Too many sequences in mock community: %s (max 250)\n" % seq_cnt)
					raise Exception()

				error_files = pipeline_core.error_calc(seq_files, mock_community_file, workdir, trace)
				
				status_stream.write("done\n")
				status_stream.flush()
			elif lexemes[0] == "chimera_check":
				status_stream.write("Running chimera check with uchime...")
				status_stream.flush()

				if len(seq_files) != 1:
					status_stream.write("Failed, multiple sequence files specified\n")
					raise Exception("Chimera check failed, multiple input files")

				seq_files = pipeline_core.chimera_check(seq_files, workdir, trace)
				
				status_stream.write("done\n")
				status_stream.flush()
			elif lexemes[0] == "refresh_mapping":
				status_stream.write("Refreshing id and sample mappings...")
				status_stream.flush()

				if len(lexemes) >= 3:
					if len(seq_files) == 1:
						seq_files[0].idmapping = lexemes[1]
						seq_files[0].sample_mapping = lexemes[2]
					else:
						status_stream.write("Failed, multiple sequence files specified\n")
						raise Exception()

				if len(lexemes) == 2:
					out_dir = lexemes[1]
				elif len(lexemes) == 4:
					out_dir = lexemes[3]
				else:
					out_dir = "filtered_mapping"

				seq_files = pipeline_core.refresh_mappings(seq_files, out_dir, workdir, trace)

				status_stream.write("done\n")
				status_stream.flush()
			elif lexemes[0] == "explode_mapping":
				status_stream.write("Expanding id and sample mapping...")
				status_stream.flush()

				if len(lexemes) >= 3:
					if len(seq_files) == 1:
						seq_files[0].idmapping = lexemes[1]
						seq_files[0].sample_mapping = lexemes[2]
					else:
						status_stream.write("Failed, multiple sequence files specified\n")
						raise Exception()

				if len(lexemes) == 2:
					out_dir = lexemes[1]
				elif len(lexemes) == 4:
					out_dir = lexemes[3]
				else:
					out_dir = "expanded_mappings"

				seq_files = pipeline_core.explode_mappings(seq_files, out_dir, workdir, trace)

				status_stream.write("done\n")
				status_stream.flush()
			elif lexemes[0] == "distance":
				status_stream.write("Computing distances...")
				status_stream.flush()
				
				mask_seq = None
				cutoff = lexemes[1]
				if len(lexemes) == 3:
					mask_seq = lexemes[2]

				matrix_files = pipeline_core.distance_matrix(seq_files, gene_name == "16s", mask_seq, cutoff, workdir, trace)

				status_stream.write("Done.\n")
				status_stream.flush()
			elif lexemes[0] == "cluster":
				status_stream.write("Clustering...")
				status_stream.flush()
				method = lexemes[1]
				step = lexemes[2]
				
				if len(lexemes) == 4:
					if len(matrix_files) == 1:
						matrix_files[0].sample_mapping = lexemes[3]
					else:
						status_stream.write("Failed, cannot specify sample file with multiple matrices")
						raise Exception()

				cluster_files = pipeline_core.cluster(matrix_files, method, step, workdir, trace)
				
				status_stream.write("done\n")
				status_stream.flush()
			elif lexemes[0] == "rep_seqs":
				status_stream.write("Finding representative sequences...")
				status_stream.flush()
				
				cutoff = lexemes[1]
				if len(lexemes) > 2:
					mask_seq = lexemes[2]
				else:
					mask_seq = None
					
				if len(lexemes) > 3:
					ref_seqs = SequenceFile(os.path.abspath(lexemes[3]))
				else:
					ref_seqs = seq_files[0]
				
				pipeline_core.rep_seqs(cluster_files, ref_seqs, cutoff, mask_seq, "rep_seqs_%s" % cutoff, workdir = workdir, trace = trace)
				
				status_stream.write("done\n")
				status_stream.flush()
			elif lexemes[0] == "jaccard_sorensen":
				status_stream.write("Running abundance corrected Jaccard and Sorensen...")
				status_stream.flush()

				if len(lexemes) > 1:
					max_cutoff = lexemes[1]
				else:
					max_cutoff = 0.15
				if len(lexemes) > 2:
					min_cutoff = lexemes[1]
				else:
					min_cutoff = 0.0
				
				pipeline_core.jaccard_sorensen(cluster_files, min_cutoff, max_cutoff, workdir=workdir, trace=trace)
				
				status_stream.write("done\n")
				status_stream.flush()
			elif lexemes[0] == "shannon_chao":
				status_stream.write("Running Shannon and Chao1...")
				status_stream.flush()
				
				pipeline_core.shannon_chao1(cluster_files, workdir=workdir, trace=trace)
				
				status_stream.write("done\n")
				status_stream.flush()
			elif lexemes[0] == "rarefaction":
				status_stream.write("Running Rarefaction...")
				status_stream.flush()
				
				pipeline_core.rarefaction(cluster_files, workdir=workdir, trace=trace)
				
				status_stream.write("done\n")
				status_stream.flush()
			else:
				raise Exception("Unknown command %s" % lexemes[0])
	except:
		if status_stream and not status_stream.closed:
			status_stream.write("Fatal error in pipeline\n")
			status_stream.close()
		
		error_mail_file = os.path.join(basedir, "pipeline_error.mail")
		error_mail = open(error_mail_file, "w")
		error_mail.write("There was a problem completing your FunGene Pipeline job.\n\n")
		
		if os.path.exists(status_file):
			for line in open(status_file):
				error_mail.write(line)
		
		error_mail.write("\n")
		error_mail.write("Please contact rdpstaff@msu.edu if you have any questions or need assistance in processing your data.\n\n")
		error_mail.write("RDP Staff\n")	
		error_mail.close()
		
		send_mail(user_email, "FunGene Pipeline job failed", error_mail_file)
		
		raise
	
	subprocess.call(["tar", "-czf", result_tar, "-C", os.path.split(workdir)[0], os.path.split(workdir)[1]])

	status_stream.write("Pipeline completed successfully\n")
	status_stream.close()

	send_mail(user_email, "FunGene Pipeline job complete", mail_file)
	
if __name__ == "__main__":
	if len(sys.argv) > 3:
		main(sys.argv[1], sys.argv[2], sys.argv[3:])
	else:
		print "USAGE: pipeline.py <options_file> <command_file> <seq_file>,[qual_file] ..."
