#!/usr/bin/python

from optparse import OptionParser, OptionGroup
import sys
import os
import shutil

import pipeline_core

def setup_controls(options, nucl_controls, prot_controls):
	if options.gene_name == "16s" or options.gene_name == "16s":
		options.is_prot = False

		if nucl_controls:
			options.is_control = True
			options.nucl_controls = os.path.abspath(nucl_controls)
		else:
			options.is_control = False
	else:
		options.is_prot = True

		if nucl_controls and prot_controls:
			options.is_control = True
			options.nucl_controls = os.path.abspath(nucl_controls)
			options.prot_controls = os.path.abspath(prot_controls)
		elif (nucl_controls and not prot_controls) or (not nucl_controls and prot_controls):
			raise ValueError("Must specify both protein and nucleotide control sequences for non 16S genes %s" % options.gene_name)	
		else:
			options.is_control = False

		options.framebot_ref = os.path.join(pipeline_core.resources_dir, options.gene_name + "/framebot.fasta")

def run_init_process(options, seq_files, trace):
	return pipeline_core.init_process(seq_files, options.tag_file, options.forward_primers, options.fedit, options.reverse_primers, options.redit, options.gene_name, options.min_length, options.min_qual, options.max_ns, options.process_notag, options.keep_primers, options.workdir, trace)

def run_pipeline(options, seq_files, run_init_process, trace):
	for seq_file in seq_files:
		seq_file.seq_file = os.path.abspath(seq_file.seq_file)
		if seq_file.qual_file:
			seq_file.qual_file = os.path.abspath(seq_file.qual_file)
	
	if len(seq_files) == 0:
		print "No files found, quitting"
		return

	print "Starting pipeline with files", seq_file

	if run_init_process:
		print "Running initial processor"
		seq_files = run_init_process(options, seq_files, trace)

	if options.is_control:
#		pipeline_core.chimera_check(options, seq_files, trace)
		print "Running error calculator"
		pipeline_core.error_calc(seq_files, options.nucl_controls, options.workdir, trace)
	
	print "Dereplicating"
	seq_files = pipeline_core.dereplicate(seq_files, workdir = options.workdir, trace = trace)
	derep_file = seq_files[0]

	if not options.is_prot and options.is_control:
		print "Removing contaminants"
		seq_files = pipeline_core.rrna16s_decontamination(seq_files, options.nucl_controls, options.decontamination_cutoff, options.workdir, trace)
#	elif options.is_prot and options.is_control:
#		blast_decontamination(options, seq_files, trace)

	if options.is_prot:
		print "Running framebot"
		seq_files = pipeline_core.framebot(seq_files, options.prot_controls, options.framebot_ref_ident, options.framebot_minlength, "glocal", False, options.workdir, trace)

	print "Aligning"
	seq_files = pipeline_core.align(seq_files, options.gene_name, options.workdir, trace)
	if options.chop:
		print "I will not be chopping, it has been disabled til further notice"
#		seq_files = simple_trimmer(options, seq_files, trace)
	
	print "Refreshing mappings"
	seq_files = pipeline_core.refresh_mappings(seq_files, "updated_mappings", options.workdir, trace)

	print "Computing distances"
	matrix_files = pipeline_core.distance_matrix(seq_files, not options.is_prot, "#=GC_RF", "0.15", options.workdir, trace)
	print "Clustering"
	clust_files = pipeline_core.cluster(matrix_files, workdir=options.workdir, trace=trace)

	print "Finding representative sequences"
	pipeline_core.rep_seqs(clust_files, seq_files[0], 0.03, "#=GC_RF", "rep_seqs_0.03", workdir=options.workdir, trace=trace)

	if len(open(clust_files[0].seq_file).readline().split(":")[1].split()) > 1:	
		print "Computing jaccard and sorensen"
		pipeline_core.jaccard_sorensen(clust_files, workdir=options.workdir, trace=trace)
	
	print "Computing shannon chao"
	pipeline_core.shannon_chao1(clust_files, workdir=options.workdir, trace=trace)
	print "Running rarefaction"
	pipeline_core.rarefaction(clust_files, workdir=options.workdir, trace=trace)

	print "Exploding mappings"
	exploded_seq_file = pipeline_core.explode_mappings(seq_files, "final_seqs", options.workdir, trace)
	print "Dealigning sequences"
	pipeline_core.dealign_seqs(exploded_seq_file, trace)

	if options.is_prot:
		print "Exploding nucl sequences"
		pipeline_core.explode_mappings([derep_file], "final_seqs_nucl", options.workdir, trace)

def main(args):
	usage="usage: %prog [options] gene_name forward_primers tag_file sequence_file"

	parser = OptionParser(usage=usage)
	parser.add_option("-f", "--force", dest="force", default=False, help="Force run, deleting previous working directory if it exists", action="store_true")
	parser.add_option("-w", "--workdir", dest="workdir", default="workdir", help="Store files in this directory")
	parser.add_option("-q", "--qual-file", dest="qual_file", help="Quality file")
	parser.add_option("", "--prot-control=", dest="prot_controls", help="Unaligned protein control sequences")
	parser.add_option("", "--nucl-control=", dest="nucl_controls", help="Unaligned nucleotide control sequences")

	init_process_opts = OptionGroup(parser, "Initial Processing", "Options used during the initial processing phase")

	init_process_opts.add_option("", "--fedit", dest="fedit", default="2", help="Maximum edit distance to forward primer")
	init_process_opts.add_option("", "--redit", dest="redit", default="0", help="Maximum edit distance to reverse primer")
	init_process_opts.add_option("-r", "--reverse-primers", dest="reverse_primers", help="Reverse primers")
	init_process_opts.add_option("", "--max-ns", dest="max_ns", default="0", help="Maxmimum number of Ns allowed in a sequence")
	init_process_opts.add_option("-l", "--min-length", dest="min_length", default="0", help="Minimum length of a sequence to pass initial processing")
	init_process_opts.add_option("", "--keep-primers", dest="keep_primers", default=False, help="Don't remove primers during initial processing", action="store_true")
	init_process_opts.add_option("", "--min-qual", dest="min_qual", default="20", help="Minimum average exponential quality score")
	init_process_opts.add_option("", "--process-notag", dest="process_notag", default=True, help="Controls whether sequences that don't match a tag go through initial processing")

	parser.add_option_group(init_process_opts)

	framebot_opts = OptionGroup(parser, "Framebot Processing", "Options used during the framebot phase")
	framebot_opts.add_option("", "--framebot-min-ident", dest="framebot_ref_ident", default=".3", help="Minimum identity to reference sequence")
	framebot_opts.add_option("", "--framebot-minlength", dest="framebot_minlength", default="100", help="Minimum length to pass")

	parser.add_option_group(framebot_opts)

	decontamination_opts = OptionGroup(parser, "Decontamination", "Options used when detecting contaminants")
	decontamination_opts.add_option("", "--decontamination-cutoff=", dest="decontamination_cutoff", help="Minimum difference in seqmatch scores to be considered a contaminant", default="0.20")

	parser.add_option_group(decontamination_opts)

	options, args = parser.parse_args()

	if len(args) != 4:
		parser.error("Incorrect number of arguments")

	options.workdir = os.path.abspath(options.workdir)
	options.gene_name = args[0]
	options.forward_primers = args[1]
	options.tag_file = os.path.abspath(args[2])

	setup_controls(options, options.nucl_controls, options.prot_controls)

	if options.qual_file:
		options.qual_file = os.path.abspath(options.qual_file)

	seq_files = [SequenceFile(os.path.abspath(args[3]), options.qual_file)]

	if os.path.exists(options.workdir):
		print "Good heavens...the workdir exists"
		sys.exit(1)

	os.mkdir(options.workdir)

	out = open(os.path.join(options.workdir, "cmd.txt"), "w")
	out.write(" ".join(sys.argv) + "\n")
	out.close()

	trace = open(os.path.join(options.workdir, "trace.txt"), "w")

	seq_files = init_process(options, seq_files, trace)
	run_pipeline(options, seq_files, True, trace)

if __name__ == "__main__":
	main(sys.argv[1:])
