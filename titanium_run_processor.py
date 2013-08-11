#!/usr/bin/python

import sys
import os
import re
import copy
import threading
import shutil
from Bio import SeqIO
import ConfigParser
sys.path.append("/work/fishjord/other_projects/new_fgp_scripts")
import pipeline
from pipeline_core import SequenceFile
import pipeline_core
import titanium_options

seq_file_pattern = ".TCA.454Reads.fna"
qual_file_pattern = ".TCA.454Reads.qual"

sample_clean_regex = re.compile("[^A-Za-z0-9\_\-\.]")

class PipelineThread ( threading.Thread ):
	def __init__(self, options, in_seqfiles, trace_stream):
		threading.Thread.__init__(self)
		self.options = copy.copy(options)
		self.in_seq_files = in_seqfiles
		self.trace_stream = trace_stream
	
	def run(self):
		pipeline.run_pipeline(self.options, self.in_seq_files, False, self.trace_stream)
		print self.options.gene_name, "processing completed"
		self.trace_stream.close()

class TitaniumDataParser:
	def __init__(self, data_stream):
		self.stream = data_stream

		self.headers = list()
		for header in  self.stream.readline().strip().split("\t"):
			self.headers.append(header.lower())

	def read_next(self):
		while True:
			line = self.stream.readline()

			if line == "":
				return None
			elif line.strip() != "":
				break

		lexemes = line.strip().split("\t")
		if len(lexemes) > len(self.headers):
			raise IOError("Line " + line + " has too many fields")

		ret = dict()

		for i in range(len(lexemes)):
			lexeme = lexemes[i].lower()
			if self.headers[i] == "bar code name":
				ret[self.headers[i]] = lexeme.replace(" ", "")
			else:
				ret[self.headers[i]] = lexeme

		ret["control"] = "control" in ret and ret["control"] == "yes"

		return ret

	def __iter__(self):
		return self

	def next(self):
		ret = self.read_next()

		if ret == None:
			raise StopIteration
		else:
			return ret

	def close():
		self.stream.close()

def relative_to_abs(path):
	if "PWD" in os.environ:
		return os.path.join(os.environ["PWD"], path)
	else:
		return relative_to_abs(path)

def do_init_process(region, lines, trace_stream, seq_file_dir, options):
	workdir = "region_%s" % region
	if not os.path.exists(workdir):
		os.mkdir(workdir)
		init_process_opts = options["general"]

		tag_file = os.path.join(workdir, "tags.txt")
		tag_stream = open(tag_file, "w")
		for line in lines:
			gene_name = line["gene"]
			if gene_name not in options:
				gene_opts = init_process_opts
				print "INFO: gene %s has no custom options, using defaults" % gene_name
			else:
				gene_opts = options[gene_name]
				print "INFO: gene %s has custom options" % gene_name
			
			extended_opts = []
			extended_opts.append("fprimer=" + line["target sequence"].upper())
			if gene_opts.use_reverse:
				extended_opts.append("rprimer=" + line["reverse primer sequence"].upper())
			else:
				extended_opts.append("rprimer=")	
			
			if init_process_opts.fedit != gene_opts.fedit:
				extended_opts.append("fedit=" + gene_opts.fedit)
			if init_process_opts.redit != gene_opts.redit:
				extended_opts.append("redit=" + gene_opts.redit)
			if init_process_opts.min_length != gene_opts.min_length:
				extended_opts.append("min_length=" + gene_opts.min_length)
			if init_process_opts.min_qual != gene_opts.min_qual:
				extended_opts.append("min_qual=" + gene_opts.min_qual)
			if init_process_opts.max_ns != gene_opts.max_ns:
				extended_opts.append("max_ns=" + gene_opts.max_ns)

			if gene_opts.gene_name == "16s":
				extended_opts.append("gene=RRNA16S")
			else:
				extended_opts.append("gene=OTHER")
		
			tag_stream.write("%s\t%s\t%s\n" % (line["bar code sequence"], line["bar code name"], ",".join(extended_opts)))
		tag_stream.close()

		init_process_opts.workdir = workdir
		init_process_opts.tag_file = tag_file

		seq_files = [SequenceFile(os.path.join(seq_file_dir, region + options["seq_file_pattern"]), os.path.join(seq_file_dir, region + options["qual_file_pattern"]))]
		pipeline.run_init_process(init_process_opts, seq_files, trace_stream)

def split_to_samples(gene_name, lines):
	sample_to_seqs = dict()
	
	for line in lines:
		sample_name = sample_clean_regex.sub("_", line["sample name"])
		seqfile = "region_" + line["region"] + "/initial_process/" + line["bar code name"] + "/" + line["bar code name"] + "_trimmed.fasta"
		qualfile = "region_" + line["region"] + "/initial_process/" + line["bar code name"] + "/" + line["bar code name"] + "_trimmed.qual"
		
		#Since initial processing can delete empty files...we have to account for this
		if not os.path.exists(seqfile):
			continue
		
		if sample_name == "":
			sample_name = "no_sample_name"
		
		if not sample_name in sample_to_seqs:
			sample_to_seqs[sample_name] = dict()
			sample_to_seqs[sample_name]["seqs"] = []
			sample_to_seqs[sample_name]["qual"] = []
		
		sample_to_seqs[sample_name]["seqs"].extend(SeqIO.parse(open(seqfile), "fasta"))
		if os.path.exists(qualfile):
			sample_to_seqs[sample_name]["qual"].append(qualfile)
			
	ret = []
	
	for sample_name in sample_to_seqs.keys():
		seqs = sample_to_seqs[sample_name]["seqs"]
		qualfiles = sample_to_seqs[sample_name]["qual"]
		if len(seqs) > 0:
			seqfile = os.path.join(gene_name, sample_name + ".fasta")
			qualfile = None
			
			out = open(seqfile, "w")
			SeqIO.write(seqs, out, "fasta")
			out.close()
			
			if len(qualfiles) > 0:
				qualfile = relative_to_abs(os.path.join(gene_name, sample_name + ".qual"))
				pipeline_core.cat_files(qualfiles, qualfile, False)
			
			ret.append(SequenceFile(relative_to_abs(seqfile), qualfile))
	
	return ret

def process_run(data_file, custom_opts_file, seq_file_dir):
	parser = TitaniumDataParser(open(data_file))
	options = titanium_options.options_from_file(custom_opts_file)
	
	region_map = dict()
	for line in parser:
		if not line["region"] in region_map:
			region_map[line["region"]] = list()
		region_map[line["region"]].append(line)

	gene_map = dict()
	trace_stream = open("region_trace.txt", "w")
	for region in region_map:
		lines = region_map[region]
		do_init_process(region, lines, trace_stream, seq_file_dir, options)

		for line in lines:
			gene_name = line["gene"]
			if not gene_name in  gene_map:
				gene = dict()
				gene["control"] = []
				gene["experimental"] = []
				gene_map[gene_name] = gene

			if line["control"]:
				gene_map[gene_name]["control"].append(line)
			else:
				gene_map[gene_name]["experimental"].append(line)
	trace_stream.close()
	
	for gene_name in gene_map.keys():
		if not gene_name in options["genes_of_interest"]:
			print "Not interested in %s" % gene_name
			continue
			
		if os.path.exists(gene_name):
			print "Already a directory for %s, skipping (manually delete if you want to rerun)" % gene_name
			continue
		
		print "Processing %s" % gene_name
		os.mkdir(gene_name)

		control_files = split_to_samples(gene_name, gene_map[gene_name]["control"])
		experimental_files = split_to_samples(gene_name, gene_map[gene_name]["experimental"])
		
		if len(control_files) > 0 or len(experimental_files) > 0:
			if len(control_files) > 0:
				workdir = relative_to_abs(os.path.join(gene_name, "control"))
				os.mkdir(workdir)

				trace_stream = open(os.path.join(workdir, "trace.txt"), "w")
				options[gene_name].is_control = True
				options[gene_name].workdir = workdir
				
				print "Spawning pipeline for %s control" % gene_name
				PipelineThread(options[gene_name], control_files, trace_stream).start()

			if len(experimental_files) > 0:
				workdir = relative_to_abs(os.path.join(gene_name, "experimental"))
				os.mkdir(workdir)

				trace_stream = open(os.path.join(workdir, "trace.txt"), "w")
				options[gene_name].is_control = False
				options[gene_name].workdir = workdir
				
				print "Spawning pipeline for %s experimental" % gene_name
				PipelineThread(options[gene_name], experimental_files, trace_stream).start()
		else:
			shutil.rmtree(gene_name)
		

if __name__ == "__main__":
	print os.getcwd()
	if len(sys.argv) != 4:
		print "USAGE: titanium_run_processor.py <run data file> <custom_config_ini> <run seq file directory>"
	else:
		process_run(sys.argv[1], sys.argv[2], sys.argv[3])		
#	process_run("run_data.txt", "/scratch/wangqion/qiong_titanium/titanium_run_04222010/20100422_reads")
