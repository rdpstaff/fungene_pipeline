#!/usr/bin/python

import sys
import pipeline
import ConfigParser
import os
import copy

def options_from_file(filename):
	config = ConfigParser.ConfigParser()
	config.read(filename)

	if not config.has_section("general"):
		raise IOError("Option file must have a general section")
	if not config.has_option("general", "genes_of_interest"):
		raise IOError("general section must have genes_of_interest option")
	if not config.has_option("general", "seq_file_pattern"):
		raise IOError("general section must have seq_file_pattern option")
	if not config.has_option("general", "qual_file_pattern"):
		raise IOError("general section must have qual_file_pattern option")
	
	opts = dict()
	opts["genes_of_interest"] = [t.strip() for t in config.get("general", "genes_of_interest").split(",")]
	opts["seq_file_pattern" ] = config.get("general", "seq_file_pattern")
	opts["qual_file_pattern" ] = config.get("general", "qual_file_pattern")
	
	process_notag = True
	default_fprimer = "fake_primer"
	default_rprimer = None
	max_ns = "0"
	
	if config.has_option("general", "default_fprimer"):
		default_fprimer = config.get("general", "default_fprimer").replace(" ", "")
		config.remove_option("general", "default_fprimer")
	if config.has_option("general", "default_rprimer"):
		default_rprimer = config.get("general", "default_rprimer").replace(" ", "")
		config.remove_option("general", "default_rprimer")
	if config.has_option("general", "process_notag"):
		process_notag = config.get("general", "process_notag") == "true"
		config.remove_option("general", "process_notag")
	if config.has_option("general", "max_ns"):
		max_ns = config.get("general", "max_ns")
		config.remove_option("general", "max_ns")
		
	config.remove_option("general", "genes_of_interest")
	config.remove_option("general", "seq_file_pattern")
	config.remove_option("general", "qual_file_pattern")
	
	opts["general"] = parse_section(Opts(), config, "general")

	for section in config.sections():
		if section == "general":
			continue
			
		if not section in opts["genes_of_interest"]:
			print "WARNING: gene %s present in options but not listed as gene of interest" % section
		opts[section] = parse_section(opts["general"], config, section)
	
	for gene_of_interest in opts["genes_of_interest"]:
		if gene_of_interest not in opts:
			print "WARNING: gene %s is listed as a gene of interest but not listed in option file, default options will be used" % gene_of_interest
			opts[gene_of_interest] = parse_section(opts["general"], None, gene_of_interest)
	
	opts["general"].process_notag = process_notag
	opts["general"].forward_primers = default_fprimer
	opts["general"].reverse_primers = default_rprimer
	opts["general"].max_ns = max_ns
	return opts

def parse_section(default_opts, config, section):
	nucl_control = None
	prot_control = None
	
	ret = copy.deepcopy(default_opts)
	ret.gene_name = section
	
	if not config or not config.has_section(section):
		pipeline.setup_controls(ret, nucl_control, prot_control)
		return ret
	
	if config.has_option(section, "nucl_control"):
		nucl_control = config.get(section, "nucl_control")
		config.remove_option(section, "nucl_control")
	if config.has_option(section, "prot_control"):
		prot_control = config.get(section, "prot_control")
		config.remove_option(section, "prot_control")
	
	pipeline.setup_controls(ret, nucl_control, prot_control)
	
	if config.has_option(section, "fedit"):
		ret.fedit = config.get(section, "fedit")
		config.remove_option(section, "fedit")
	if config.has_option(section, "redit"):
		ret.redit = config.get(section, "redit")
		config.remove_option(section, "redit")
	if config.has_option(section, "min_length"):
		ret.min_length = config.get(section, "min_length")
		config.remove_option(section, "min_length")
	if config.has_option(section, "min_qual"):
		ret.min_qual = config.get(section, "min_qual")
		config.remove_option(section, "min_qual")
	if config.has_option(section, "max_ns"):
		ret.max_ns = config.get(section, "max_ns")
		config.remove_option(section, "max_ns")
	if config.has_option(section, "keep_primers"):
		ret.keep_primers = config.get(section, "keep_primers").lower() == "true"
		config.remove_option(section, "keep_primers")
	if config.has_option(section, "gene_name"):
		ret.gene_name = config.get(section, "gene_name")
		config.remove_option(section, "gene_name")		
		
	if config.has_option(section, "framebot_ref_ident"):
		ret.framebot_ref_ident = config.get(section, "framebot_ref_ident")
		config.remove_option(section, "framebot_ref_ident")
	if config.has_option(section, "framebot_minlength"):
		ret.framebot_minlength = config.get(section, "framebot_minlength")
		config.remove_option(section, "framebot_minlength")
	if config.has_option(section, "decontam_cutoff"):
		ret.decontam_cutoff = config.get(section, "decontam_cutoff")
		config.remove_option(section, "decontam_cutoff")
		
	if config.has_option(section, "chop"):
		ret.chop = config.get(section, "chop").lower() == "true"
		config.remove_option(section, "chop")
		
		if config.has_option(section, "model_chop_start"):
			ret.model_chop_start = int(config.get(section, "model_chop_start"))
			config.remove_option(section, "model_chop_start")
		else:
			ret.model_chop_end = 0
		if config.has_option(section, "model_chop_end"):
			ret.model_chop_end = int(config.get(section, "model_chop_end"))
			config.remove_option(section, "model_chop_end")
		else:
			raise IOError("Must specify model chop end if you turn on chopping")
	else:
		ret.chop = False
		
	if config.has_option(section, "use_reverse"):
		ret.use_reverse = config.get(section, "use_reverse") == "true"
		config.remove_option(section, "use_reverse")
	
	if len(config.items(section)) != 0:
		raise IOError("Unknown options %s in section %s" % (config.items(section), section))
	
	return ret

class Opts:
	def __init__(self, gene_name = "not_a_gene", fedit = 2, redit = 0, use_reverse_primer = True, keep_primers = False, max_ns = 0, min_length = 300, min_qual = 20, framebot_ref_ident = 0.3, framebot_minlength = 100, nucl_control = None, prot_control = None, decontam_cutoff = 0.1, process_notag=True):
		self.gene_name = gene_name

		self.fedit = str(fedit)
		self.redit = str(redit)
		self.min_length = str(min_length)
		self.max_ns = str(max_ns)
		self.min_qual = str(min_qual)
		self.keep_primers = keep_primers

		self.framebot_ref_ident = str(framebot_ref_ident)
		self.framebot_minlength = str(framebot_minlength)
		self.decontamination_cutoff = str(decontam_cutoff)

		self.use_reverse = use_reverse_primer
		self.process_notag = process_notag

		pipeline.setup_controls(self, nucl_control, prot_control)
