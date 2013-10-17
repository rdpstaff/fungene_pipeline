#!/usr/bin/python

import sys
import os
import shutil
import subprocess
import ConfigParser
import inspect
import glob
from Bio import SeqIO

"""
Setup preamble

Will look for a file (config.ini) in the directory containing this 
script file (will look in the right place even if imported from
a script in another directory).  

The config file should be an ini file with a pipeline section that 
contains the properties resource_dir, blastx_db, blastn_db, 
cmalign_cmd, hmmalign_cmd, blast_cmd, formatdb_cmd.

Some basic sanity checking is done to ensure the executables exist
and the classpath contains at least some jars.
"""

config_file = os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe()))[0], "config.ini"))
config_file.replace("/export/home", "/home") #stupid nfs

config = ConfigParser.ConfigParser()
config.read(config_file)

sys.stderr.write("Config file: %s\n" % config_file)

#setup paths
resources_dir = config.get("pipeline", "resource_dir")
blastx_db = config.get("pipeline", "blastx_db")
blastn_db = config.get("pipeline", "blastn_db")

#setup the commands
parse_error_analysis = config.get("pipeline", "parse_error_analysis_cmd")
cmalign = config.get("pipeline", "cmalign_cmd")
hmmalign = config.get("pipeline", "hmmalign_cmd")
blast = config.get("pipeline", "blast_cmd")
usearch = config.get("pipeline", "usearch_cmd")
formatdb = config.get("pipeline", "formatdb_cmd")

qsub_path = config.get("pipeline", "gridware_env_path")

distribute_jobs = (config.get("pipeline", "distribute_jobs") == "true")

for cmd in [cmalign, hmmalign, blast, formatdb]:
	if not os.access(cmd, os.X_OK):
		raise ValueError("%s doesn't exist or isn't executable" % cmd)

#setup the class paths to the required classes
init_process_class = "edu.msu.cme.rdp.initprocess.InitialProcessorMain"
decontamination_class = "edu.msu.cme.rdp.chimerabot.DecontaminationBot"
framebot_class = "edu/msu/cme/rdp/framebot/cli/FramebotMain"
chimerabot_class = "edu/msu/cme/rdp/chimerabot/ChimeraBot"
error_class = "edu.msu.cme.rdp.alignment.errorcheck.CompareErrorType"
align_merger_class = "edu/msu/cme/rdp/alignment/AlignmentMerger"
cluster_main_class = "edu.msu.cme.pyro.cluster.ClusterMain"
jaccard_sorensen_class = "edu/msu/cme/rdp/abundstats/cli/AbundMain"
shannon_chao_class = "edu/msu/cme/rdp/abundstats/ShannonChao"
rarefaction_class = "edu/msu/cme/rdp/rarefaction/Rarefaction"
aligner_stats_class = "edu/msu/cme/rdp/pyro/stats/AlignerStats"
cluster_stats_class = "edu/msu/cme/rdp/pyro/stats/ClusterStats"

"""
	Our classes!
"""
class SequenceFile:
	"""
	SequenceFile is a wrapper to hold the sequence and quality file
	Although to be fair a lot of non-sequence files are stored in
	seq_file (cluster results, etc) so it is a bit of a misnomer, but
	what can you do?
	"""
	
	def __init__(self, seq_file, qual_file = None, idmapping = None, sample_mapping = None):
		self.seq_file = seq_file
		self.qual_file = qual_file
		self.idmapping = idmapping
		self.sample_mapping = sample_mapping

	def __repr__(self):
		return "%s %s" % (self.seq_file, self.qual_file)

class Command:
	"""
	Class to wrap a command array, stdout target file, and stderr target file
	used for telling run_cmds what commands to run
	"""
	def __init__(self, cmd, stdout=None, stderr=None):
		self.cmd = cmd
		self.stdout = stdout
		self.stderr = stderr
		
class BlastResult:
	"""
	Given a tab (-m 8 or -m 9) blast result file line 
	(non comment line), parses out the fields
	"""
	def __init__(self, line):
		line = line.split("\t")
		if len(line) != 12:
			raise Exception("Invalid tabular blast line '" + tab_line + "'")
			
		self.qid = line[0]
		self.sid = line[1]
		self.ident = float(line[2])
		self.length = int(line[3])
		self.mismatches = int(line[4])
		self.gap_openings = int(line[5])
		self.qstart = int(line[6])
		self.qend = int(line[7])
		self.sstart = int(line[8])
		self.send = int(line[9])
		self.eval = float(line[10])
		self.bits = float(line[11])

"""
Utility functions
"""

def run_commands(cmds, trace, distributed=False, workdir=os.getcwd()):
	"""
	Runs a list of Command classes either serially in the current thread or
	submits the jobs to gridware (if distributed=True), returns when all
	jobs complete
	"""
	jids = []
	for cmd in cmds:
		if distributed and distribute_jobs:
			qsub = ["qsub", "-terse", "-b", "y", "-wd", workdir, "-v", "PATH=%s" % qsub_path]
			if cmd.stdout:
				qsub.extend(["-o", cmd.stdout])

			qsub.extend(cmd.cmd)

			trace.write(" ".join(qsub) + "\n")
			trace.flush()
			qsub_stdout = subprocess.Popen(qsub, stdout=subprocess.PIPE).communicate()[0].strip()
			#print qsub_stdout.strip()

			jids.append(qsub_stdout)
		else:
			out = sys.stdout

			trace.write(" ".join(cmd.cmd))
			if cmd.stdout != None:
				out = open(cmd.stdout, "w")
				trace.write(" > " + cmd.stdout)
			trace.write("\n")
			trace.flush()

			subprocess.call(cmd.cmd, stdout=out, cwd=workdir)

			if out != sys.stdout:
				out.close()

	if distributed and len(jids) != 0:
		qsub = ["qsub", "-terse", "-sync", "y", "-b", "y", "-wd", workdir, "-hold_jid", ",".join(jids), "echo"]
#		subprocess.check_call(qsub)
		"""So there is this annoying problem where sometimes the queue master becomes unresponsive...so basically if it exits it's done...tabun"""
		subprocess.call(qsub, stdout=subprocess.PIPE)
		trace.write(" ".join(qsub) + "\n")
		trace.flush()

def split_seq_file(f, max_seqs, workdir, suffix):
	"""
	Given a sequence file and the maximum number of sequences
	splits the input file in to n files that contain at most
	max_seqs seqs where n = ceil(seqs_in_f / max_seqs) with
	the file name suffix + fileno
	"""
	file_count = 0
	ret_files = []
	seqs = []
	for seq in SeqIO.parse(open(f), "fasta"):
		seqs.append(seq)
		if len(seqs) > max_seqs:
			out_file = os.path.join(workdir, str(file_count) + "_" + suffix)
			out = open(out_file, "w")
			SeqIO.write(seqs, out, "fasta")
			out.close()
			seqs = []
			ret_files.append(out_file)
			file_count = file_count + 1

	if len(seqs) > 0:
		out_file = os.path.join(workdir, str(file_count) + "_" + suffix)
		out = open(out_file, "w")
		SeqIO.write(seqs, out, "fasta")
		out.close()
		file_count = file_count + 1
		ret_files.append(out_file)

	return ret_files

def cat_files(in_files, out, delete_files=True):
	"""
	concatinates a list of files (in_files) and writes the output
	to out (uses the cat command)
	"""
	if len(in_files) == 0:
		raise Exception("Attemtping to cat together no files!")

	for f in in_files:
		missing_files = []
		if not os.path.exists(f):
			missing_files.append(f)

		if len(missing_files) > 0:
			raise Exception("Catting will fail, there are files that don't exist: %s" % ", ".join(missing_files))

	out_stream = open(out, "w")
	
	if len(in_files) > 1000:
		join_num = 0
		join_files = []

		for i in range(0, len(in_files), 1000):
			join_file = "%s_join%s" % (out, join_num)
			join_num += 1
			join_files.append(join_file)

			cmd = ["cat"]
			cmd.extend(in_files[i:i + 1000])

			join_stream = open(join_file, "w")
			subprocess.check_call(cmd, stdout=join_stream)
			join_stream.close()

		if delete_files:
			for f in in_files:
				os.remove(f)

		in_files = join_files

	cmd = ["cat"]
	cmd.extend(in_files)

	subprocess.call(cmd, stdout=out_stream)
	out_stream.close()
	
	if not os.path.exists(out):
		raise Exception("Catting of %s to %s failed" % (",".join(in_files), out))
		
	if delete_files:
		for f in in_files:
			os.remove(f)

def check_unique_files(seq_files):
	"""
	Given a list of files makes sure the file names
	are all unique
	
	returns true if they are unique, false otherwise
	"""
       	file_names = set()

	for seq_file in seq_files:
		file_names.add(seq_file.seq_file)

	return len(file_names) == len(seq_files)

"""
Processing functions
"""

def init_process(seq_files, tag_file, forward_primers, fedit = 2, reverse_primers = None, redit = 0, gene_name="OTHER", min_length = 150, min_qual = 0, max_ns = 0, process_notag = True, keep_primers = False, workdir = os.getcwd(), trace = sys.stderr):
	"""
	Runs initial process with the supplied parameters
	
	required args = seq_files, tag_file, forward_primers
	optional args = fedit(2), reverse_primers(None), redit(0),
					gene_name(OTHER), min_length(150), min_qual(0),
					max_ns(0), process_notag(True), keep_primers(False),
					workdir(cwd), trace(stderr)
	"""
	init_proc_dir_name = "initial_process"
	if len(seq_files) != 1:
		raise ValueError("Initial process can only handle one sequence file at this time")
	
	seq_file = seq_files[0]
	
	cmd = ["cafe", "-Xmx1g", init_process_class, "--forward-primers", forward_primers, "--min-length", min_length, "--max-ns", max_ns, "--max-forward", fedit, "--outdir", workdir, "--seq-file", seq_file.seq_file, "--tag-file", tag_file, "--result-dir-name", init_proc_dir_name, "--min-qual", min_qual]
	if reverse_primers:
		cmd.extend(["--reverse-primers", reverse_primers, "-max-reverse", redit])
	if keep_primers:
		cmd.append("--keep-primer")
	if not process_notag:
		cmd.append("--skip-notag")

	cmd.append("--gene-name")
	if gene_name == "RRNA_16S_BACTERIA" or gene_name == "RRNA_16S_ARCHAEA" or gene_name == "RRNA_28S":
		cmd.append(gene_name)
	else:
		cmd.append("OTHER")
	
	if seq_file.qual_file:
		cmd.extend(["--qual-file", seq_file.qual_file])

	run_commands([Command(cmd)], trace)
	shutil.rmtree(os.path.join(workdir, "tag_sort"))

	seq_files = []
	init_proc_dir = os.path.join(workdir, init_proc_dir_name)
	for f in os.listdir(init_proc_dir):
		if f == "NoTag":
			continue

		if os.path.isdir(os.path.join(init_proc_dir, f)):
			stem = os.path.join(os.path.join(init_proc_dir, f), f) + "_trimmed"
			if os.path.exists(stem + ".qual"):
				seq_files.append(SequenceFile(stem + ".fasta", stem + ".qual", seq_file.idmapping, seq_file.sample_mapping))
			else:
				seq_files.append(SequenceFile(stem + ".fasta", None, seq_file.idmapping, seq_file.sample_mapping))
	
	return seq_files	

def error_calc(in_seq_files, nucl_ref_file, workdir = os.getcwd(), trace = sys.stderr):
	errors_dir = os.path.join(workdir, "error_summary")
	split_dir = os.path.join(errors_dir, "splits")
	os.mkdir(errors_dir)
	os.mkdir(split_dir)

	cmds = []
	parse_error_analysis_cmds = []
	assemble_map = dict()

	for seq_file in in_seq_files:
		seq_file_name = os.path.split(seq_file.seq_file)[1].split(".")[0]
		split_in_files = split_seq_file(seq_file.seq_file, 250, split_dir, seq_file_name + ".fasta")
		
		split_align_out_files = []
		split_mismatch_out_files = []
		split_indel_out_files = []
		split_qual_out_files = []
		
		for split in split_in_files:
			split_align_out_file = os.path.join(split_dir, os.path.split(split)[1] + "_alignments.txt")
			split_mismatch_out_file = os.path.join(split_dir, os.path.split(split)[1] + "_mismatch.txt")
			split_indel_out_file = os.path.join(split_dir, os.path.split(split)[1] + "_indels.txt")
			split_qual_out_file = os.path.join(split_dir, os.path.split(split)[1] + "_qual.txt")

			split_align_out_files.append(split_align_out_file)
			split_mismatch_out_files.append(split_mismatch_out_file)
			split_indel_out_files.append(split_indel_out_file)
			split_qual_out_files.append(split_qual_out_file)

			cmd = ["cafe", "-Xmx1g", error_class, nucl_ref_file, split]
			if seq_file.qual_file:
				cmd.extend([seq_file.qual_file])
			
			cmds.append(Command(cmd))

		pairwise_file = os.path.join(errors_dir, seq_file_name + "_pairwise.aln")
		mismatch_file = os.path.join(errors_dir, seq_file_name + "_mismatch.txt")
		indels_file = os.path.join(errors_dir, seq_file_name + "_indel.txt")
		qual_file = os.path.join(errors_dir, seq_file_name + "_qual.txt")
		
		assemble_map[pairwise_file] = split_align_out_files
		assemble_map[mismatch_file] = split_mismatch_out_files
		assemble_map[indels_file] = split_indel_out_files
		
		error_summary_file = os.path.join(errors_dir, seq_file_name + "_error_summary.txt")
		
		if seq_file.qual_file:
			assemble_map[qual_file] = split_qual_out_files
			parse_error_analysis_cmds.append(Command([parse_error_analysis, "-q", qual_file, pairwise_file, mismatch_file, indels_file, nucl_ref_file], error_summary_file))
		else:
			parse_error_analysis_cmds.append(Command([parse_error_analysis, pairwise_file, mismatch_file, indels_file, nucl_ref_file], error_summary_file))

	run_commands(cmds, trace, True, split_dir)

	for k in assemble_map.keys():
		cat_files(assemble_map[k], k)
		
	run_commands(parse_error_analysis_cmds, trace, True, workdir + "/../")

	shutil.rmtree(split_dir)
	
def rrna16s_decontamination(in_seq_files, ref_seq_file, cutoff = 0.1, workdir = os.getcwd(), trace = sys.stderr):
	"""
	Runs a seqmatch 16s decontamination answering the question for each sequence "Is there a closer sequence
	in the rdp database than in my reference sequence set?" and filters the ones that are flagged as contaminants
	
	required args = in_seq_files, ref_seq_file
	optional args = cutoff(0.1), workdir(cwd), trace(stderr)
	
	"""
	if not check_unique_files(in_seq_files):
		raise Exception("I won't overwrite output files, two input files have the same name. " + ", ".join(in_seq_files))

	decontam_dir = os.path.join(workdir, "decontamination")
	split_dir = os.path.join(decontam_dir, "splits")
	os.mkdir(decontam_dir)
	os.mkdir(split_dir)

	cutoff = str(cutoff)

	out_seq_files = []
	cmds = []
	assemble_map = dict()

	for seq_file in in_seq_files:
		seq_file_name = os.path.split(seq_file.seq_file)[1].split(".")[0]

		split_in_files = split_seq_file(seq_file.seq_file, 1000, split_dir, seq_file_name + ".fasta")

		split_out_files = []
		failed_out_files = []
		stdout_files = []

		for split in split_in_files:
			split_out_file = os.path.join(split_dir, os.path.split(split)[1] + ".out")
			failed_out_file = split_out_file + ".failed"
			stdout_file = split_out_file + "_stdout.txt"

			split_out_files.append(split_out_file)
			failed_out_files.append(failed_out_file)
			stdout_files.append(stdout_file)

			cmd = ["cafe", "-Xmx1g", decontamination_class, ref_seq_file, cutoff, split, split_out_file, failed_out_file]
			cmds.append(Command(cmd, stdout_file))

		out_seq_file = os.path.join(decontam_dir, "decontam_" + seq_file_name + ".fasta")
		out_seq_files.append(SequenceFile(out_seq_file, seq_file.qual_file, seq_file.idmapping, seq_file.sample_mapping))

		assemble_map[os.path.join(decontam_dir, seq_file_name + "_failed.txt")] = failed_out_files
		assemble_map[out_seq_file] = split_out_files
		assemble_map[os.path.join(decontam_dir, seq_file_name + "_stdout.txt")] = stdout_files

	run_commands(cmds, trace, True, split_dir)

	for k in assemble_map.keys():
		cat_files(assemble_map[k], k)

	shutil.rmtree(split_dir)

	return out_seq_files

def blast_decontamination(in_seq_files, prot_controls, cutoff = 50, workdir = os.getcwd(), trace = sys.stderr):
	"""
	Flags anything as contaminants with a better -cutoff- bits saved score in the ncbi nr db using
	blastx
	
	required args = in_seq_files, prot_controls
	optional args = cutoff(50), workdir(cwd), trace(stderr)
	
	"""
	if not check_unique_files(in_seq_files):
		raise Exception("I won't overwrite output files, two input files have the same name. " + ", ".join(in_seq_files))

	blast_dir = os.path.join(workdir, "blast_deconamination")
	split_dir = os.path.join(blast_dir, "splits")
	os.mkdir(blast_dir)
	os.mkdir(split_dir)
	ref_db = os.path.join(blast_dir, "blastx_reference.fasta")

	shutil.copyfile(prot_controls, ref_db)

	run_commands([Command([formatdb, "-i", ref_db, "-p", "T", "-l", "/dev/null"])], trace)

	cmds = []
	assemble_map = dict()
	decontam_results = dict()
	
	for seq_file in in_seq_files:
		seq_file_name = os.path.split(seq_file.seq_file)[1].split(".")[0]

		split_in_files = split_seq_file(seq_file.seq_file, 250, split_dir, seq_file_name + ".fasta")
		split_out_files = []
		ref_blast_files = []

		for split in split_in_files:
			split_out_file = os.path.join(split_dir, os.path.split(split)[1] + ".out")
			ref_blast_file = split_out_file + "_refs"
			split_out_files.append(split_out_file)
			ref_blast_files.append(ref_blast_file)

			cmd = [blast, "-p", "blastx", "-d", ref_db, "-i", split, "-o", ref_blast_file, "-m", "8", "-v", "10"]
			cmds.append(Command(cmd))
			cmd = [blast, "-p", "blastx", "-d", blastx_db, "-i", split, "-o", split_out_file, "-m", "8", "-v", "10"]
			cmds.append(Command(cmd))

		assemble_map[os.path.join(blast_dir, seq_file_name + ".txt")] = split_out_files
		assemble_map[os.path.join(blast_dir, seq_file_name + "_references.txt")] = ref_blast_files
		decontam_results[os.path.join(blast_dir, seq_file_name)] = (seq_file.seq_file, os.path.join(blast_dir, seq_file_name + ".txt"), os.path.join(blast_dir, seq_file_name + "_references.txt"))

	run_commands(cmds, trace, True, split_dir)

	for k in assemble_map.keys():
		cat_files(assemble_map[k], k)

	shutil.rmtree(split_dir)
	
	out_seq_files = []
	for result_stem in decontam_results.keys():
			orig_seq_file, nr_result_file, control_blast_file = decontam_results[result_stem]
		
			control_results = dict()
			nr_results = dict()
			
			for line in open(nr_result_file):
				if line.strip() == "" or line[0] == "#":
					continue
					
				result = BlastResult(line)
				if not result.qid in nr_results:
					nr_results[result.qid] = result
			
			for line in open(control_blast_file):
				if line.strip() == "" or line[0] == "#":
					continue
					
				result = BlastResult(line)
				if not result.qid in control_results:
					control_results[result.qid] = result
			
				out = open(result_stem + "_summary.txt", "w")
				passed_seqs = []
				for seq in SeqIO.parse(open(orig_seq_file), "fasta"):
					seqid = str(seq.id)
								
					"""If it didn't hit a control sequence it's gotta be pretty bad..."""
					if seqid in control_results:
						control_result = control_results[seqid]
						if seqid in nr_results:
							nr_result = nr_results[seqid]
							diff = control_result.bits - nr_result.bits
				
							out.write("%s\t%s\t%f\t%f\t%s\t%f\t%f\t%f\n" % (seqid, control_result.sid, control_result.ident, control_result.bits, nr_result.sid, nr_result.ident, nr_result.bits, diff))
				
							if diff < cutoff:
								passed_seqs.append(seq)
						else:
							out.write("%s\t%s\t%f\t%f\tNo NR Hit\n" % (seqid, control_result.sid, control_result.ident, control_result.bits))
							passed_seqs.append(seq)
					else:
						out.write("%s\tNo control hit\n" % seqid)
				out.close()
						
				out = open(result_stem + "_passed.fasta", "w")
				SeqIO.write(passed_seqs, out, "fasta")
				out.close()
			
			out_seq_files.append(SequenceFile(result_stem + "_passed.fasta", seq_file.qual_file, seq_file.idmapping, seq_file.sample_mapping))
					
	return out_seq_files

def framebot(in_seq_files, prot_ref_file, min_ident=0.4, min_length=50, alignment_model="glocal", blast_failed_seqs = False, workdir=os.getcwd(), trace=sys.stderr):
	"""
	Runs framebot on the given sequence files
	
	required args = in_seq_files, prot_ref_file
	optional args = min_ident(.4), min_length(50), alignment_model(glocal)
					blast_failed_seqs(False), workdir(cwd), trace(stderr)
	"""
	if not check_unique_files(in_seq_files):
		raise Exception("I won't overwrite output files, two input files have the same name. " + ", ".join(in_seq_files))

	framebot_dir = os.path.join(workdir, "framebot")
	split_dir = os.path.join(framebot_dir, "splits")
	os.mkdir(framebot_dir)
	os.mkdir(split_dir)

	out_seq_files = []

	cmds = []
	blast_cmds = []
	assemble_map = dict()

	if not os.path.exists(prot_ref_file):
		raise Exception("Protein reference file for frame bot, %s, doesn't exist" % prot_ref_file)

	for seq_file in in_seq_files:
		seq_file_name = os.path.split(seq_file.seq_file)[1].split(".")[0]

		seq_file.qual_file = None

		split_in_files = split_seq_file(seq_file.seq_file, 50, split_dir, seq_file_name + ".fasta")

		framebot_files = []
		framebot_failed_files = []
		prot_files = []
		nucl_files = []
		nucl_failed_files = []
		stdout_files = []
		qual_files = []

		for split in split_in_files:
			stem = os.path.join(split_dir, os.path.split(split)[1])
			prot_out_file = stem + "_corr_prot.fasta"
			framebot_file = stem + "_framebot.txt"
			framebot_failed_file = stem + "_failed_framebot.txt"
			nucl_file = stem + "_corr_nucl.fasta"
			qual_file = stem + "_corr_nucl.qual"
			nucl_failed_file = stem + "_failed_nucl.fasta"
			framebot_stdout = stem + "_stdout"

			prot_files.append(prot_out_file)
			framebot_files.append(framebot_file)
			framebot_failed_files.append(framebot_failed_file)
			nucl_files.append(nucl_file)
			nucl_failed_files.append(nucl_failed_file)
			stdout_files.append(framebot_stdout)

			cmd = ["cafe", "-Xmx1g", framebot_class, "--alignment-mode", alignment_model, "--identity-cutoff", min_ident, "--length-cutoff", min_length, "--result-stem", stem]
			if not prot_ref_file.endswith(".idx"):
				cmd.append("--no-metric-search")
			if seq_file.qual_file:
				cmd.extend(["--quality-file", seq_file.qual_file])
				qual_files.append(qual_file)
			
			cmd.extend([prot_ref_file, split])
			cmds.append(Command(cmd, framebot_stdout))

		framebot_out = os.path.join(framebot_dir, seq_file_name + "_framebot.txt")
		framebot_failed_out = os.path.join(framebot_dir, seq_file_name + "_failed_framebot.txt")
		prot_seq_file = os.path.join(framebot_dir, seq_file_name + "_prot_corr.fasta")
		nucl_seq_file = os.path.join(framebot_dir, seq_file_name + "_nucl_corr.fasta")
		if seq_file.qual_file:
			qual_file = os.path.join(framebot_dir, seq_file_name + "_nucl_corr.qual")
		nucl_failed_file = os.path.join(framebot_dir, seq_file_name + "_nucl_failed.fasta")
		stdout_file = os.path.join(framebot_dir, seq_file_name + "_stdout.txt")

		out_seq_files.append(SequenceFile(prot_seq_file, None, seq_file.idmapping, seq_file.sample_mapping))

		assemble_map[prot_seq_file] = prot_files
		assemble_map[framebot_out] = framebot_files
		assemble_map[framebot_failed_out] = framebot_failed_files
		assemble_map[stdout_file] = stdout_files
		assemble_map[nucl_seq_file] = nucl_files
		assemble_map[nucl_failed_file] = nucl_failed_files
		if seq_file.qual_file:
			assemble_map[qual_file] = qual_files
		
		blast_cmds.append(Command([blast, "-p", "blastn", "-d", blastn_db, "-i", nucl_failed_file, "-o", os.path.join(framebot_dir, seq_file_name + "_failed_blastn.txt"), "-m", "8", "-v", "10"]))
		blast_cmds.append(Command([blast, "-p", "blastx", "-d", blastx_db, "-i", nucl_failed_file, "-o", os.path.join(framebot_dir, seq_file_name + "_failed_blastx.txt"), "-m", "8", "-v", "10"]))

	run_commands(cmds, trace, True, split_dir)

	for k in assemble_map.keys():
		cat_files(assemble_map[k], k)
	
	if blast_failed_seqs:
		run_commands(blast_cmds, trace, True, split_dir)

	shutil.rmtree(split_dir)

	return out_seq_files

def chimera_check(in_seq_files, workdir = os.getcwd(), trace = sys.stderr):
	if not check_unique_files(in_seq_files):
		raise Exception("I won't overwrite output files, two input files have the same name. " + ", ".join(in_seq_files))

	
	chimera_dir = os.path.join(workdir, "chimera_check")
	uchime_infile = os.path.join(chimera_dir, "uchime_in.fasta")
	chimera_report_file = os.path.join(chimera_dir, "uchime_report.txt")
	non_chimeras_file = os.path.join(chimera_dir, "non_chimeric.fasta")
	chimera_alignment_file = os.path.join(chimera_dir, "uchime_alignments.txt")

	os.mkdir(chimera_dir)
	if len(in_seq_files) != 1 or not in_seq_files[0].idmapping:
		raise ValueError("Chimera check currently only works with one dereplicated file")

	seq_file = in_seq_files[0]

	"""
	So uchime requires this weird input file format, a fasta file with a 'size annotation', basically a ;size=n;
	but it has to be added to the SEQID not the friggin' label, so I have to completely hack this in to make it
	work with our dereplication tool
	"""
	#first we need to count the abundances
	amplicon_abundances = {}
	for line in open(seq_file.idmapping):
		lexemes = line.strip().split(" ")
		if len(lexemes) != 2:
			continue

		lexemes = lexemes[1].split(",")
		amplicon_abundances[lexemes[0]] = len(lexemes)

	#Now we perform an abomination unto nature by hacking apart the fasta file
	uchime_in = open(uchime_infile, "w")
	for seq in SeqIO.parse(open(seq_file.seq_file), "fasta"):
		uchime_in.write(">{0};size={1};\n{2}\n".format(seq.id, amplicon_abundances[seq.id], str(seq.seq).upper()))
	uchime_in.close();

	cmd = [usearch, "-uchime_denovo", uchime_infile, "-uchimeout", chimera_report_file, "-nonchimeras", non_chimeras_file, "-uchimealns", chimera_alignment_file]

	run_commands([Command(cmd)], trace)
	tmp_file = os.path.join(chimera_dir, "tmp.fasta")
	out = open(tmp_file, "w")
	for seq in SeqIO.parse(open(non_chimeras_file), "fasta"):
		out.write(">{0}\n{1}\n".format(str(seq.id).split(";")[0], seq.seq))
	out.close()
	shutil.move(tmp_file, non_chimeras_file)

	ret = SequenceFile(non_chimeras_file, None, seq_file.idmapping, seq_file.sample_mapping)

	return [ret]

def align(in_seq_files, gene_name, workdir = os.getcwd(), trace = sys.stderr):
	"""
	Aligns the supplied sequence files using a known model (16s, or various functional genes)
	
	required args = in_seq_files, gene_name
	optional args = workdir(cwd), trace(stderr)
	"""
	if not check_unique_files(in_seq_files):
		raise Exception("I won't overwrite output files, two input files have the same name. " + ", ".join(in_seq_files))

	if gene_name == "16s" or gene_name == "archaea16s" or gene_name == "fungi28s":
		align_cmd = cmalign
		model = os.path.join(os.path.join(resources_dir, gene_name), "model.cm")		
	else:
		align_cmd = hmmalign
		model = os.path.join(os.path.join(resources_dir, gene_name), "model.hmm")

	align_dir = os.path.join(workdir, "alignment")
	split_dir = os.path.join(align_dir, "splits")
	os.mkdir(align_dir)
	os.mkdir(split_dir)

	out_seq_files = []
	cmds = []
	assemble_map = dict()

	for seq_file in in_seq_files:
		seq_file_name = os.path.split(seq_file.seq_file)[1]
		
		if "." in seq_file_name:
			seq_file_name = "".join(seq_file_name[:seq_file_name.rfind(".")])
		seq_dir = os.path.join(split_dir, seq_file_name)
		os.mkdir(seq_dir)

		split_in_files = split_seq_file(seq_file.seq_file, 250, split_dir, seq_file_name + ".fasta")

		for split in split_in_files:
			split_out_file = os.path.join(seq_dir, os.path.split(split)[1] + ".out")

			cmd = [align_cmd]
			if gene_name == "16s" or gene_name == "archaea16s" or gene_name == "fungi28s":
				cmd.append("-g").append("--noprob")
			else:
				cmd.append("--allcol")

			cmd.extend(["-o", split_out_file, model, split])
			cmds.append(Command(cmd))

		out_seq_file = os.path.join(align_dir, seq_file_name + "_aligned.fasta")
		out_seq_files.append(SequenceFile(out_seq_file, seq_file.qual_file, seq_file.idmapping, seq_file.sample_mapping))
		assemble_map[out_seq_file] = seq_dir

	run_commands(cmds, trace, True, split_dir)

	merge_cmds = []
	for k in assemble_map.keys():
		cmd = ["cafe", "-Xmx1g", align_merger_class, assemble_map[k], k]
		merge_cmds.append(Command(cmd))

	for out_seq_file in out_seq_files:
		cmd = ["cafe", "-Xmx1g", aligner_stats_class]
		
		if out_seq_file.idmapping:
			cmd.extend(["-i", out_seq_file.idmapping])
		if out_seq_file.sample_mapping:
			cmd.extend(["-s", out_seq_file.sample_mapping])
		
		cmd.extend([out_seq_file.seq_file, out_seq_file.seq_file])
		
		merge_cmds.append(Command(cmd))

	run_commands(merge_cmds, trace)
		

	shutil.rmtree(split_dir)

	return out_seq_files

def dereplicate(in_seq_files, unaligned = True, mask_seq = None, prefix = "all_seqs", workdir = os.getcwd(), trace = sys.stderr):
	"""
	Dereplicates the given files in to a single nonredundant file
	resulting SequenceFile has the idmapping and sample_mapping fields
	populated, optionally masking sequences (mask sequence MUST be supplied
	to dereplicate multiple aligned files)
	
	required args = in_seq_files
	optional args = unaligned(True), mask_seq(None), workdir(cwd), trace(stderr)
	"""

	derep_file = os.path.join(workdir, "%s_derep.fasta" % prefix)
	qual_file = os.path.join(workdir, "%s_derep.qual" % prefix)
	id_file = os.path.join(workdir, "%s.ids" % prefix)
	sample_file = os.path.join(workdir, "%s.samples" % prefix)

	if unaligned:
		if mask_seq:
			sys.stderr.write("Specified mask sequence and unaligned...ignoring mask sequence\n")
		mode = "--unaligned"
	else:
		if mask_seq:
			mode = "--model-only=%s" % mask_seq
		else:
			mode = "--aligned"
	
	cmd = ["cafe", "-Xmx2g", cluster_main_class, "derep", mode, "-o", derep_file, id_file, sample_file]

	qual_files = []
	for seq_file in in_seq_files:
		if seq_file.qual_file and os.path.exists(seq_file.qual_file):
			qual_files.append(seq_file.qual_file)
		cmd.append(seq_file.seq_file)
	
	run_commands([Command(cmd)], trace)
	if len(qual_files) > 0:
		cat_files(qual_files, qual_file, False)
		ret = SequenceFile(derep_file, qual_file, id_file, sample_file)
	else:
		ret = SequenceFile(derep_file, None, id_file, sample_file)

	return [ret]

def distance_matrix(in_seq_files, is_nucl, mask_seq = None, cutoff = 0.15, workdir = os.getcwd(), trace = sys.stderr):
	"""
	Creates a distance matrix for each of the input files
	each input file MUST have an id mapping (also must have
	a sample mapping if you plan to cluster)
	
	required args = in_seq_files, is_nucl
	optional args = mask_seq(None), workdir(cwd), trace(stderr)
	"""
	dist_dir = os.path.join(workdir, "dist_matrix")
	os.mkdir(dist_dir)
	
	if cutoff:
		cutoff = str(cutoff)
		
	ret = []
	for seq_file in in_seq_files:
		if not seq_file.idmapping:
			raise ValueError("Sequence file must have id mapping")
	
		matrix_file = os.path.join(dist_dir, "%s_matrix.bin" % (os.path.split(seq_file.seq_file)[1]))
	
		cmd = ["cafe", "-Xmx2g", cluster_main_class, "dmatrix", "--id-mapping", seq_file.idmapping, "--in", seq_file.seq_file, "--outfile", matrix_file]
		if mask_seq:
			cmd.extend(["--mask", mask_seq])
		if not is_nucl:
			cmd.extend(["-l", "25"])
		if cutoff:
			cmd.extend(["--dist-cutoff", cutoff])
		
		run_commands([Command(cmd)], trace)
	
		ret.append(SequenceFile(matrix_file, None, seq_file.idmapping, seq_file.sample_mapping))
		
	return ret

def cluster(in_files, clust_method = "complete", step = 0.01, workdir = os.getcwd(), trace = sys.stderr):
	"""
	clusters 
	"""
	clust_dir = os.path.join(workdir, "clustering")
	os.mkdir(clust_dir)
	step = str(step)

	ret = []
	for matrix_file in in_files:
		if not matrix_file.idmapping or not matrix_file.sample_mapping:
			raise ValueError("Sequence file must have id and sample mapping")

		clust_file = os.path.join(clust_dir, os.path.split(matrix_file.seq_file)[1].replace("_matrix.bin", "") + "_" + clust_method + ".clust")
		cmd = ["cafe", "-Xmx2g", cluster_main_class, "cluster", "--method", clust_method, "--id-mapping", matrix_file.idmapping, "--sample-mapping", matrix_file.sample_mapping, "--dist-file", matrix_file.seq_file, "--outfile", clust_file, "--step", step]
		
		run_commands([Command(cmd)], trace)

		ret.append(SequenceFile(clust_file, matrix_file.qual_file, matrix_file.idmapping, matrix_file.sample_mapping))

		os.remove(matrix_file.seq_file)	#these things are huge and not very useful to the user...
		
	cmd = ["cafe", "-Xmx2g", cluster_stats_class, clust_dir]
	cmd.extend([x.seq_file for x in ret])
	run_commands([Command(cmd)], trace)

	return ret

def rep_seqs(in_clust_files, aligned_seq_file, cutoff, mask_seq = None, out_dir = "representative_seqs", workdir = os.getcwd(), trace = sys.stderr):
	rep_seq_dir = os.path.join(workdir, out_dir)
	os.mkdir(rep_seq_dir)
	
	if len(in_clust_files) != 1:
		raise Exception("Representative sequences can only be found for a single cluster file...")

	clust_file = in_clust_files[0].seq_file
	cutoff = str(cutoff)
	
	cmd = ["cafe", "-Xmx2g", cluster_main_class, "rep-seqs", "--out", rep_seq_dir]
	if mask_seq:
		cmd.append("--mask-seq=%s" % mask_seq)
	if aligned_seq_file.idmapping:
		cmd.extend(["--id-mapping", aligned_seq_file.idmapping])
	
	cmd.extend([in_clust_files[0].seq_file, cutoff, aligned_seq_file.seq_file])
	run_commands([Command(cmd)], trace)

	ret = []
	for f in glob.glob(os.path.join(rep_seq_dir, "*.fasta")):
		ret.append(SequenceFile(f, None))

	return ret

def refresh_mappings(in_seq_files, out_dir="filtered_mapping", workdir = os.getcwd(), trace = sys.stderr):
	if len(in_seq_files) != 1:
		raise Exception("Can't refresh mapping with multiple sequence files")
	
	seq_file = in_seq_files[0]

	if not seq_file.idmapping or not seq_file.sample_mapping:
		raise Exception("Both id mapping and sample mapping must be present to refresh mappings")

	mapping_dir = os.path.join(workdir, out_dir)
	os.mkdir(mapping_dir)
	
	new_id_mapping = os.path.join(mapping_dir, "filtered_ids.txt")
	new_sample_mapping = os.path.join(mapping_dir, "filtered_samples.txt")
	cmd = ["cafe", "-Xmx2g", cluster_main_class, "refresh-mappings", seq_file.seq_file, seq_file.idmapping, seq_file.sample_mapping, new_id_mapping, new_sample_mapping]
	
	run_commands([Command(cmd)], trace)
	
	seq_file.idmapping = os.path.abspath(new_id_mapping)
	seq_file.sample_mapping = os.path.abspath(new_sample_mapping)

	return in_seq_files
	
def explode_mappings(in_seq_files, out_dir="explode_mappings", workdir = os.getcwd(), trace = sys.stderr):
	if len(in_seq_files) != 1:
		raise Exception("Can't refresh mapping with multiple sequence files")
	
	seq_file = in_seq_files[0]

	if not seq_file.idmapping or not seq_file.sample_mapping:
		raise Exception("Both id mapping and sample mapping must be present to explode mappings")
	
	explode_dir = os.path.join(workdir, out_dir)
	os.mkdir(explode_dir)
	
	cmd = ["cafe", "-Xmx2g", cluster_main_class, "explode-mappings", "-w", "-o", explode_dir, seq_file.idmapping, seq_file.sample_mapping, seq_file.seq_file]
	
	run_commands([Command(cmd)], trace)
	
	ret = []
	for f in os.listdir(explode_dir):
		ret.append(SequenceFile(os.path.join(explode_dir, f), seq_file.qual_file))
	return ret

def jaccard_sorensen(cluster_files, min_cutoff = 0, max_cutoff = 0.15, out_dir="jaccard_and_sorensen", workdir = os.getcwd(), trace = sys.stderr):
	if len(cluster_files) != 1:
		raise Exception("Cannot run jaccard and sorensen with multiple cluster files")
	
	cluster_file = cluster_files[0]
	min_cutoff = str(min_cutoff)
	max_cutoff = str(max_cutoff)
	
	js_dir = os.path.join(workdir, out_dir)
	os.mkdir(js_dir)
	
	cmd = ["cafe", "-Xmx2g", jaccard_sorensen_class, "--jaccard", "--sorensen", "--lower-cutoff", min_cutoff, "--upper-cutoff", max_cutoff, "-r", js_dir]

	if os.path.exists("/usr/bin/R"):
		cmd.extend(["--R-location", "/usr/bin/R"])

	cmd.append(cluster_file.seq_file)
	
	run_commands([Command(cmd)], trace)
	
	return None

def shannon_chao1(cluster_files, out_dir="shannon_chao1", workdir = os.getcwd(), trace = sys.stderr):
	if len(cluster_files) != 1:
		raise Exception("Cannot run shannon&chao1 with multiple cluster files")
	
	cluster_file = cluster_files[0]
	
	shannon_chao_dir = os.path.join(workdir, out_dir)
	out_file = os.path.split(cluster_file.seq_file)[1]
	out_file = out_file[:out_file.rfind(".")] + "_shannon_chao1.txt"
	out_file = os.path.join(shannon_chao_dir, out_file)
	os.mkdir(shannon_chao_dir)
	
	cmd = ["cafe", "-Xmx2g", shannon_chao_class, cluster_file.seq_file, out_file]
	
	run_commands([Command(cmd)], trace)
	
	return None

def rarefaction(cluster_files, out_dir="rarefaction", workdir = os.getcwd(), trace = sys.stderr):
	if len(cluster_files) != 1:
		raise Exception("Cannot run rarefaction with multiple cluster files")
	
	cluster_file = cluster_files[0]
	
	rarefaction_dir = os.path.join(workdir, out_dir)
	os.mkdir(rarefaction_dir)
	
	cmd = ["cafe", "-Xmx2g", rarefaction_class, cluster_file.seq_file, rarefaction_dir, "true"]
	
	run_commands([Command(cmd)], trace)
	
	return None
	
def dealign_seqs(in_seq_files, trace = sys.stderr):
	ret = []
	for seq_file in in_seq_files:
		out_file = seq_file.seq_file
		out_file = os.path.join(os.path.split(out_file)[0], "unaligned_" + os.path.split(out_file)[1])
		
		out = open(out_file, "w")
		for seq in SeqIO.parse(open(seq_file.seq_file), "fasta"):
			if seq.id[0] == "#":
				continue
			out.write(">%s\n%s\n" % (seq.id, str(seq.seq).replace("-", "").replace(".", "").lower()))
		out.close()
		ret.append(SequenceFile(out_file, seq_file.qual_file, seq_file.idmapping, seq_file.sample_mapping))
		
	return ret


