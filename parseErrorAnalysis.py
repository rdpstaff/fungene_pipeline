#! /usr/bin/python

import os
import sys
from optparse import OptionParser
from Bio import SeqIO
import re
sys.path.append("/work/wangqion/python_scripts/")
import seq_trimmer_model

base_regex = re.compile("[^A-Za-z]")


class SequenceMatch:
	def __init__(self, query, match):
		self.query = query
		self.match = match
		self.mismatches = list()
		self.indels = list()
		self.avg_qscore = -1
		
	def process_qual(self, line):
		self.avg_qscore = float(line.strip().split()[1])
		
	def process_mismatch(self, line):
	##ref_pos from mismatch and indel starts from 1
		lexemes = line.strip().split()
		mismatch = dict()
		mismatch["r_char"] = lexemes[3]
		mismatch["q_char"] = lexemes[4]
		mismatch["query_pos"] = lexemes[5]
		mismatch["ref_pos"] = lexemes[6]
		mismatch["qscore"] = lexemes[7]
		
		self.mismatches.append(mismatch)
		
	def process_indel(self, line):
	##ref_pos from mismatch and indel starts from 1
		lexemes = line.strip().split()
		indel = dict()
		indel["ref_homo_count"] = lexemes[3]
		indel["query_homo_count"] = lexemes[4]
		indel["badchar"] = lexemes[5]
		indel["query_pos"] = lexemes[6]
		indel["ref_pos"] = lexemes[7]
		indel["qscore"] = lexemes[8]
		self.indels.append(indel)
	
def read_files(match_file, mismatch_file, indel_file, qscore_file):
	seq_dict = dict()
	for line in open(match_file):
		line = line.strip()
		if ( line != "" and line[0] == ">"):
			lexemes = line.split()
			seq_dict[lexemes[1]] = SequenceMatch(lexemes[1], lexemes[2])
			
	#
	for line in open(mismatch_file):
		line = line.strip()
		if ( line != ""):
			lexemes = line.split()
			seq_dict[lexemes[0]].process_mismatch(line) 
	#		
	for line in open(indel_file):
		line = line.strip()
		if ( line != ""):
			lexemes = line.split()
			seq_dict[lexemes[0]].process_indel(line) 
			
	#	
	if qscore_file != None:
		for line in open(qscore_file):
			line = line.strip()
			if ( line != ""):
				lexemes = line.split()
				seq_dict[lexemes[0]].process_qual(line) 
			
	return seq_dict		

def removeBadseq(bad_idfile, seq_dict):
	
	for line in open(bad_idfile):
		line = line.strip()
		if ( line != ""):
			lexemes = line.split()
			if ( lexemes[0] in seq_dict):
				del seq_dict[lexemes[0]]

def removeFailedChopseq(chop_seq_dict, seq_dict):
	for seqID in seq_dict.keys():
		if ( seqID not in chop_seq_dict):
			del seq_dict[seqID]		
				
				
def process_std_seq(infile):
## need the definition for the sequences
	std_dict = dict()
	for line in open(infile):
		line = line.strip()
		if ( line != "" and line[0] == ">"):
			lexemes = line.split()			
			if (len(lexemes) <2 ):
				definition = lexemes[0].replace(">", "")
			else:
				definition = ""
				for index in range( 1, len(lexemes)):
					definition +=  " " + lexemes[index]

			std_dict[lexemes[0].replace(">", "")] = definition 
	return std_dict		

		
def process_model_pos_map(infile):
	## unalign_pos (start from 0, need to change to 1 since the ref_pos from mismatch and indel starts from 1), model_pos
	## format: seqid frame unalign_pos model_pos
	
	seq_modelpos_dict = dict()
	for line in open(infile):
		line = line.strip()
		if ( line != ""):
			lexemes = line.split("\t")
			modelpos_dict = seq_modelpos_dict.get(lexemes[0], dict())
			adjusted_unalign_pos = int(lexemes[1]) -1
			modelpos_dict[adjusted_unalign_pos] = lexemes[2]
			seq_modelpos_dict[lexemes[0]] = modelpos_dict
				
	return seq_modelpos_dict				
			
# calculate the total number of seqs matched to the std seqs
def get_totalseq_count(seq_dict, qscore_cutoff):
	total_seq = 0
	for seqID in seq_dict.keys():
		seq = seq_dict[seqID]
		if seq.avg_qscore > qscore_cutoff or seq.avg_qscore == -1:
			total_seq += 1
	return total_seq


	
# print the match count to the std seqs		
def get_match_count(seq_dict, qscore_cutoff):
	match_count_dict = dict()
	for seqID in seq_dict.keys():
		seq = seq_dict[seqID]
		## if there is no quality score info, we just count it
		if seq.avg_qscore > qscore_cutoff or seq.avg_qscore == -1:
			match_count_dict[seq.match] = match_count_dict.get(seq.match,0) +1
	
	return match_count_dict


## print the copies of  the standard sequences from the same organism
def get_std_copy(std_dict):
	name_dict = dict()
	print "\n### standard sequences ###"
	print "standard_id\tdefinition"
	
	for id in sorted(std_dict.keys()):
		print "%s\t%s" %(id, std_dict[id])
		name = std_dict[id]
		name_dict[name] = name_dict.get(name, 0) +1
	
	print "\n### standard sequence copies ###"
	print "standard_definition\tcount"	
	
	for name in sorted(name_dict.keys()):
		print "%s\t%s" %(name, name_dict.get(name))
		


def get_total_match_count(seq_dict, std_dict):
	q0_match_count_dict = get_match_count(seq_dict, 0)	
	totalseq_count = get_totalseq_count(seq_dict, 0)
	print "\n### best reference match count for seq with average QScore >=0 ###"
	print "standard_seqID\tdefinition\tQScore_0_count\tQScore_0_pct"
	for id in sorted(q0_match_count_dict.keys()):
		print "%s\t%s\t%s\t%s" %(id, std_dict[id], q0_match_count_dict[id],  float(100*q0_match_count_dict.get(id,0))/float(totalseq_count))

##compare the match count based on different qscore cutoff
def compare_match_count(seq_dict, std_dict, qscore_cutoff):
	q0_match_count_dict = get_match_count(seq_dict, 0)	
	totalseq_count = get_totalseq_count(seq_dict, 0)
	#qcutoff_match_count_dict = get_match_count(seq_dict, qscore_cutoff)	
	print "\n### percent sequences passed the Qscore cutoff comparing to the ones passing QScore 0 ###" 
	header = "Qscore_cutoff";
	for id in sorted(q0_match_count_dict.keys()):
		header += "\t" + id 
	print "%s" %(header)
	
	for cutoff in range ( 15 , qscore_cutoff):
		qcutoff_match_count_dict = get_match_count(seq_dict, cutoff)	
		outstring = str(cutoff)
		for id in sorted(q0_match_count_dict.keys()):
			outstring += "\t" + str(float(100*qcutoff_match_count_dict.get(id,0))/float(q0_match_count_dict[id]))
		print "%s" %(outstring)
		

# print the number_errors, seq count		
def get_error_count(seq_dict):
	count_dict = dict()
	total_mismatches = 0
	total_indels = 0
	
	for seqID in seq_dict.keys():
		seq = seq_dict[seqID]
		mismatch_indels = len(seq.mismatches) + len(seq.indels)
		count_dict[mismatch_indels] = count_dict.get(mismatch_indels,0) +1
		total_mismatches += len(seq.mismatches)
		total_indels += len(seq.indels)
	
	total_seqs = get_totalseq_count(seq_dict, 0)
	print "\n### total mismatches and indels ###"
	print "Total seqs\tTotal Mismatches\tTotal Indels"
	print "%s\t%s\t%s\t" %(total_seqs, total_mismatches, total_indels)
	print "\n### mismatches and indels sequence count###"
	print "no_mismatch_indels\tcount\tpercent_seq\tpercent_error"
	for key in sorted(count_dict.keys()):
		print "%s\t%s\t%s\t%s" %(key, count_dict[key], float(100*count_dict[key])/float( len(seq_dict.keys()) ), float(100*key*count_dict[key])/float(total_mismatches+total_indels ))	
		
# group by each standard sequence, print the number_errors, seq count		
def get_error_count_by_std(seq_dict):
	match_dict = dict()
	total_mismatch_indels = 0
	seq_error_dict = dict()
	for seqID in seq_dict.keys():
		seq = seq_dict[seqID]
		mismatch_indels = len(seq.mismatches) + len(seq.indels)
		match_dict[seq.match]= match_dict.get(seq.match, dict())
		count_dict = match_dict.get(seq.match)
		count_dict[mismatch_indels] = count_dict.get(mismatch_indels,0) +1
		total_mismatch_indels += mismatch_indels
		seq_error_dict[seq.match] = int(seq_error_dict.get(seq.match, 0)) + mismatch_indels
	
	print "\n### mismatches and indels sequence count group by standard sequence###"
	print "standard_seqID\ttotal_no_mismatch_indels\tpercent_error"
	for std_ID in seq_error_dict.keys():
		print "%s\t%s\t%s" %(std_ID, seq_error_dict.get(std_ID), float(100*seq_error_dict.get(std_ID))/float(total_mismatch_indels))
	print 
	print "standard_seqID\tno_mismatch_indels\tcount\tpercent_error"
	for std_ID in sorted(match_dict.keys()):
		print ""
		count_dict = match_dict.get(std_ID)
		for key in sorted(count_dict.keys()):
			if ( int(seq_error_dict.get(std_ID)) == 0):
				print "%s\t%s\t%s\t%s" %(std_ID, key, count_dict[key], 0)
			else:
				print "%s\t%s\t%s\t%s" %(std_ID, key, count_dict[key], float(100*key*count_dict[key])/float( int(seq_error_dict.get(std_ID))))
		
		
def get_hotspot(seq_dict, seq_modelpos_dict):
	mismatch_count_dict = dict()
	indel_count_dict = dict()
	total_mismatches = 0
	total_indels = 0
	for seqID in seq_dict.keys():
		seq = seq_dict[seqID]
		for mismatch in seq.mismatches:
			if ( int(mismatch["ref_pos"]) in seq_modelpos_dict.get(seq.match)):				
				modelpos = int(seq_modelpos_dict.get(seq.match).get( int(mismatch["ref_pos"])))
			else:
				modelpos = -1
			mismatch_count_dict[modelpos] = mismatch_count_dict.get(modelpos, 0) +1
			total_mismatches += 1
		
		for indel in seq.indels:
			if ( int(indel["ref_pos"]) in seq_modelpos_dict.get(seq.match)) :
				modelpos = int(seq_modelpos_dict.get(seq.match).get(int(indel["ref_pos"])) )
			else:
				modelpos = -1
			indel_count_dict[modelpos] = indel_count_dict.get(modelpos, 0) +1
			total_indels += 1
			
	print "\n### mismatch hot spots###"
	print "std_model_pos\tcount\tpercent mismatches\tcumulative mismatches"
	total = 0
	for key in sorted(mismatch_count_dict.keys()):
		total += mismatch_count_dict[key]
		print "%s\t%s\t%s\t%s" %(key, mismatch_count_dict[key], float(100*mismatch_count_dict[key])/float(total_mismatches ), float(100*total)/float( total_mismatches))	
	
	print "\n### indel hot spots###"
	print "std_model_pos\tcount\tpercent indels\tcumulative indels"
	total = 0
	for key in sorted(indel_count_dict.keys()):
		total += indel_count_dict[key]
		print "%s\t%s\t%s\t%s" %(key, indel_count_dict[key], float(100*indel_count_dict[key])/float( total_indels ), float(100*total)/float( total_indels) )	
		
## remove the mismatch and indels occurs outside the allowed model_pos range
def remove_dontcare_error(seq_dict, seq_modelpos_dict, start_pos, end_pos):
	for seqID in seq_dict.keys():
		seq = seq_dict[seqID]
		temp_mismatches = list()

		for index in range( len(seq.mismatches) ):
			mismatch = seq.mismatches[index]
			## there are cases where the nucleotides don't code a amino acid and don't have corresponding alignment position
			modelpos = int(seq_modelpos_dict.get(seq.match).get( int(mismatch["ref_pos"]), -1))
			if ( modelpos >= start_pos and modelpos <= end_pos):	
				temp_mismatches.append(mismatch)
			#else:
			#	print "remove mismatch  %s\t%s\t%s" %(seqID, modelpos, mismatch)
		seq.mismatches = temp_mismatches
		
		temp_indels = list();
		for index in range( len(seq.indels) ):
			indel = seq.indels[index]
			modelpos = int(seq_modelpos_dict.get(seq.match).get( int(indel["ref_pos"]), -1))
			if ( modelpos >= start_pos and modelpos <= end_pos):
				temp_indels.append(indel)
			#else :
			#	print "remove indel %s\t%s\t%s" %(seqID, modelpos, indel)
		seq.indels = temp_indels

#Mismatch mapped to each standard sequence
## refseq --> ref_pos --> r_char + q_char --> count
def get_mismatch_stdseq(seq_dict, seq_modelpos_dict):
	mismatch_dict = dict()
	totalseq_count_dict = get_match_count(seq_dict, 0)	
	for seqID in seq_dict.keys():
		seq = seq_dict[seqID]
		for mismatch in seq.mismatches:
			mismatch_dict[seq.match]= mismatch_dict.get(seq.match, dict())	
			refpos_dict = mismatch_dict.get(seq.match)	
			refpos_dict[mismatch["ref_pos"]] = refpos_dict.get(mismatch["ref_pos"], dict())
			mismatch_char_dict = refpos_dict.get(mismatch["ref_pos"])	
			concat_key = mismatch["r_char"] + "\t" + mismatch["q_char"]
			mismatch_char_dict[concat_key] = mismatch_char_dict.get(concat_key, 0) +1
					

	print "\n### mismatch map to standard sequence ###"
	print "standard_seqID\tstd_unalign_pos\tstd_model_pos\tr_char\tq_char\tcount\tpercent"
	for std_ID in sorted(mismatch_dict.keys()):
		print ""
		refpos_dict = mismatch_dict.get(std_ID)	
		for ref_pos in sorted(refpos_dict.keys()):
			modelpos = seq_modelpos_dict.get(std_ID).get( int(ref_pos))
			mismatch_char_dict = refpos_dict.get(ref_pos)
			for mismatch_chars in sorted(mismatch_char_dict.keys()):	
				print "%s\t%s\t%s" %(std_ID, mismatch_char_dict.get(mismatch_chars), totalseq_count_dict.get(std_ID))
				print "%s\t%s\t%s\t%s\t%s\t%s" %(std_ID, ref_pos, modelpos, mismatch_chars, mismatch_char_dict.get(mismatch_chars), float(100*mismatch_char_dict.get(mismatch_chars))/float(totalseq_count_dict.get(std_ID)))	
	
 	
#indel mapped to each standard sequence
## refseq --> ref_pos --> count
def get_indel_stdseq(seq_dict, seq_modelpos_dict):
	indel_dict = dict()
	totalseq_count_dict = get_match_count(seq_dict, 0)
	for seqID in seq_dict.keys():
		seq = seq_dict[seqID]
		for indel in seq.indels:
			indel_dict[seq.match] = indel_dict.get(seq.match, dict())
			refpos_dict = indel_dict.get(seq.match)
			refpos_dict[indel["ref_pos"]] = refpos_dict.get(indel["ref_pos"], 0) +1
		
	print "\n### indels map to standard sequence ###"
	print "standard_seqID\tstd_unalign_pos\tstd_model_pos\tcount\tpercent"
	for std_ID in sorted(indel_dict.keys()):
		print ""
		refpos_dict = indel_dict.get(std_ID)	
		for ref_pos in sorted(refpos_dict.keys()):
			modelpos = seq_modelpos_dict.get(std_ID).get( int(ref_pos))
			print "%s\t%s\t%s\t%s\t%s" %(std_ID, ref_pos, modelpos, refpos_dict.get(ref_pos), (float(100*refpos_dict.get(ref_pos))/float(totalseq_count_dict.get(std_ID))) )	

## base substitutions errors
def get_base_sub_error(seq_dict):
	count_dict = dict()
	for seqID in seq_dict.keys():
		seq = seq_dict[seqID]
		for mismatch in seq.mismatches:
			concat_key = mismatch["r_char"] + "\t" + mismatch["q_char"]
			count_dict[concat_key] = count_dict.get(concat_key, 0) +1
	print "\n### base substitutions errors ###"
	print "standard_base\tquery_base\tcount"
	for mismatch_chars in sorted(count_dict.keys()):
		print "%s\t%s" %( mismatch_chars, count_dict.get(mismatch_chars))	
			
## Insertion errors ##
def get_indel_error(seq_dict):
	insertion_count_dict = dict()
	deletion_count_dict = dict()
	for seqID in seq_dict.keys():
		seq = seq_dict[seqID]
		for indel in seq.indels:
			concat_key = indel["ref_homo_count"] + "\t" + indel["query_homo_count"]
			if (int(indel["ref_homo_count"]) < int(indel["query_homo_count"])):
				insertion_count_dict[concat_key] = insertion_count_dict.get(concat_key, 0) +1
			else:
				deletion_count_dict[concat_key] = deletion_count_dict.get(concat_key, 0) +1
	print "\n### Insertion errors ###"
	print "ref_homo_count\tquery_homo_count\tcount"
	for key in sorted(insertion_count_dict.keys()):
		print "%s\t%s" %( key, insertion_count_dict.get(key))
		
	print "\n### deletion errors ###"
	print "ref_homo_count\tquery_homo_count\tcount"
	for key in sorted(deletion_count_dict.keys()):
		print "%s\t%s" %( key, deletion_count_dict.get(key))	

## of sequences with # of mismatch + indels binned by Qscore
def get_error_by_qscore(seq_dict):
	qscore_dict = dict()
	for seqID in seq_dict.keys():
		seq = seq_dict[seqID]
		qscore = int(seq.avg_qscore)
		qscore_dict[qscore] = qscore_dict.get(qscore, dict())
		count_dict = qscore_dict.get(qscore)
		count_dict["no_seqs"] = count_dict.get("no_seqs", 0) +1
		count_dict["no_mis_indels"] = count_dict.get("no_mis_indels", 0) + len(seq.mismatches) + len(seq.indels)
		
	print "\n### number of sequences with # of mismatch + indels binned by Qscore ###"
 	print "qscore\tno_seqs\tno_mis_indels\terror_per_seq"
 	for qscore in sorted(qscore_dict.keys()):
 		count_dict = qscore_dict.get(qscore)
 		print "%s\t%s\t%s\t%s" %(qscore, count_dict["no_seqs"], count_dict["no_mis_indels"], float(count_dict["no_mis_indels"])/count_dict["no_seqs"] )


#avg qscore histogram
def qscore_histogram(seq_dict):
	count_dict = dict()
	sum = 0
	num_seqs = 0
	
	for seqID in seq_dict.keys():
		seq = seq_dict[seqID]
		qscore = int(seq.avg_qscore)
		count_dict[qscore] = count_dict.get(qscore,0) +1
		sum += seq.avg_qscore
		num_seqs += 1
	
	print "\n### Q score ###"
	if ( num_seqs== 0):
		print "Average Q score: NA" 
	else:
		print "Average Q score: %s" %( (sum/num_seqs))
	print "\n### Q score histogram ###"
	print "Qscore\tcount\tpercent"
	total = 0
	for key in sorted(count_dict.keys()):
		total += count_dict[key]
		print "%s\t%s\t%s" %(key, count_dict[key], float(100*count_dict[key])/float( len(seq_dict.keys()) ))	

## plot number of errors and number of seqs passed the Q score filter
def qscore_seqpassed(seq_dict):
	qscore_dict = dict();
	min_display_qscore = 10;
	max_qscore = 40;
	max_error = 10;
	
	for seqID in seq_dict.keys():
		seq = seq_dict[seqID]
		mismatch_indels = len(seq.mismatches) + len(seq.indels)
		q = 0
		while ( q <= max_qscore):
			if ( seq.avg_qscore >= q):
				qscore_dict[q] = qscore_dict.get(q, dict())
				error_dict = qscore_dict.get(q)
				if ( mismatch_indels < max_error):
					error_dict[mismatch_indels] = error_dict.get(mismatch_indels, 0) +1
				else :
					error_dict[max_error] = error_dict.get(max_error, 0) +1
			q += 1
			
	print "\n## percent of seqs with the specified error that passed the Q score"
	header = "Qscore"
	for e in range(0, max_error):
		header = header + "\tE" + str(e)
	header += "\tE>=" + str(max_error)
	print "%s" %(header)
	
	q0_error_dict = qscore_dict.get(0)
	q = min_display_qscore
	while ( q <= max_qscore):
		if q not in qscore_dict.keys():
			q +=1
			continue;
		error_dict = qscore_dict.get(q)
		val = str(q)
		for e in range(0, max_error+1):
			if e in error_dict.keys():
				val += "\t" + str( float(qscore_dict.get(q).get(e)) *100/ float(q0_error_dict.get(e)) ) 
			else: 
				val += "\t" + str(0)
		print "%s" %(val)
		q +=1

'''
	print "\n## number of seqs with the specified error that passed the Q score"
	print "%s" %(header)
	q0_error_dict = qscore_dict.get(0)
	q = min_display_qscore
	while ( q <= max_qscore):
		if q not in qscore_dict.keys():
			q +=1
			continue;
		error_dict = qscore_dict.get(q)
		val = str(q)
		for e in range(0, max_error+1):
			if e in error_dict.keys():
				val += "\t" + str( float(qscore_dict.get(q).get(e))  ) 
			else:
				val += "\t" + str(0)
		print "%s" %(val)
		q +=1
'''

## calculate error for seqs after chopping based on the alignment
def	calSeqError(seq_dict, chop_seq_dict):
	count_dict = dict()
	total_number_seqs = 0
	
	for seqID in chop_seq_dict.keys():
		if (seqID not in seq_dict.keys() ):
			continue
		if (seqID.startswith("#=") ):	
			continue
		
		total_number_seqs += 1

		seq = seq_dict[seqID]
		chopped_len = len(re.subn(base_regex, "", str(chop_seq_dict.get(seqID)))[0])
		mismatch_indels = len(seq.mismatches) + len(seq.indels)
		# step 0.1%
		pct_error = 1000*float(mismatch_indels) /chopped_len
		#print "found %s\t%s\t%s\t%s\t%s" %(seqID, pct_error, mismatch_indels, chopped_len, int(pct_error) )
		
		count_dict[int(pct_error)] = count_dict.get(int(pct_error),0) +1
		
	print "\n### error/seq after chopping the seqs ###"
	print "percent_error\tpercent_seq_passed"
	cum = 0
	for key in sorted(count_dict.keys()):
		cum += count_dict[key]
		## back to percent
		print "%s\t%s" %(float(key)/float(10), float(cum)/float(total_number_seqs))	

## remove seq with certain % error
def rmSeqwithError(seq_dict, chop_seq_dict, error_cutoff):
	for seqID in chop_seq_dict.keys():
		if (seqID not in seq_dict.keys() ):
			continue
		seq = seq_dict[seqID]
		chopped_len = len(re.subn(base_regex, "", str(chop_seq_dict.get(seqID)))[0])
		mismatch_indels = len(seq.mismatches) + len(seq.indels)
		error = float(mismatch_indels) /chopped_len
		
		if ( error > error_cutoff):
			del seq_dict[seqID]
	
def getFileName(fileNameWithPath):
	lexemes = fileNameWithPath.split('/')
	return lexemes[len(lexemes)-1]
	
if __name__ == "__main__":
	usage="usage: %prog [options] pairwise.txt mismatch.txt indel.txt standard_nucl_seqs.fa"
	
	parser = OptionParser(usage=usage)
	parser.add_option("-q", "--qual", dest="quality_file",help="quality output file from errorcheck program", metavar="FILE")
	parser.add_option("-i", "--ignore", dest="ignore_ids",help="ignored sequence ids", metavar="FILE")
	parser.add_option("-c", "--ignore_chimera", dest="ignore_chimera_ids",help="ignored chimera sequence ids", metavar="FILE")
	parser.add_option("-m", "--model_pos_mapping", dest="model_pos_mapping",help="model position mapping file for standard seqs", metavar="FILE")
	parser.add_option("-s", "--start", dest="start_model_pos", help="start DNA model position to be included")
	parser.add_option("-e", "--end", dest="end_model_pos", help="end DNA model position to be included")
	parser.add_option("-a", "--align", dest="aligned_query_file", help="to calculate errors in the query file")

	(options, args) = parser.parse_args()
	if len(args) != 4:
		parser.error("Incorrect number of arguments")
	argsStr = ''
	for arg in args:
		argsStr += getFileName(arg) + " "
	print "arguments: %s" %( argsStr)

	qual_file = None
	if options.quality_file:
		qual_file = options.quality_file
	if options.start_model_pos:
		if not options.end_model_pos:
			options.end_model_pos = sys.maxint 
	if options.end_model_pos:
		if not options.start_model_pos:
			options.start_model_pos = -sys.maxint -1

	seq_dict = read_files(args[0], args[1], args[2], qual_file);
	std_dict = process_std_seq(args[3])
	if options.model_pos_mapping:
		print "read model position mapppig file: %s" %( getFileName(options.model_pos_mapping))
		seqs_model_pos_dict = process_model_pos_map(options.model_pos_mapping)
	if options.ignore_ids:
		print "remove ignored ids from file: %s" %(getFileName(options.ignore_ids))
		removeBadseq(options.ignore_ids, seq_dict)
	if options.ignore_chimera_ids:
		print "remove ignored chimera ids from file: %s" %(getFileName(options.ignore_chimera_ids))
		
		removeBadseq(options.ignore_chimera_ids, seq_dict)
	if options.start_model_pos or options.end_model_pos:
		print "remove positons outside of this range: \t%s\t%s" %(options.start_model_pos, options.end_model_pos)
		if not options.model_pos_mapping :
			parser.error("start_model_pos or end_model_pos requires a model_pos_mapping file")
		remove_dontcare_error(seq_dict, seqs_model_pos_dict, int(options.start_model_pos), int(options.end_model_pos))
        	
	if options.aligned_query_file:
		if not options.start_model_pos and not options.end_model_pos :
			parser.error("requires start_model_pos and end_model_pos to calculate the error in the sequences using the alignment file to slice")
		else:	
			print "chop the sequences form start to end: \t%s\t%s" %(options.start_model_pos, options.end_model_pos)	
			seq_before_chop = get_totalseq_count(seq_dict, 0)
			chop_seq_dict= seq_trimmer_model.chop(options.aligned_query_file, int(options.start_model_pos), int(options.end_model_pos))
			removeFailedChopseq(chop_seq_dict, seq_dict)
			
			seq_after_chop = get_totalseq_count(seq_dict, 0)
			
			calSeqError(seq_dict, chop_seq_dict)
			rmSeqwithError(seq_dict, chop_seq_dict, 0.03)
			
			seq_after_rmerror = get_totalseq_count(seq_dict, 0)

			print "\ntotal number of seqs before chop\tafter chop\tafter remove 3%error"
			print "%s\t%s\t%s" %(seq_before_chop, seq_after_chop, seq_after_rmerror)
			q0_match_count_dict = get_match_count(seq_dict, 0)	
			totalseq_count = get_totalseq_count(seq_dict, 0)
			print "\n### best reference match count for seq " 
			print "standard_seqID\tdefinition\tseq_count\tseq%" 
			for id in sorted(q0_match_count_dict.keys()):
				print "%s\t%s\t%s\t%s" %(id, std_dict[id], q0_match_count_dict[id], float(100*q0_match_count_dict.get(id,0))/float(totalseq_count))

			sys.exit() 
		
	print ""
	get_std_copy(std_dict)
	get_error_count(seq_dict)
	get_total_match_count(seq_dict, std_dict)
	compare_match_count(seq_dict, std_dict,40)
	qscore_seqpassed(seq_dict)
	get_error_count_by_std(seq_dict)

	qscore_histogram(seq_dict)
	get_error_by_qscore(seq_dict)
	get_base_sub_error(seq_dict)
	get_indel_error(seq_dict)
	if options.model_pos_mapping :
		get_hotspot(seq_dict,seqs_model_pos_dict) 
		get_mismatch_stdseq(seq_dict, seqs_model_pos_dict)
		get_indel_stdseq(seq_dict, seqs_model_pos_dict)
		
