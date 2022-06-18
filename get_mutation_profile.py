#!/usr/bin/env	python3

"""
Given a list of mutations and a fasta sequence obtain a table with the respective mutation profile

By Veronica Mixao
@INSA
"""

import os
import sys
#from shutil import which
import argparse
import textwrap
from Bio import SeqIO, AlignIO
import pandas


# functions	----------

def get_ref_coords(seq):
	""" get a dictionary with d[alignment position] = reference position
	input: alignment
	output: dictionary coords[ref] = alignment
	"""
	
	coords = {}
	i = 0 # ref
	j = 0 # del count
	k = 0 # align
	
	for nucl in seq:
		k += 1
		if nucl == "-":
			if j == 0:
				j = 1 # del count
			else:
				j += 1
			code = str(i) + "." + str(j)
			coords[code] = k
		else:
			i += 1
			j = 0
			coords[i] = k
			
	return coords
	
	
def mut_profile(sequences, before, after, ref, ref_seq, coords, profiles, mutation_df):
	""" obtain a dataframe with a summary of the mutation profile
	input: fasta and mutations
	output: dataframe """
	
	prof = profiles.split(",")
	ref_seq_nogaps = ref_seq.ungap("-")
	
	if "POS" in mutation_df.columns and "ALT" in mutation_df.columns and "REF" in mutation_df.columns:
		mutations = mutation_df["POS"].values.tolist()
		ref_nucl = mutation_df["REF"].values.tolist()
		alt_nucl = mutation_df["ALT"].values.tolist()
		if "ID" in mutation_df.columns:
			samples = mutation_df[mutation_df.columns[0]].values.tolist()
			info = {"sample": samples, "ref_position": mutations, "ref": ref_nucl, "alt": alt_nucl, "motif_ref": [], "observed_profile": [],  "profile_of_interest": []}
		else:
			info = {"ref_position": mutations, "ref": ref_nucl, "alt": alt_nucl, "motif_ref": [], "observed_profile": [],  "profile_of_interest": []}
		
		for index, row in mutation_df.iterrows():
			pos = int(row["POS"])
			pos_0 = pos - 1
			ref_motif = ""
			for p in range(int(pos) - int(before), int(pos) + (int(after) + 1)):
				ref_motif += ref_seq_nogaps[p-1]
			info["motif_ref"].append(ref_motif)
				
			profile1_REF = str(ref_seq_nogaps[pos-2]) + str(ref_seq_nogaps[pos-1])
			profile1_SEQ = str(ref_seq_nogaps[pos-2]) + str(row["ALT"])
			profile1 = profile1_REF + ">" + profile1_SEQ
		
			profile2_REF = str(ref_seq_nogaps[pos-1]) + str(ref_seq_nogaps[pos])
			profile2_SEQ = str(row["ALT"]) + str(ref_seq_nogaps[pos])
			profile2 = profile2_REF + ">" + profile2_SEQ
				
			info["observed_profile"].append(profile1 + " or " + profile2)
				
			profile_of_interest = []
			if profile1 in prof:
				profile_of_interest.append(profile1)
			if profile2 in prof:
				profile_of_interest.append(profile2)
				
			if len(profile_of_interest) == 0:
				info["profile_of_interest"].append("other")
			else:
				info["profile_of_interest"].append(" or ".join(profile_of_interest))
		
		for col in mutation_df.columns:
			if col != "POS" and col != "REF" and col != "ALT" and col != "ID":
				if col not in info.keys():
					info[col] = mutation_df[col].values.tolist()
	else:
		if len(sequences) == 1: # alignment was not provided
			print("Only 1 sequence provided and I cannot find alternative alleles in the mutation list! Cannot continue!!!")
			sys.exit()
		else: 
			positions = mutation_df[mutation_df.columns[0]].values.tolist()
			info = {"sample": [], "ref_position": [], "ref": [], "alt": [], "motif_ref": [], "motif_sample": [], "observed_profile": [],  "profile_of_interest": []}
			
			for pos in positions:
				pos_0 = int(coords[pos]) - 1
				for record in sequences:
					if record.id != ref:
						seq = record.seq.upper()
						info["sample"].append(record.id)
						info["ref_position"].append(pos)
						info["ref"].append(ref_seq[pos_0])
						info["alt"].append(seq[pos_0])
				
						ref_motif = ""
						seq_motif = ""
						for p in range(int(pos) - int(before), int(pos) + (int(after) + 1)):
							ref_motif += ref_seq_nogaps[p-1]
							seq_motif += seq[int(coords[p])-1]
						info["motif_ref"].append(ref_motif)
						info["motif_sample"].append(seq_motif)
				
						profile1_REF = str(ref_seq_nogaps[pos-2]) + str(ref_seq_nogaps[pos-1])
						profile1_SEQ = str(seq[int(coords[pos-1])-1]) + str(seq[int(coords[pos])-1])
						profile1 = profile1_REF + ">" + profile1_SEQ
		
						profile2_REF = str(ref_seq_nogaps[pos-1]) + str(ref_seq_nogaps[pos])
						profile2_SEQ = str(seq[int(coords[pos])-1]) + str(seq[int(coords[pos+1])-1])
						profile2 = profile2_REF + ">" + profile2_SEQ
				
						info["observed_profile"].append(profile1 + " or " + profile2)
				
						profile_of_interest = []
						if profile1 in prof:
							profile_of_interest.append(profile1)
						if profile2 in prof:
							profile_of_interest.append(profile2)
				
						if len(profile_of_interest) == 0:
							info["profile_of_interest"].append("other")
						else:
							info["profile_of_interest"].append(" or ".join(profile_of_interest))
			
	df = pandas.DataFrame(info)
	
	return df


def summary(mx, out):
	""" summarize the profiles of interest detected
	input: pandas matrix
	output: txt file """
	
	if "sample" in mx.columns:
		for sample in set(mx[mx.columns[0]].values.tolist()):
			data = mx[mx[mx.columns[0]].astype(str) == sample] # filter the dataframe
			observations = data["profile_of_interest"].values.tolist()
			counter = {}
			for obs in set(observations):
				counter[obs] = observations.count(obs)
				info2report = []
				if len(counter.keys()) > 0:
					for v in sorted(counter, key=counter.get, reverse=True):
						rel_freq = float(counter[v]/len(observations))
						statistics = str(v) + " (" + str(round(rel_freq * 100,1)) + "%)" 
						info2report.append(statistics)
					joint = ", ".join(info2report) + " (n = " + str(len(observations)) + ")"
				else:
					joint = " - "
			print("\tPatterns of interest found in " + str(sample) + ": " + joint)
			
		data2report = {"ref_position": [], "ref": [], "alt": [], "motif_ref": [], "motif_sample": [], "observed_profile": [], "profile_of_interest": []}
		keys = ["alt", "motif_sample", "observed_profile", "profile_of_interest"]
		checked = []
		for mut in mx["ref_position"].values.tolist():
			if mut not in checked:
				data = mx[mx["ref_position"].astype(str) == str(mut)]
				data2report["ref_position"].append(str(mut))
				data2report["ref"].append(data["ref"].values.tolist()[0])
				data2report["motif_ref"].append(data["motif_ref"].values.tolist()[0])
				
				for parameter in keys:
					if parameter in data.columns:
						observations = data[parameter].values.tolist()
						counter = {}
						for obs in set(observations):
							counter[obs] = observations.count(obs)
							info2report = []
							if len(counter.keys()) > 0:
								for v in sorted(counter, key=counter.get, reverse=True):
									rel_freq = float(counter[v]/len(observations))
									statistics = str(v) + " (" + str(round(rel_freq * 100,1)) + "%)" 
									info2report.append(statistics)
								joint = ", ".join(info2report) + " (n = " + str(len(observations)) + ")"
							else:
								joint = " - "
						data2report[parameter].append(joint)
				
					else:
						data2report[parameter].append("-")
				checked.append(mut)	
		
		data2report_df = pandas.DataFrame(data2report)
		data2report_df.to_csv(out + "_report.tsv", index = False, header=True, sep ="\t")
					
	else:
		data = mx # filter the dataframe
		observations = data["profile_of_interest"].values.tolist()
		counter = {}
		for obs in set(observations):
			counter[obs] = observations.count(obs)
			info2report = []
			if len(counter.keys()) > 0:
				for v in sorted(counter, key=counter.get, reverse=True):
					rel_freq = float(counter[v]/len(observations))
					statistics = str(v) + " (" + str(round(rel_freq * 100,1)) + "%)" 
					info2report.append(statistics)
				joint = ", ".join(info2report) + " (n = " + str(len(observations)) + ")"
			else:
				joint = " - "
		print("\tPatterns of interest found: " + joint)
		
		
# main	----------

if __name__ == "__main__":
    
	# argument options
    
	parser = argparse.ArgumentParser(prog="partitioning_HC.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                          get_mutation_profile.py                            #
									#                                                                             #
									###############################################################################  
									                            
									This script was developed to rapidly obtain the sequence context flanking
									SNPs of interest and determine their 2bp-mutational profile/signature (e.g.
									APOBEC3-mediated viral genome editing).
									
									This script can be run by providing different combinations of inputs
									
									OPTION 1
									
									Input 1: TSV file with the columns POS REF ALT (i.e. 1-indexed reference 
									position, reference allele and alternative allele)
									Input 2: Fasta file including the reference genome
									
									Output 1: TSV file with the mutation context and profile
									
									
									OPTION 2
									
									Input 1: TSV file with the columns ID POS REF ALT (i.e. sample ID, 1-indexed 
									reference position, reference allele and alternative allele)
									Input 2: Fasta file including the reference genome
									
									Output 1: TSV file with the mutation context and profile for each sample present 
									in the TSV input
									Output 2: TSV file with a summary report for each mutation including the
									different patterns observed and their respective frequency
									
									OPTION 3
									
									Input 1: Single-column file with a list of 1-indexed reference positions of 
									interest
									Input 2: Multiple Sequence Alignment (fasta) including the reference genome
									
									Output 1: TSV file with the mutation context and profile for each sample present 
									in the alignment
									Output 2: TSV file with a summary report for each mutation including the
									different patterns observed and their respective frequency
									
									
									NOTE: IN OPTIONS 1 AND 32, THE ORDER OF THE COLUMNS IN THE INPUT 1 IS NOT
									IMPORTANT, BUT THEIR NAME IS (ID, POS, REF, ALT)!!!
									
									-----------------------------------------------------------------------------"""))
	
	group0 = parser.add_argument_group("Mutation profile", "Provide input/output specifications")
	group0.add_argument("-f", "--fasta", dest="fasta", required=True, type=str, help="[MANDATORY] Input sequence file (fasta)")
	group0.add_argument("-m", "--mutation_list", dest="mutation", required=True, type=str, help="[MANDATORY] Input mutation list that can be: 1) single-column file with 1-based reference position\
						information (in this case the fasta file must be a multiple sequence alignment of all the sequences of interest); OR 2) tsv file with the columns POS, REF, and ALT \
						where POS = 1-based reference position. If you want to include information for more than one sample per position, add also the column 'ID' (note that the order of the \
						columns is not important but their name is!)")
	group0.add_argument("-r", "--reference", dest="ref", type=str, required=True, help="[MANDATORY] Reference sequence name")
	group0.add_argument("-b", "--before", dest="before", type=int, default=5, help="[OPTIONAL] Number of nucleotides to report BEFORE the mutation (default = 5)")
	group0.add_argument("-a", "--after", dest="after", type=int, default=5, help="[OPTIONAL] Number of nucleotides to report AFTER the mutation (default = 5)")
	group0.add_argument("-p", "--profiles", dest="profiles", type=str, default="GA>AA,TC>TT", help="[OPTIONAL] Comma-separated list of 2bp-mutational profiles of interest (upper-case!). \
						Default = 'GA>AA,TC>TT'")
	group0.add_argument("-o", "--output", dest="output", type=str, default="Mutation_profile", help="[OPTIONAL] Tag for output file name. Default = Mutation_profile")

	args = parser.parse_args()
	
	# read fasta file
	
	print("Loading the fasta sequence...")
	sequences = AlignIO.read(args.fasta, "fasta")
	print("\tLoaded " + str(len(sequences)) + " sequences.")
	
	# get reference sequence
	
	print("Defining reference sequence...")
	reference = args.ref
	ref_seq = ""
	for record in sequences:
		if record.id == reference:
			ref_seq = record.seq.upper()
	
	if ref_seq == "":
		print("Could not find the reference name in the fasta provided!!! Cannot continue!")
		sys.exit()
	
	# read the mutation list
	
	print("Reading the mutation list...")
	mutation_df = pandas.read_table(args.mutation)
	
	# get coordinate correspondence
	
	print("Get reference and alignment position correspondence...")
	coords = get_ref_coords(ref_seq)
			
	# get profile information
	
	print("Get mutation profile...")
	mx = mut_profile(sequences, int(args.before), int(args.after), reference, ref_seq, coords, args.profiles, mutation_df)
	mx.to_csv(args.output + ".tsv", index = False, header=True, sep ="\t")

	# check percentage of profiles of interest
	
	print("Get summary of the detected profiles of interest.")
	summary(mx, args.output)
	
	"""
	# generate snipit plot
	
	if len(sequences) > 1:
		if which("snipit") is not None:
			print("Running snipit (https://github.com/aineniamh/snipit)...")
			os.system("snipit " + args.fasta + " -r " + args.ref + " -o " + args.output)
	"""
	
	print("Done!")
