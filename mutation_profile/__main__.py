#!/usr/bin/python

# __main__.py

import sys
import os
import argparse
import textwrap
from Bio import SeqIO, AlignIO
import pandas


from mutation_profile.mutation_profile import get_ref_coords, mut_profile, summary

version = "0.2.1"

def main():
    
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
	group0.add_argument("-f", "--fasta", dest="fasta", required=False, type=str, help="[MANDATORY] Input sequence file (fasta)")
	group0.add_argument("-m", "--mutation_list", dest="mutation", required=False, type=str, help="[MANDATORY] Input mutation list that can be: 1) single-column file with 1-based reference position\
						information (in this case the fasta file must be a multiple sequence alignment of all the sequences of interest); OR 2) tsv file with the columns POS, REF, and ALT \
						where POS = 1-based reference position. If you want to include information for more than one sample per position, add also the column 'ID' (note that the order of the \
						columns is not important but their name is!)")
	group0.add_argument("-r", "--reference", dest="ref", type=str, required=False, help="[MANDATORY] Reference sequence name")
	group0.add_argument("-b", "--before", dest="before", type=int, default=5, help="[OPTIONAL] Number of nucleotides to report BEFORE the mutation (default = 5)")
	group0.add_argument("-a", "--after", dest="after", type=int, default=5, help="[OPTIONAL] Number of nucleotides to report AFTER the mutation (default = 5)")
	group0.add_argument("-p", "--profiles", dest="profiles", type=str, default="GA>AA,TC>TT", help="[OPTIONAL] Comma-separated list of 2bp-mutational profiles of interest (upper-case!). \
						Default = 'GA>AA,TC>TT'")
	group0.add_argument("-o", "--output", dest="output", type=str, default="Mutation_profile", help="[OPTIONAL] Tag for output file name. Default = Mutation_profile")
	group0.add_argument("-v", "--version", dest="version", action="store_true", help="Print version and exit")

	args = parser.parse_args()
	
	# check if version	----------
	
	if args.version:
		print("version:", version)
		sys.exit()

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
	
	
	print("Done!")


if __name__ == "__main__":
	
	main()
