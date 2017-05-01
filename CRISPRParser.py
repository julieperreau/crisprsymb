#!/usr/bin/python

# CRISPRParser.py - Written by Julie Perreau in April 2017

# This script is part 3/3 of a series:
# 1. GenomeDownload.py: Download genomes from NCBI
# 2. CRISPRFinder.py: Submit genomes to the CRISPRFinder server and download the results
# 3. CRISPRParser.py: Parse CRISPRFinder results to make a summary table (CRISPR_Counts.tsv)

# Directions for use:
# Place this script in the same directory as your CRISPRFinder.py result files (*_crispr.txt)
# Your folder should include a file named "Genome_Filenames.txt"
# 	Genome_Filenames.txt is used to associate your *_crispr.txt file names with their corresponding genome names
#	It should be formatted with this information on each line:
#	*_crispr.txt \t desired genome name \n

# Genome_Filenames.txt example:
# Bifidobacterium_asteroides_319161.fna	Bifidobacterium asteroides strain Bin7
# Frischella_perrara_237181_crispr.txt	Frischella perrara strain PEB0191
# Snodgrassella_alvi_1196083.59_contigs.2_crispr.txt	Snodgrassella alvi HK9x
	
# Usage: for i in *_crispr.txt; do python CRISPRParser.py $i EmailAddress; done
# EmailAddress is the address you use for accessing NCBI
# Output: CRISPR_Counts.tsv

import os, sys, re, time, signal
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio import Entrez
from Bio.SeqUtils import six_frame_translations

filename=sys.argv[1]
email_address=sys.argv[2]
Entrez.email = email_address

genome=filename.strip('_crispr.txt')
spacer=""
array_total=0
confirmed_total=0
questionable_total=0
arraystatus=""
genome_name=""
crisprID=""
SpacerTotal_Questionable=0
SpacerTotal_Confirmed=0

# MAKING A DESCRIPTIVE FILE OF THE CRISPR ARRAYS IDENTIFIED BY CRISPRFINDER --------------

# if not os.path.exists("CRISPR_Sequences"):
#     os.makedirs("CRISPR_Sequences")
# 
# with open(filename) as a: # Opening the *_crispr.txt results file
# 	with open("Genome_Filenames.txt", "r") as p: # Opening the genome name file
# 		for line in p:
# 			if filename in line:
# 				genome_name=line.split("\t")[1].strip("\n")
# 	crisprcount=0 # Count the number of questionable and confirmed crispr arrays
# 	length_bp=0 # Count the length of the genome
# 	scaffold_count=0 # Count the number of scaffolds in the genome
# 	spacerlist=[]
# 	sequence_dico={}
# 	for i in a:
# 		if i.startswith('Sequence:'): # Record the sequence number
# 			seqID=i.split(': ')[1].strip('\n')
# 			scaffold_count+=1
#  		if i.startswith('Sequence description :'):
#  			descriptor=i.split(': ')[1].strip("\n")
# 		if i.startswith('Confirmed CRISPRs ='): # Record # of questionable & confirmed CRISPR arrays
# 			confirmed=re.split('[= ]',i)[4].strip('\n')
# 			questionable=re.split('[= ]',i)[9].strip('\n')
# 			confirmed_total=confirmed_total+int(confirmed)
# 			questionable_total=questionable_total+int(questionable)
# 			array_total=questionable_total+confirmed_total # Add up the total # of CRISPR arrays
# 		if i.startswith('Length (bp): '):
# 			length_seq=int(i.split(":")[1].strip(" "))
# 			length_bp=length_bp+length_seq
# 		if i.startswith("Questionable "): # Record the status of the current array
# 			arraystatus="Questionable"
# 		if i.startswith("Confirmed CRISPRs ("):
# 			arraystatus="Confirmed"
# 		if i.startswith("CRISPR id"):
# 			spacerlist=[]
# 			crisprID=i.split(": ")[1].strip('\n')
# 		if i.startswith('DR consensus'):
# 			repeat=i.split(': ')[1].strip('\n') # Record the repeat sequence
# 		if i.startswith('DR length'):
# 			repeatlength=i.split(': ')[1].split(' Number')[0].strip('\n').strip(' ')
# 			numberspacers=int(i.split(': ')[2].strip('\n').strip(' '))
# 			if arraystatus=="Questionable":
# 				SpacerTotal_Questionable=SpacerTotal_Questionable+numberspacers
# 			if arraystatus=="Confirmed":
# 				SpacerTotal_Confirmed=SpacerTotal_Confirmed+numberspacers
# 			spacer=""
# 			crisprcount+=1
# 		if re.match('\d',i):
# 			b=i.split(" ")
# 			if len(b)==4:
# 				spacer=b[2] # Record the sequence
# 				spacerlist.append(spacer) # Make a list of the spacer sequences for the array
# 			else:
# 				pass
# 			if len(crisprID)>0:
# 				sequence_dico[crisprID]=spacerlist
# 			else:
# 				continue
# 
# SpacerTotal=SpacerTotal_Questionable+SpacerTotal_Confirmed
# 
# if os.path.isfile('CRISPR_Counts.tsv') == False:
# 	with open("CRISPR_Counts.tsv","w") as countfile:
# 		countfile.write("Genome\tGenome_Size\tScaffolds\tTotal_CRISPR_arrays\tTotal_CRISPR_spacers\tConfirmed_CRISPR_arrays\tConfirmed_CRISPR_spacers\tQuestionable_CRISPR_arrays\tQuestionable_CRISPR_spacers\n")
# if os.path.getsize("CRISPR_Counts.tsv") > 0:
# 	with open("CRISPR_Counts.tsv","a") as countfile:
# 		countfile.write(str(genome_name)+"\t"+str(length_bp)+"\t"+str(scaffold_count)+"\t"+str(array_total)+"\t"+str(SpacerTotal)+"\t"+str(confirmed_total)+"\t"+str(SpacerTotal_Confirmed)+"\t"+str(questionable_total)+"\t"+str(SpacerTotal_Questionable)+"\n")

# TO MAKE A FILE FOR EACH GENOME THAT CONTAINS THE SPACER SEQUENCES ----------------------


# for i in sequence_dico.iteritems():
# 	with open("./CRISPR_Sequences/"+genome+"_crispr_sequences.fna", "w") as seqfile:
# 		for crisprID in sequence_dico:
# 			spcounter=0
# 			for sp in sequence_dico[crisprID]:
# 				spcounter+=1
# 				seqfile.write(">"+str(descriptor)+", "+str(crisprID)+":"+str(spcounter)+"\n")
# 				seqfile.write(sp+"\n")


# TO BLAST AGAINST THE ONLINE NCBI NT DATABASE -------------------------------------------

e_val=10

# Make a file to store BLAST hits
if os.path.isfile('BLASTX_Hits.tsv') == False:
	with open("BLASTX_Hits.tsv","w") as countfile:
		countfile.write("Genome\tQuery\tSequence_Hit\tLength\tE-Value\n")

# BLASTX (translated nuc query to protein db)
print "Running BLAST for: "+genome
fasta_string = open("CRISPR_Sequences/"+genome+"_crispr_sequences.fna").read()
result_handle = NCBIWWW.qblast("blastx", "nr", fasta_string) # Add entrez_query="txid10239[orgn]" if searching for viruses
records= NCBIXML.parse(result_handle)
for record in records:
	for alignment in record.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < e_val:
				if os.path.getsize("BLASTX_Hits.tsv") > 0:
					with open("BLASTX_Hits.tsv","a") as countfile:
						countfile.write(str(genome)+"\t"+record.query+"\t"+str(alignment.title)+"\t"+str(alignment.length)+"\t"+str(hsp.expect)+"\n")
result_handle.close()
records.close()

# BLAST the spacer sequences against the NCBI nr database
# print "Running BLAST for: "+genome
# fasta_string = open("CRISPR_Sequences/"+genome+"_crispr_sequences.fna").read()
# result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string) # Add entrez_query="txid10239[orgn]" if searching for viruses
# records= NCBIXML.parse(result_handle)
# for record in records:
# 	for alignment in record.alignments:
# 		for hsp in alignment.hsps:
# 			if hsp.expect < e_val:
# 				if os.path.getsize("BLASTN_NT_Hits.tsv") > 0:
# 					with open("BLASTN_NT_Hits.tsv","a") as countfile:
# 						countfile.write(str(genome)+"\t"+record.query+"\t"+str(alignment.title)+"\t"+str(alignment.length)+"\t"+str(hsp.expect)+"\n")
# result_handle.close()
# records.close()
# 
# BLAST the spacer sequences against the NCBI wgs database
# fasta_string = open("CRISPR_Sequences/"+genome+"_crispr_sequences.fna").read()
# result_handle = NCBIWWW.qblast("blastn", "wgs", fasta_string)
# records= NCBIXML.parse(result_handle)
# for record in records:
# 	for alignment in record.alignments:
# 		for hsp in alignment.hsps:
# 			if hsp.expect < e_val:
# 				if os.path.getsize("BLASTN_WGS_Hits.tsv") > 0:
# 					with open("BLASTN_WGS_Hits.tsv","a") as countfile:
# 						countfile.write(str(genome)+"\t"+record.query+"\t"+str(alignment.title)+"\t"+str(alignment.length)+"\t"+str(hsp.expect)+"\n")
# result_handle.close()
# records.close()
