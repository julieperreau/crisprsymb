#!/usr/bin/python

# Script for downloading genomes from NCBI
# Julie Perreau - March 8 2017

# Example command: python GenomeDownload.py "Gilliamella apicola"

from __future__ import division
import sys, csv, os, itertools, operator, numpy as np
from collections import Counter

from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord

spp=sys.argv[1] # Input the species name you want to search for
sppname=spp.strip("\"")

Entrez.email = "jmaperreau@gmail.com"
Entrez.tool = "MyLocalScript"

if not os.path.exists("Downloaded_Genomes/"): # Create a folder for downloaded genomes
 	os.makedirs("Downloaded_Genomes/")

handle = Entrez.read(Entrez.esearch(db="assembly", term=sppname+'[Organism]',retmax=10000))
genomeId = handle['IdList'][0:] # Puts the GenBank identifiers for every assembly into a list 

for i in genomeId:
	linker = Entrez.read(Entrez.elink(db="nucleotide",dbfrom="assembly",id=i,linkname="assembly_nuccore_insdc")) # elink links the assembly ID to the nuccore IDs for the individual genomes
	numbIDs=len(linker[0]["LinkSetDb"][0]["Link"]) # Counts the number of contigs
	for j in range(0,numbIDs):
		try:
			GI_ID = linker[0]["LinkSetDb"][0]["Link"][j]["Id"]
			record = Entrez.efetch(db="nucleotide", id=GI_ID, rettype="fasta", retmode="text")
 			fullname=os.getcwd()+"/Downloaded_Genomes/"
			sppname=sppname.strip("\n").replace(".","").replace("=","").replace("+","").replace(" ", "_").replace("-","_")
			local_file=open(fullname+sppname+"_"+i+".fna", 'a')
			local_file.write(record.read())
			record.close()
			local_file.close()
		except IndexError:
			continue