#!/usr/bin/python

# CRISPRFinder.py - Written by Julie Perreau in April 2017

# This script is part 2/3 of a series:
# 1. GenomeDownload.py: Download genomes from NCBI
# 2. CRISPRFinder.py: Submit genomes to the CRISPRFinder server and download the results
# 3. CRISPRParser.py: Parse CRISPRFinder results to make a summary table (CRISPR_Counts.tsv)

# This script submits genomes from your directory to CRISPRFinder (http://crispr.i2bc.paris-saclay.fr/Server/)
# and saves the .txt output to a new directory called "CRISPR Results"

# Directions for use:
# Must have this script in the same directory as your genome files (*.fna)
# Must have chromedriver installed and know the chromedriver path on your computer
#	Can install with homebrew package manager: brew install chromedriver

# Usage: for i in *.txt; do python CRISPRFinder.py path/to/chromedriver; done
# Example chromedriver path: '/Users/Username/scripts/chromedriver'

import os, sys
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.by import By
from selenium.webdriver.support import ui
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException

fastafile = sys.argv[1]
chromedriverpath = sys.argv[2]
genome = fastafile.strip(".fna")
browser = webdriver.Chrome(executable_path = chromedriverpath)
if not os.path.exists("CRISPR_Results"): # Create a folder for downloaded genomes
 	os.makedirs("CRISPR_Results")
fullname=os.getcwd()+"/CRISPR_Results/"
local_file=fullname+genome+"_crispr.txt"

driver = webdriver.Chrome(chromedriverpath)
driver.get("http://crispr.i2bc.paris-saclay.fr/Server/")

driver.find_element_by_name("fname").send_keys(os.getcwd()+'/'+fastafile)
driver.find_element_by_name("submit").click()

element = WebDriverWait(browser, 30)
results = driver.find_element_by_tag_name("body").text

with open(local_file, 'w+') as t:        
    t.write('\n'.join(results.split("\n")[3:]))

t.close()
driver.quit()
browser.quit()
