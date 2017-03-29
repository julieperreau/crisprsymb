#!/usr/bin/python

# Script for submitting genomes from a folder to CRISPRFinder (http://crispr.i2bc.paris-saclay.fr/Server/)
# Saves the .txt output of a CRISPRFinder search to a new folder, titled "CRISPR Results"
# Julie Perreau - March 8 2017

# Example usage: for i in *.txt; do python CRISPRFinder.py

import os, sys
from selenium import webdriver
from selenium.webdriver.support import ui
from selenium.webdriver.support.wait import WebDriverWait    
from selenium.webdriver.common.keys import Keys

fastafile = sys.argv[1]
genome = fastafile.strip(".fna")
browser = webdriver.Chrome(executable_path = '/Users/JuliePerreau/scripts/chromedriver')
if not os.path.exists("CRISPR_Results"): # Create a folder for downloaded genomes
 	os.makedirs("CRISPR_Results")
fullname=os.getcwd()+"/CRISPR_Results/"
local_file=fullname+genome+"_crispr.txt"

def page_is_loaded(driver):
    return driver.find_element_by_tag_name("body") != None

driver = webdriver.Chrome("/Users/JuliePerreau/scripts/chromedriver")
driver.get("http://crispr.i2bc.paris-saclay.fr/Server/")

wait = ui.WebDriverWait(driver, 10)
wait.until(page_is_loaded)

driver.find_element_by_name("fname").send_keys(os.getcwd()+'/'+fastafile)
driver.find_element_by_name("submit").click()
results = driver.find_element_by_xpath("//div[@class='content']").text

with open(local_file, 'w+') as t:        
    t.write('\n'.join(results.split("\n")[3:]))

t.close()	
driver.quit()
browser.quit()