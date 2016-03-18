#!/usr/bin/env python
#
# Test system setup
#
import datetime
import subprocess

def check_system(command):
    return subprocess.call("type " + command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE) == 0


def main():
    run_date_time = datetime.datetime.strftime(datetime.datetime.now(), '%d/%m/%y %H:%M:%S')
    print("+ - - - - - - - - - - - - - - - - - - - - +\n" + 
          "           annotate_sequence-v2.0          \n" + 
          "+ - - - - - - - - - - - - - - - - - - - - +\n")
    print("Script run on: " + run_date_time + "\n")
   # print("Arguments in use:\n" + 
   #       "        --vcf " + vcf_file + "\n" + 
   #       "        --fasta " + fasta + "\n" + 
   #       "        --chr " + chromosome + "\n" + 
   #       "        --pos " + position_snp + "\n" + 
   #       "        --alt " + alternate_allele + "\n" + 
   #       "        --mask " + mask_character + "\n" + 
   #       "        --size " + size_region + "\n" + 
   #       "        --output " + output_file + "\n")
    print("\n")
    exists_bcftools = check_system("bcftools")
    assert exists_bcftools == True, "\n\nError: Could not run bcftools, please install bcftools or add it to your PATH.\n"
    exists_bioawk = check_system("bioawk")
    assert exists_bioawk == True, "\n\nError: Could not run bioawk, please install bioawk or add it to your PATH.\n"


main()
