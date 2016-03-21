#!/usr/bin/env python
#
# Annotate the sequences. 
#

import argparse
import subprocess
import sys
import shlex

SIZE = 200


iupac = {'A_G': 'R',
      'G_A': 'R',
      'C_T': 'Y',
      'T_C': 'Y',
      'G_C': 'S',
      'C_G': 'S',
      'G_T': 'K',
      'T_G': 'K',
      'A_C': 'M',
      'C_A': 'M',
      'A_T': 'W',
      'T_A': 'W'
    }

def get_region(position, size=200):
    region = str(int(position) - size) + "-" + str(int(position)+ size)
    return region

def extract_snps(all_vcf, chromosome, position):
    """
            Extract SNPs flanking the regions of the SNP.
    """
    region = get_region(position)

    command = ["bcftools", "view", all_vcf, chromosome + ":" + region]
    #print command
    try:
        output = subprocess.check_output(command)
    except subprocess.CalledProcessError:
        sys.stderr.write("Could not run bcftools\n")
        sys.exit(1)
    with open(position + ".vcf", 'w') as f:
        f.write(output)
    return output
    

def get_sequence(fasta, chromosome, position, alt_allele):
    """
            Extrat sequence surrounding the SNP.
    """
    region = str(int(position) - 200)
    if "chr" in chromosome:
        chromosome = chromosome.split('chr')[1]
    length = "401"
    command= """bioawk -c fastx '{{ if("{0}" == $name){{ print substr($seq,{1}, {2})}}}}' {3}""".format(chromosome, region, length, fasta)
    try:
        output = subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError:
        sys.stderr.write("Could not run bioawk\n")
        sys.exit(1)
    output = output.strip()
    return(output) 

def change_table(ref, alt):
    combine = ref + "_" + alt
    code = iupac[combine]
    return code

def annotate_nearby(all_vcf, position_snp, sequence):
    """
        Annotate the nearby SNPs and indels
    """
    sequence = list(sequence)
    for line in all_vcf.split('\n'):
        if "#" not in line:
            if line == '':
                break
            l_s = line.split('\t')
            vcf_position = int(l_s[1])
            ref = l_s[3]
            alt = l_s[4]
            # TODO check this is consistent with the results in the 1-indexed and 0-indexed lists.
            #print vcf_position
            print(position_snp)
            position =  int(vcf_position)  - (int(position_snp) - 200)
            #print sequence[position]
            print(ref, alt)
            if len(ref) > 1 or len(alt) > 1:
                alt_s = alt.split(',')
                alt_s = [len(o) for o in alt]
                idx = alt_s.index(min(alt_s))
                print(idx)
                length_alt_s = alt_s[idx]
                print(length_alt_s)
                dist = len(ref) - length_alt_s 
                if dist > 0:
                    for i in range(length_alt_s, length_alt_s + dist):
                        sequence[position + i] = 'N' 
                else:
                    print('wow')
            else:
                change_snp = change_table(ref, alt)
                print(change_snp)
                print(sequence[position])
                if position == 200:
                    continue
                sequence[position] =  change_snp

    return ''.join(sequence)

def finalise_design(output, alt_allele):
    output = output[0:200] + '[' + output[200] + '/' + alt_allele + ']' + output[201:]
    return output

def write_snp_design(sequence, position, chromosome):
    """
        Write the SNP design
    """
    with open(chromosome + "." + position, 'w') as f:
        f.write("> " + chromosome + ":" + position +'\n')
        f.write(sequence)
        f.write('\n')

def main():
    """
            Annotate the sequence for geneworks.
    """
    parser = argparse.ArgumentParser(description="Annotate sequence for geneworks.")
    parser.add_argument("all_vcf")
    parser.add_argument('-f','--fasta', dest='fasta', help='Fasta file', required=True)
    parser.add_argument('-c','--chr', dest="chromosome", help="Chromosome to extract from", required=True)
    parser.add_argument('-p','--pos', dest="position", help="Position of the SNP", required=True)
    parser.add_argument('-a','--alt', dest='alternate_allele', help='What the alternate allele is', required=True)
    args = parser.parse_args()
    subset_vcf = extract_snps(args.all_vcf, args.chromosome, args.position)
    sequence = get_sequence(args.fasta, args.chromosome, args.position, args.alternate_allele)
    sequence = annotate_nearby(subset_vcf, args.position, sequence)
    sequence = finalise_design(sequence, args.alternate_allele)
    write_snp_design(sequence,args.position, args.chromosome)  


if __name__=="__main__":
    main()
