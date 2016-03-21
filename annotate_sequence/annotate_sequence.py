#!/usr/bin/env python
#
# Annotate a sequence with variants from a .vcf file. 
#

import argparse
import subprocess
import sys
import shlex
import datetime

iupac = {
    "A_G": "R", "G_A": "R",
    "C_T": "Y", "T_C": "Y",
    "G_C": "S", "C_G": "S",
    "G_T": "K", "T_G": "K",
    "A_C": "M", "C_A": "M",
    "A_T": "W", "T_A": "W"
    }

snp_dictionary = {}

def check_system(command):
    return subprocess.call("type " + command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE) == 0


def get_region(position_snp, size_region):
    region = str(int(position_snp) - int(size_region)) + "-" + str(int(position_snp) + int(size_region))
    return region


def extract_snps(vcf_file, chromosome, position_snp, size_region, output_file):
    """
            Extract SNPs in the region flanking the SNP of interest.
    """
    region = get_region(position_snp, size_region)

    command = ["bcftools", "view", vcf_file, chromosome + ":" + region]

    try:
        output = subprocess.check_output(command, stderr = subprocess.PIPE)
    except subprocess.CalledProcessError:
        sys.stderr.write("\nError: Could not run bcftools\n" + 
                         "    Are you using a bgzipped vcf? (extension .vcf.gz)\n" + 
                         "    Have you made an index file for your vcf? (extension .vcf.gz.tbi)\n")
        sys.exit(1)
    with open(output_file + ".vcf", "w") as f:
        f.write(output)
    return output
    

def get_sequence(fasta, chromosome, position_snp, alt_allele, size_region):
    """
            Extract sequence surrounding the SNP.
    """
    region = str(int(position_snp) - int(size_region))
    if "chr" in chromosome:
        chromosome = chromosome.split("chr")[1]
    length = str((int(size_region) * 2) + 1)
    command = """bioawk -c fastx '{{ if("{0}" == $name){{ print substr($seq,{1}, {2})}}}}' {3}""".format(chromosome, region, length, fasta)
    try:
        output = subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError:
        sys.stderr.write("\nError: Could not run bioawk\n")
        sys.exit(1)
    output = output.strip()
    return output  


def change_table(ref, alt, mask_character):
    combine = ref + "_" + alt
    code = iupac[combine]
    if mask_character == "iupac":
        return code
    else:
        return "N"


def annotate_nearby(subset_vcf, position_snp, sequence, mask_character, size_region):
    """
        Annotate the nearby SNPs and indels
    """
    sequence = list(sequence)
    for line in subset_vcf.split("\n"):
        if "#" not in line:
            if line == "":
                break
            line_split = line.split("\t")
            vcf_position = int(line_split[1])
            ref = line_split[3]
            alt = line_split[4]
            position = int(vcf_position) - (int(position_snp) - int(size_region))
            snp_dictionary[int(vcf_position)] = [ref, alt]
            if len(ref) > 1 or len(alt) > 1:
                alt_s = alt.split(",")
                alt_s = [len(o) for o in alt]
                idx = alt_s.index(min(alt_s))
                length_alt_s = alt_s[idx]
                dist = len(ref) - length_alt_s 
                if dist > 0:
                    for i in range(length_alt_s, length_alt_s + dist):
                        sequence[position + i] = "N" 
            else:
                change_snp = change_table(ref, alt, mask_character)
                if position == size_region:
                    continue
                sequence[position] = change_snp

    return "".join(sequence)


def finalise_design(output, alt_allele, size_region):
    output = output[0:int(size_region)] + "[" + output[int(size_region)] + "/" + alt_allele + "]" + output[int(size_region)+1:]
    return output

def write_files(vcf_file, fasta, chromosome, position_snp, alternate_allele,
                mask_character, size_region, output_file, sequence):
    """
        Write out the result files
    """
    run_date_time = datetime.datetime.strftime(datetime.datetime.now(), '%d/%m/%y %H:%M:%S')
    region = get_region(position_snp, size_region)
    snps_found = len(snp_dictionary)
    """
        Write a log file of the arguments and outputs of this script
    """
    with open(output_file + ".log", "w") as f:
        f.write("+ - - - - - - - - - - - - - - - - - - - - +\n" + 
                "           annotate_sequence-v2.0          \n" + 
                "+ - - - - - - - - - - - - - - - - - - - - +\n" + 
                "\n")
        f.write("Script run on: " + run_date_time + "\n")
        f.write("Arguments in use:\n" + 
                "        --vcf " + str(vcf_file) + "\n" + 
                "        --fasta " + str(fasta) + "\n" + 
                "        --chr " + str(chromosome) + "\n" + 
                "        --pos " + str(position_snp) + "\n" + 
                "        --alt " + str(alternate_allele) + "\n" + 
                "        --mask " + str(mask_character) + "\n" + 
                "        --size " + str(size_region) + "\n" + 
                "        --output " + str(output_file) + "\n")
        f.write("\n")
        f.write("The region within " + chromosome + ":" + region + " was extracted from the " + vcf_file + "file.\n" + 
                "   A new file called " + chromosome + ":" + position_snp + ".vcf was created.\n")
        f.write("\n")
        f.write(str(snps_found) + " SNPs/variants were found within the extracted region.\n")
        for snp,alleles in snp_dictionary.items():
            f.write("\n" + str(chromosome) + ":" + str(snp) + ";")
            f.write(" [" + alleles[0] + "/" + alleles[1] + "]\n")
        f.write("\n")
        f.write("+ - - - - - - - - - - - - - - - - - - - - +\n\n")
        f.write("Sequence Output:\n")
        f.write("\n")
        f.write("> " + chromosome + ":" + position_snp +"\n")
        f.write(sequence)
        f.write("\n")

    """
        Write a fasta format file of the annotated sequence
    """
    with open(output_file + ".fasta", "w") as f:
        f.write("> " + chromosome + "_" + position_snp + "\n")
        f.write(sequence)
        f.write("\n")


def main():
    """
            Annotate the sequence with SNPs from a vcf file - ready for custom TaqMan design / GeneWorks.
    """
    parser = argparse.ArgumentParser(description = "Annotate a sequence with SNPs from a vcf file - ready for custom TaqMan design / GeneWorks.")
    parser.add_argument("-v", "--vcf", dest = "vcf_file", help = "a bgzipped VCF file (extension .vcf.gz), this file also requires an index file (extension .vcf.gz.tbi)", required = True)
    parser.add_argument("-f", "--fasta", dest = "fasta", help = "Fasta file", required = True)
    parser.add_argument("-c", "--chr", dest = "chromosome", help = "Chromosome to extract region from; check whether your vcf uses 'chr1' or '1' formatting", required = True)
    parser.add_argument("-p", "--pos", dest = "position_snp", help = "Position of the SNP", required = True)
    parser.add_argument("-a", "--alt", dest = "alternate_allele", help = "What the alternate allele is", required = True)
    parser.add_argument("-m", "--mask", dest = "mask_character", help = "Do you want nearby SNPs labelled with iupac codes or Ns?, Default = 'iupac', Alternate = 'N'", required = False)
    parser.add_argument("-s", "--size", dest = "size_region", help = "How much sequence either side of you SNP do you want? Default = 200bp")
    parser.add_argument("-o", "--output", dest = "output_file", help = "What do you want to call the output file name? Default = chr.pos") 
    args = parser.parse_args()
    run_date_time = datetime.datetime.strftime(datetime.datetime.now(), '%d/%m/%y %H:%M:%S')
    if args.output_file == None:
        args.output_file = args.chromosome + "." + args.position_snp
    if args.size_region == None:
        args.size_region = 200
    else:
        args.size_region = int(args.size_region)
    if args.mask_character == None:
        args.mask_character = "iupac"
    sys.stderr.write("\n" + 
                     "+ - - - - - - - - - - - - - - - - - - - - +\n" + 
                     "           annotate_sequence-v2.0          \n" + 
                     "+ - - - - - - - - - - - - - - - - - - - - +\n")
    sys.stderr.write("Script run on: " + run_date_time + "\n")
    sys.stderr.write("Arguments in use:\n" + 
                     "        --vcf " + str(args.vcf_file) + "\n" + 
                     "        --fasta " + str(args.fasta) + "\n" + 
                     "        --chr " + str(args.chromosome) + "\n" + 
                     "        --pos " + str(args.position_snp) + "\n" + 
                     "        --alt " + str(args.alternate_allele) + "\n" + 
                     "        --mask " + str(args.mask_character) + "\n" + 
                     "        --size " + str(args.size_region) + "\n" + 
                     "        --output " + str(args.output_file) + "\n")
    sys.stderr.write("\n")
    assert args.mask_character == "iupac" or args.mask_character == "N", "\n--mask requires either iupac or N as the input argument.\n"
    exists_bcftools = check_system("bcftools")
    assert exists_bcftools == True, "\nError: Could not run bcftools, please install bcftools or add it to your PATH.\n"
    exists_bioawk = check_system("bioawk")
    assert exists_bioawk == True, "\nError: Could not run bioawk, please install bioawk or add it to your PATH.\n"
    subset_vcf = extract_snps(args.vcf_file, args.chromosome, args.position_snp, args.size_region, args.output_file)
    sequence = get_sequence(args.fasta, args.chromosome, args.position_snp, args.alternate_allele, args.size_region)
    sequence = annotate_nearby(subset_vcf, args.position_snp, sequence, args.mask_character, args.size_region)
    sys.stderr.write("\n" + str(len(snp_dictionary)) + " SNPs/variants were found within the extracted region.\n")
    assert len(snp_dictionary) > 0, "\nError: no SNPs were found in the vcf file;\n    - check you have specified the right SNP position\n    - check you have used the right build version\n    - check how your vcf file specifies chromsome ID (chr1 vs. 1)\n"
    for snp,alleles in snp_dictionary.items():
        sys.stderr.write("\n" + str(args.chromosome) + ":" + str(snp) +
                         ";" + " [" + str(alleles[0]) + "/" + str(alleles[1]) + "]")
    sys.stderr.write("\n")
    sequence = finalise_design(sequence, args.alternate_allele, args.size_region)
    write_files(args.vcf_file, args.fasta, args.chromosome, args.position_snp, args.alternate_allele,
                args.mask_character, args.size_region, args.output_file, sequence)
    sys.stderr.write("\nannotate_sequence-v2.0 was successfully run.\n" + 
                     "    A file called " + str(args.output_file) + ".fasta was created.\n" + 
                     "    A file called " + str(args.output_file) + ".log was created.\n" + 
                     "    A file called " + str(args.output_file) + ".vcf was created.\n" + 
                     "\n" + 
                     "+ - - - - - - - - - - - - - - - - - - - - +\n")


if __name__=="__main__":
    main()