# Annotate a sequence segment with SNPs/variants
## This script is intended for use when designing custom TaqMan probes or Sequenom / GeneWorks assays.

## To install this script package:

 - Download a .zip file of this repository
 - To install the package unzip the file (put it somewhere save eg. home/Executables) and cd into it (in Terminal)
 - In Terminal type:

```
        python2.7 setup.py install --user
```

 - Check that python is also in your bash PATH
    - In Terminal type:

```
        nano ~/.bashrc
```
    - In the window that opens check there is a line that says:

```
        export PATH=$PATH:/Users/[yourhomefolder]/Library/Python/2.7/bin
```
    - Save any edits you make in this window then type:

```
        ln -s ~/.bashrc ~/.bash_profile
        source ~/.bashrc
```

 - Now this script is set up as a command line program for you to use!

## To run this script:

 - You will need:
    - A .vcf.gz file to extract SNPs/variants of interest from and the corresponding index file (.vcf.gz.tbi) 
    - A .fasta file of the same genome build your vcf uses to specify SNP/variant positions
    - The chromosome and position of your SNP of interest
    - The base corresponding to the alternate allele

 - Command arguments:
    - Required arguments:
```
        -v, --vcf VCF_FILE
                        a bgzipped VCF file (extension .vcf.gz), this file
                        also requires an index file (extension .vcf.gz.tbi)
        -f, --fasta FASTA
                        Fasta file
        -c, --chr CHROMOSOME
                        Chromosome to extract region from; check whether your
                        vcf uses 'chr1' or '1' formatting
        -p, --pos POSITION_SNP
                        Position of the SNP 
        -a, --alt ALTERNATE_ALLELE
                        What the alternate allele is
```
    - Optional arguments:

```
        -m, --mask MASK_CHARACTER
                        Do you want nearby SNPs labelled with iupac codes or
                        Ns?, Default = 'iupac', Alternate = 'N'
        -s, --size SIZE_REGION
                        How much sequence either side of you SNP do you want?
                        Default = 200
        -o, --output OUTPUT_FILE
                        What do you want to call the output file name? Default
                        = chr.pos
```

 - Example Usage and Output:

```
annotate_sequence --vcf polyReSeq.vcf.gz --fasta human_build37.fasta --chr chr11 --pos 64323183 --alt T --size 300 --output 11_64323183_poly

+ - - - - - - - - - - - - - - - - - - - - +
           annotate_sequence-v2.0          
+ - - - - - - - - - - - - - - - - - - - - +

Script run on: 21/03/16 17:15:50

Arguments in use:

        --vcf polyReSeq.vcf.gz
        --fasta human_build37.fasta
        --chr chr11
        --pos 64323183
        --alt T
        --mask iupac
        --size 300
        --output 11_64323183_poly

4 SNPs/variants were found within the extracted region.

chr11:64323080; [C/T]
chr11:64323312; [A/C]
chr11:64322891; [G/C]
chr11:64323183; [C/T]

annotate_sequence-v2.0 was successfully run.

    A file called 11_64323183_poly.fasta was created.
    A file called 11_64323183_poly.log was created.
    A file called 11_64323183_poly.vcf was created.
```
