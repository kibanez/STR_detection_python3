This repository contains scripts and tools for analyzing RepeatExpansions of short tandem repeats (STRs).

# STR_pileup.py:
Visualization of the pileup of reads aligning to an STR locus for one or more samples.
## Description:
  Genotyping STR loci is a difficult problem due to multiple reasons including flanks with similar genomic sequences, difficulty in aligning reads sequencing errors and sample contamination. ExpansionHunter is a tool that can automatically predict the genotypes of STR loci based on read alignments. However, due to the difficult nature of the problem, there may be errors in the predicted genotype. STR_pileup.py, is a standalone Python script that creates a visualization of all read alignments to an STR locus to enable visual inspection of the genotype predictions made by ExpansionHunter. The script can visualize one or more samples, including trios to compare genotype calls in multiple samples. For each sample the script requires genotypes predicted by ExpansionHunter in JSON format and the read alignments file (log file) reported by ExpansionHunter in YAML format.
## Usage:
`python3 STR_pileup.py --help`
### Use case 1: Single sample
`python3 STR_pileup.py --json JSON_FILE --read_align READ_ALIGN_FILE [--color] [--repeat_id REPEAT_ID]`
### Use case 2: Comparison of genotypes from multiple samples
`python3 STR_pileup.py --read_align_list READ_ALIGN_FILE_LIST [--color] [--repeat_id REPEAT_ID]`
## Requirements:
* Python3
* Matplotlib
## Inputs:
### Use case 1: Single sample
1. `--json FILE`: JSON output of genotype calls from ExpansionHunter
2. `--read_align FILE`: Read alignment file (log file) generated be ExpansionHunter in YAML format.
### Use case 2: Comparison of genotypes from multiple samples
1. `--read_align_list FILE`: A 3 column text (or tsv) file containing the list of EH output for all samples:
* Column1: Sample name,
* Column2: JSON file from EH output,
* Column3: Read alignment (log file) from EH output  


The pileups for the samples are shown in the order reported in the input file.  
**Note:** If the paths to the JSON and log file are not absolute, they are chosen relative to the location of the sample list file.
### Additional options:
[Optional] `--repeat_id REPEAT_ID`: Create pileup only for STR locus ID matching the repeat_id in the JSON file. Default: plots pileups for all loci in JSON file.  
[Optional] `--color`: Flag stating that each base-pair should be colored in IGV format. Default: Matching bases in black, mismatching bases in red.
## Output:
2 images per locus: \<PREFIX\>\_\<REPEAT_ID\>.png and \<PREFIX\>\_\<REPEAT_ID\>.pdf.  
The output files are created in the current working directory. PREFIX matches the filename of the input JSON file or the sample list file.  
Specifically:
### Use case 1: Single sample
If the input JSON is "/\<path to parent directory\>/sample.json" and REPEAT_ID is C9ORF72, then the output files will be produced in the current working directory with file names sample_C9ORF72.pdf and sample_C9ORF72.png.
### Use case 2: Comparison of genotypes from multiple samples
If the input JSON is "/\<path to parent directory\>/samplelist.txt" or "/\<path to parent directory\>/samplelist" and the REPEAT_ID is C9ORF72, then the output files will be produced in the current working directory with file names samplelist_C9ORF72.pdf and samplelist_C9ORF72.png.
