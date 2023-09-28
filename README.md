# NuMP: Nucleocytovirotica Metabolic Profiler
A program to quickly and thouroughly visualize the metabolic and functional potential of Nucleocytovirotica genomes in a dataset. NuMP is based off of analysis done in Ha, A. D., Moniruzzaman, M., & Aylward, F. O. (2021). High transcriptional activity and diverse functional repertoires of hundreds of giant viruses in a coastal marine system. MSystems, 6(4), e00293-21.

## Dependencies
1. HMMER (v3.3.2)
2. Python (v3+) with pandas, numpy, and Bio
3. R(v4.3.1) with tidyverse installed
4. Prodigal

## Installation

Simply clone this repository to get the program up and running.
`git clone https://github.com/BenMinch/NuMP` 

### Getting the database
NuMP uses the PFAM hmm database. You will need to download that https://www.ebi.ac.uk/interpro/download/Pfam/ . Make sure you download the PFAM-A-models file. After downloading it, rename the file Pfam.hmm and store the database in a folder called 'hmm'. You will also need to press the file with hmmpress `hmmpress Pfam.hmm`. 

# Running the Program

### Inputs
1. -i: Input directory with genome files in fasta format.
2. -o: Name of desired output directory
3. -t: A tab-separated taxonomy file with a column for genome and one for family. These two names must be on the file.
4. -m: Mode selection. Can select "all" for each genome to be profiled individually, or "family" for the annotations to be grouped by family. You can also choose "both" here.

### Outputs
1. Annomazing_final.csv: All protein annotations for all of your genomes.
2. metabolic_data.csv: All of the selected metabolic genes with information from which genome it came from, the taxonomy, the gene, and the pathway. The protein name is also present in query_id so you can easily recover the gene and do further read mapping with this data.
3. bubble_family/genomes.pdf: A bubble plot with the functional profiling of your genomes or families.
4. Prodigal folder: all protein and gene predictions for your genomes.

### Example run

`python NuMP.py -i Genomes -o Genomes_NuMP -t genome_taxonomy.tsv -mode both`

# Copywright
NuMP Copyright (C) 2023 Benjamin Minch

This program is free software: you can redistribute it and/or modify it under the terms of the MIT License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT License for more details.
