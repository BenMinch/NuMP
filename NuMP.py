###NuMP: Nucleocytovirotica Metabolic Profiling###

import pandas as pd
import numpy as np
import sys, os, argparse, re, subprocess, shutil, glob

argparser = argparse.ArgumentParser(description='''NuMP: Nucleocytovirotica Metabolic Profiling''')
argparser.add_argument('-i', '--input', help='Input directory containing the .fasta genome files', required=True)
argparser.add_argument('-o', '--output', help='Output directory', required=True)
argparser.add_argument('-t', '--taxonomy', help='Taxonomy file with column for genome and one for family (must be tsv)', required=True)
argparser.add_argument('-m', '--mode', help='Mode: "all" for all genomes, "family" for grouping by family classification, or "both"', required=True)

args = argparser.parse_args()

#Create output directory
if os.path.exists(args.output):
    shutil.rmtree(args.output)
os.makedirs(args.output)

input_dir = args.input
taxonomy_file = args.taxonomy
output_dir = args.output
mode = args.mode

#Prodigal to predict proteins
#make prodigal directory
prodigal_dir = output_dir + '/prodigal'
if not os.path.exists(prodigal_dir):
    os.makedirs(prodigal_dir)
prodigal= 'python scripts/prodigal_launcher.py '+ input_dir + ' ' + prodigal_dir
subprocess.call(prodigal, shell=True)

#Add filename to every fasta header
for filename in os.listdir(prodigal_dir):
    if filename.endswith('.faa'):
        filepath = os.path.join(prodigal_dir, filename)
        with open(filepath, 'r') as file:
            lines = file.readlines()
        with open(filepath, 'w') as file:
            for line in lines:
                if line.startswith('>'):
                    header = '>' + filename + '_' + line[1:]
                    file.write(header)
                else:
                    file.write(line)



#Combine all proteins into one file
cat_all= 'cat ' + prodigal_dir + '/*.faa > ' + output_dir + '/all_proteins.faa'
subprocess.call(cat_all, shell=True)

#Run Annomazing on the proteins 
annomazing_dir = output_dir + '/annomazing'
if not os.path.exists(annomazing_dir):
    os.makedirs(annomazing_dir)

annomazing= 'python scripts/AnnoMazing.py -i ' + output_dir + '/all_proteins.faa -o annomazing -db hmm/Pfam.hmm -annot resources/pfam_annotations.csv'
subprocess.call(annomazing, shell=True)


##Run the R script to parse the file and graph it
r_script_genome = 'Rscript scripts/metabolicplotter.r annomazing_final.csv ' + output_dir + ' '+ mode + ' ' + taxonomy_file
subprocess.call(r_script_genome, shell=True)

#move annomazing file to output directory
move= 'mv annomazing_final.csv ' + output_dir
subprocess.call(move, shell=True)