'''
This workflow is for mapping fastq files to a fasta 
file with genes in to calculate coverage.
The reads are mapped with BWA and and filtered.
The coverage for each gene is calculated are compared 
for each IND.


Author: CMK (heavily inspired from Elise Lucotte and Laurits Skov)
Date: may. 2018 

Changelog: 
feb     20th    cleanup
feb     23th    follow present style
feb     24th    move templates to seperate file,  upgrade gwf to 1.2.1 and python to v3
feb     27th    filter updated
mar 14th                new parametric setup with a list of dictionaries and no controller
mar 18th        successfull batching


 ~ todo ~
* record timeline/core/memusage etc.
'''
import sys, os
from os.path import basename # More paths?
from gwf import Workflow
from workflow_templates import *
import json # For pretty printing trees
import subprocess # For running bash in python
from datetime import datetime # To set dates in titles
import re # to escape symbols in strings


# Define the root gwf object
gwf = Workflow()


def get_individuals(path):
        '''
        Takes a path leading to a sample_info.txt-file
        and returns a list of individuals
        '''
        with open(path, 'r') as file:
                return [line.split('\t')[0] for line in file][1:]
                #                                                ^ We are only interested in the names of the individuals (first column)
                #                                                                                         ^ We are only interested in the rows after the labels (labels are on row 0)


def get_filenames(individuals, path):
        '''
        Takes a list of individuals
        and returns a dictionary for pairs of genome-files
        '''
        dico_names = {} # dictionary for unique individuals
        for individual in individuals:  # for each individual in the sample_info.txt-file:
                dico_names[individual] = []     # create an empty list for each individual
                # file = working_dir+'../Names_files/' + individual + '.txt'    # path to the link-file (every individual has its own file)
                with open(path+individual+'.txt') as file:      # close afterwards
                        lines = [line.strip() for line in file] # transfer the strippedlines to a list, so we can later pair them together. The stripping removes newlines present in the end of each line.
                        for pair in zip(lines[0::2], lines[1::2]):      # zipping the strided lines, puts them in pairs [(1,2), (2,3), ...]
                                dico_names[individual].append(pair)     # append the file pairs to the list in the dictionary for each individual
        return dico_names


# This is a list of batches. Each batch is i dictionary with various keys and their values.
# KEY                                           VALUE
# title                                         The name of the batch
# chromosome                            The name of the chromosome being worked on
# rel_reference                         The relative path to the reference (artificial chromosome)
# rel_sample_info_file          The relative path to the sample info file, containing a list of individuals
# rel_individual_genomes        The relative path to the directory containing a file for each individual, each file containing a list og filenames which are the haploid genome builds
# abs_genome_dir                        The absolute path to the directory containing the files listed in the rel_individual_genomes directory's files.


# chimp
batches = [{"title": "batch_x6_chimp",
        "chromosome": "x",
        "description": "k\u00f8rsel fredag d 30. apr med ar. chrom. for chimp",
        "rel_reference": "ac_chimp_x.fa",
        "rel_sample_info_file": "sample_info_chimp.txt",
        "rel_individual_genomes": "genomes/",
        "abs_genome_dir": "/home/cmkobel/MutationRates/faststorage/NEW_PIPELINE/TrimFASTQ/TrimmedFASTQ/"},
        
        {"title": "batch_x6.3_gorilla",
        "chromosome": "x",
        "description": "gorilla med chimps opsætning ændret",
        "rel_reference": "ac_gorilla_x.fa",
        "rel_sample_info_file": "sample_info_gorilla.txt",
        "rel_individual_genomes": "genomes/",
        "abs_genome_dir": "/home/cmkobel/MutationRates/faststorage/NEW_PIPELINE/TrimFASTQ/TrimmedFASTQ/"},
        
        {"title": "batch_y6.2_chimp",
        "chromosome": "y",
        "description": "ændret opsætning til chimp y",
        "rel_reference": "ac_chimp_y.fa",
        "rel_sample_info_file": "sample_info_chimp.txt",
        "rel_individual_genomes": "genomes/",
        "abs_genome_dir": "/home/cmkobel/MutationRates/faststorage/NEW_PIPELINE/TrimFASTQ/TrimmedFASTQ/"}]


# Protocol:
for batch in batches: # Måske batch bare skal hedde b ?
        
        # if batch['chromosome'] == 'x':
        #       batch['rel_reference'] ='artificial_chr_x.fa'
        # elif batch['chromosome'] == 'y':
        #       batch['rel_reference'] ='artificial_chr_y.fa'

        individuals = get_individuals(batch['rel_sample_info_file'])
        dico_names = get_filenames(individuals, batch['rel_individual_genomes'])

        # 0: initialize
        gwf.target_from_template(batch['title']+'_0_initialize', initialize(batch['title'], json.dumps(batch, sort_keys=False, indent=4)))

        # 1: index
        gwf.target_from_template(batch['title'] + '_1_index', index_genome(batch['title'], batch['rel_reference']))

        for individual in dico_names:
                BAM_files=[] # for collecting BAM-files later.
                for num, pair in enumerate(dico_names[individual]):
                        # 2: Mapping
                        gwf.target_from_template(batch['title'] + '_2_map_'+individual+str(num), bwa_map_pe(batch['title'], batch['rel_reference'], batch['abs_genome_dir']+pair[0], batch['abs_genome_dir']+pair[1], individual+str(num)))
                        BAM_files.append(individual+str(num)+'_sort_dedup.bam') # collect the filenames for use in the merge point
                print('BAM_files:', BAM_files)

                # 3: Merging the bam files # merge the bam-files for each individual
                gwf.target_from_template(batch['title'] + '_3_Merge_BAMS_' + individual, merge_bams_new(batch['title'], individual = individual, infiles = BAM_files, outfile = individual+'_merged.bam', input = individual+str(num) + '_sort_dedup.bam')) # denne bruges jo ikke??

                # 4: Filtering the reads        # I think this is about the quality? like - removing bad quality reads?
                gwf.target_from_template(batch['title'] + '_4_filter_bam'+individual, filter_bam_file(batch['title'], individual))

                # 5: Get coverage for each individual
                gwf.target_from_template(batch['title'] + '_5_get_coverage'+individual, get_coverage(batch['title'], individual))

                # 6: Calculate CNV
                gwf.target_from_template(batch['title'] + '_6_calc_cnv' + individual, get_cnv(batch['title'], individual, batch['chromosome']))