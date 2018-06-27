'''
Author: CMK 2018
        
'''

import sys, os
#from os.path import basename # More paths?
from gwf import Workflow
from workflow_templates import *
import json # For pretty printing trees
#import subprocess # For running bash in python
#from datetime import datetime # To set dates in titles
#import re # to escape symbols in strings


# Define the root gwf object
gwf = Workflow()


def pretty_print_dict(dict):
        '''
        Pretty prints a dictionary to get an overview
        '''
        return("dictionary: " + json.dumps(dict, sort_keys=False, indent=4))


def get_file_tree(path):
        '''
        Make a dictionary with files for each individual
        '''

        dict = {}
        with open(path, 'r') as file:
                for line in file:
                        err_id = line.strip()
                        dict[err_id] = [err_id + "_pass_1.fastq", err_id + "_pass_2.fastq", ]
                #print(pretty_print_dict(dict))
                return dict

process = {}
process[0] = False      #init
process[1] = False      #idx ac
process[2] = True
process[3] = False
process[4] = False
process[5] = False
process[6] = False



batches = [{'title': 'batch_1_test',
    'chromosome': 'x',
    'description': 'testing testing',
    
    'rel_batchdir': '../batches/',
    'rel_ac': '../data/ac/ac_chimp_x.fa', # formerly rel_reference
    
    'rel_ind_list': '../data/individuals.txt',
    'rel_fastq_dir': '../data/fastq/',
    }]




for b in batches:
    # initialize
    if process[0]:
        gwf.target_from_template(b['title'] + '_0_init', initialize(b['rel_batchdir'], b['title'], json.dumps(b, sort_keys=False, indent=4)))
    

    # shorthand for working dir for the rest of this loop.
    b['wd'] = b['rel_batchdir'] + b['title'] + '/'
    # ../batches/batch_1_test
    
    #index acs    
    if process[1]:
        gwf.target_from_template(b['title'] + '_1_idx_ac', index_genome(b['wd'], b['rel_ac'])) # I'm not sure if title is even used..

    individuals_tree = get_file_tree(b['rel_ind_list'])
    for individual, fastq_files in individuals_tree.items():
        #print(individual, fastq_files[0], fastq_files[1])

        if process[2]:
            gwf.target_from_template(b['title'] + '_2_map_' + individual, bwa_map_pe(
                b['wd'],
                b['rel_ac'],
                b['rel_fastq_dir'] + fastq_files[0],
                b['rel_fastq_dir'] + fastq_files[1],
                individual))










            



        






# def old_get_individuals(path):
#         '''
#         Takes a path leading to a sample_info.txt-file
#         and returns a list of individuals
#         '''
#         with open(path, 'r') as file:
#                 return [line.split('\t')[0] for line in file][1:]
#                 #                                                ^ We are only interested in the names of the individuals (first column)
#                 #                                                                                         ^ We are only interested in the rows after the labels (labels are on row 0)


# def old_get_filenames(individuals, path):
#         '''
#         Takes a list of individuals
#         and returns a dictionary for pairs of genome-files
#         '''
#         dico_names = {} # dictionary for unique individuals
#         for individual in individuals:  # for each individual in the sample_info.txt-file:
#                 dico_names[individual] = []     # create an empty list for each individual
#                 # file = working_dir+'../Names_files/' + individual + '.txt'    # path to the link-file (every individual has its own file)
#                 with open(path+individual+'.txt') as file:      # close afterwards
#                         lines = [line.strip() for line in file] # transfer the strippedlines to a list, so we can later pair them together. The stripping removes newlines present in the end of each line.
#                         for pair in zip(lines[0::2], lines[1::2]):      # zipping the strided lines, puts them in pairs [(1,2), (2,3), ...]
#                                 dico_names[individual].append(pair)     # append the file pairs to the list in the dictionary for each individual
#         return dico_names


# # This is a list of batches. Each batch is i dictionary with various keys and their values.
# # KEY                                           VALUE
# # title                                         The name of the batch
# # chromosome                            The name of the chromosome being worked on
# # rel_reference                         The relative path to the reference (artificial chromosome)
# # rel_sample_info_file          The relative path to the sample info file, containing a list of individuals
# # rel_individual_genomes        The relative path to the directory containing a file for each individual, each file containing a list og filenames which are the haploid genome builds
# # abs_genome_dir                        The absolute path to the directory containing the files listed in the rel_individual_genomes directory's files.


# # chimp
# batches = [{"title": "batch_x6_chimp",
#         "chromosome": "x",
#         "description": "k\u00f8rsel fredag d 30. apr med ar. chrom. for chimp",
#         "rel_reference": "ac_chimp_x.fa",
#         "rel_sample_info_file": "sample_info_chimp.txt",
#         "rel_individual_genomes": "genomes/",
#         "abs_genome_dir": "/home/cmkobel/MutationRates/faststorage/NEW_PIPELINE/TrimFASTQ/TrimmedFASTQ/"},
        
#         {"title": "batch_x6.3_gorilla",
#         "chromosome": "x",
#         "description": "gorilla med chimps opsætning ændret",
#         "rel_reference": "ac_gorilla_x.fa",
#         "rel_sample_info_file": "sample_info_gorilla.txt",
#         "rel_individual_genomes": "genomes/",
#         "abs_genome_dir": "/home/cmkobel/MutationRates/faststorage/NEW_PIPELINE/TrimFASTQ/TrimmedFASTQ/"},
        
#         {"title": "batch_y6.2_chimp",
#         "chromosome": "y",
#         "description": "ændret opsætning til chimp y",
#         "rel_reference": "ac_chimp_y.fa",
#         "rel_sample_info_file": "sample_info_chimp.txt",
#         "rel_individual_genomes": "genomes/",
#         "abs_genome_dir": "/home/cmkobel/MutationRates/faststorage/NEW_PIPELINE/TrimFASTQ/TrimmedFASTQ/"}]


# # Protocol:
# for batch in batches: # Måske batch bare skal hedde b ?
        
#         # if batch['chromosome'] == 'x':
#         #       batch['rel_reference'] ='artificial_chr_x.fa'
#         # elif batch['chromosome'] == 'y':
#         #       batch['rel_reference'] ='artificial_chr_y.fa'

#         individuals = get_individuals(batch['rel_sample_info_file'])
#         dico_names = get_filenames(individuals, batch['rel_individual_genomes'])

#         # 0: initialize
#         gwf.target_from_template(batch['title']+'_0_initialize', initialize(batch['title'], json.dumps(batch, sort_keys=False, indent=4)))

#         # 1: index
#         gwf.target_from_template(batch['title'] + '_1_index', index_genome(batch['title'], batch['rel_reference']))

#         for individual in dico_names:
#                 BAM_files=[] # for collecting BAM-files later.
#                 for num, pair in enumerate(dico_names[individual]):
#                         # 2: Mapping
#                         gwf.target_from_template(batch['title'] + '_2_map_'+individual+str(num), bwa_map_pe(batch['title'], batch['rel_reference'], batch['abs_genome_dir']+pair[0], batch['abs_genome_dir']+pair[1], individual+str(num)))
#                         BAM_files.append(individual+str(num)+'_sort_dedup.bam') # collect the filenames for use in the merge point
#                 print('BAM_files:', BAM_files)

#                 # 3: Merging the bam files # merge the bam-files for each individual
#                 gwf.target_from_template(batch['title'] + '_3_Merge_BAMS_' + individual, merge_bams_new(batch['title'], individual = individual, infiles = BAM_files, outfile = individual+'_merged.bam', input = individual+str(num) + '_sort_dedup.bam')) # denne bruges jo ikke??

#                 # 4: Filtering the reads        # I think this is about the quality? like - removing bad quality reads?
#                 gwf.target_from_template(batch['title'] + '_4_filter_bam'+individual, filter_bam_file(batch['title'], individual))

#                 # 5: Get coverage for each individual
#                 gwf.target_from_template(batch['title'] + '_5_get_coverage'+individual, get_coverage(batch['title'], individual))

#                 # 6: Calculate CNV
#                 gwf.target_from_template(batch['title'] + '_6_calc_cnv' + individual, get_cnv(batch['title'], individual, batch['chromosome']))