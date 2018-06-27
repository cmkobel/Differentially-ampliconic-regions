# Author: CMK


import os

# 0: initialize
def initialize(batchdir, title, description):
    '''
    make a dir for each individual and possibly paste the description of it
    '''
    inputs = []
    outputs = [] # dirs are not outputs, so no need to fill anything here
    options = {'cores': 1, 'memory': 100, 'walltime': '00:01:00'}
    spec = '''
        mkdir {batchdir}{title}
        cd {batchdir}{title}
        echo $(date)'\n{description}\n' > batch_parameters.txt'''.format(
            batchdir = batchdir,
            title = title,
            description = description)
    return inputs, outputs, options, spec   


# 1: Index the ref. genome
def index_genome(wd, refgenome):
    #refgenome_stem = os.path.splitext(refgenome)[0]
    refgenome_stem = os.path.splitext(refgenome)[0].split('/')[-1]
    print(refgenome_stem) #debug
    
    inputs = [refgenome]
    outputs = [wd + '/' + refgenome_stem + extension for extension in ['.amb', '.ann', '.pac', '.sa', '.bwt', '.dict', '.fai']]
    options = {'cores': 1, 'memory': 2000, 'walltime': '00:25:00'}
    #options = {}
    spec =  '''sleep 1;
        cd {wd}; cp ../{refgenome} .;
        source /com/extra/bwa/0.7.5a/load.sh;
        bwa index -p {refgenome_stem} -a bwtsw {refgenome_stem}.fa
        '''.format(
        wd = wd,
        refgenome_stem = refgenome_stem,
        refgenome = refgenome)
    #spec =  """mkdir {title}; cd {title}; touch test; """.format(title=title)

    return inputs, outputs, options, spec



#2
def bwa_map_pe(batch_wd, refgenome, read1, read2, individual):
    '''
    Maps the hap. genomes to the reference

    1st sambamba view -F : Filters for the paired reads
    2nd sambamba view -F : specify the reference genome for writing the output
    sambamba sort : sorts the reads per coordinate (default)
    sambamba markdup : removes the duplicate reads in bam file
    sambamba flagstat : gives a flag stat to BAM files
    sambamba view -F : filter the bam files
    '''
    refgenome_stem = os.path.splitext(refgenome)[0].split('/')[-1]

    inputs = [read1,
              read2,
              refgenome]

    for extension in ['.amb', '.ann', '.pac']: inputs.append(batch_wd + '/' + refgenome_stem + extension) # kan godt flattenes pænere
    # ooutputs = [
    #   title+'/'+individual+'_sorted.bam',
    #   title+'/'+individual+'_sorted.bam.bai',
    #   title+'/'+individual+'_unsorted.bam',
    #   title+'/'+individual+'_sort_dedup.bam',
    #   title+'/'+individual+'_sort_dedup.bam.bai',
    #   title+'/'+individual+'_sort_dedup.bam.flagstat',
    #   title+'/'+individual+'.COMPLETED']

    outputs = [batch_wd + '/' + individual+extension for extension in [
        '_sorted.bam',
        '_sorted.bam.bai',
        '_unsorted.bam',
        '_sort_dedup.bam',
        '_sort_dedup.bam.bai',
        '_sort_dedup.bam.flagstat',
        '.COMPLETED']]

    options = {'cores': 4, 'memory': 2000, 'walltime': '00:01:00'}
    #options = {'cores': 16, 'memory': 24000, 'walltime': '24:00:00'}

    spec = '''
        source /com/extra/bwa/0.7.5a/load.sh
        source /com/extra/sambamba/0.5.1/load.sh
        cd {batch_wd}
        mkdir {individual}
        #cd {individual}

        echo ../{read1}
        echo ../{refgenome}
        echo {refgenome_stem}


        bwa mem -M -t 16 -a {refgenome_stem} ../{read1} ../{read2} \
        | sambamba view -f bam -F "proper_pair" -S -t 8 /dev/stdin \
        | sambamba view -f bam -T {refgenome_stem}.fa /dev/stdin > {individual}/{individual}_unsorted.bam
          
        # sambamba sort -t 16 -m 24GB --out={individual}_sorted.bam --tmpdir=/scratch/$GWF_JOBID/ {individual}_unsorted.bam
        # sleep 10

        #rm -f {individual}_unsorted.bam
        # sambamba markdup -t 16 {individual}_sorted.bam --tmpdir=/scratch/$GWF_JOBID/ {individual}_sort_dedup.bam
        # sleep 10

        #rm -f {individual}_sorted.bam
        #rm -f {individual}_sorted.bam.bai
        # sambamba flagstat -t 16 {individual}_sort_dedup.bam > {individual}_sort_dedup.bam.flagstat
        sleep 10'''.format(
            #title = title,
            batch_wd = batch_wd, 
            refgenome = refgenome,
            refgenome_stem = refgenome_stem,
            read1 = read1,
            read2 = read2,
            individual = individual)

    return inputs, outputs, options, spec


# 3
# Jeg tror ikke dette trin er nødvendigt, da der kun er en enkelt 
def merge_bams_new(title, individual, infiles, outfile, input):
    inputs = [title+'/'+i for i in infiles]
    outputs = [title+'/'+individual+'/'+outfile] # must be a list
    options = {'cores': 4, 'memory': 8000, 'walltime': '02:45:00'}
    spec = '''
        source /com/extra/sambamba/0.5.1/load.sh
        cd {title}

        mkdir {individual}

        #culprit:
        sambamba merge -t 4 {individual}/{outfile} {inbams}


        sleep 30
        #rm -f {inbams}
        echo "All done at "$(date) > {individual}.COMPLETED'''.format(
            title = title,
            outfile  = outfile,
            inbams = ' '.join(infiles),
            individual = individual)

    #print(spec)
    return inputs, outputs, options, spec


# 4
def filter_bam_file(title, individual):
    '''
    Filter bam files

    This functions filters the BAM file so we only get the reads that map uniquely and are paired
    '''
    inputs = [title+'/'+individual+'/'+individual+'_merged.bam'] # er det nødvendigt at kalde den individual to gange?
    outputs = [title+'/'+individual+'/'+individual+'_filtered.bam',
               title+'/'+individual+'/'+individual+'_filtered.bam.bai']
    options = {'cores': 1, 'memory': 8000, 'walltime': '02:00:00'}
    spec = '''
        cd {title}
        source /com/extra/sambamba/0.5.1/load.sh
        sambamba view -F "not (duplicate or secondary_alignment or unmapped) and mapping_quality >= 50 and cigar =~ /100M/ and [NM] < 3" -f bam ./{ind}/{ind}_merged.bam > ./{ind}/{ind}_filtered.bam
        sambamba index ./{ind}/{ind}_filtered.bam 
        #rm -f {ind}_merged.bam'''.format(title=title, ind=individual)
    #print(spec)
    return inputs, outputs, options, spec


# 5
def get_coverage(title, individual):
    """
    saves the files *_cov.txt inside each individuals directory
    This file contains the number of reads touching each position in the ar.chromosome.
    """
    inputs = [title+'/'+individual+'/'+individual+'_filtered.bam']
    outputs = [title+'/'+individual+'/'+individual+'_cov.txt']
    options = {'cores': 4, 'memory': 4000, 'walltime': '04:00:00'}
    spec = '''
        cd {title}
        source /com/extra/sambamba/0.5.1/load.sh
        source /com/extra/samtools/1.3/load.sh
        cd {ind}

        samtools depth {ind}_filtered.bam > {ind}_cov.txt'''.format(title=title, ind=individual)
    return inputs, outputs, options, spec


# 6:
def get_cnv(title, individual, chrom):
    inputs = [title + '/' + individual + '/' + individual + '_cov.txt']
    #ooutputs = [inputs[0]+'_pd_median.txt']
    outputs = [title + '/cn_medians/' + chrom + '_' + individual + '_cn_median.csv']
    options = {'cores': 4, 'memory': 64000, 'walltime': '04:00:00'}
    spec = '''
    if [ ! -d {dir} ]; then mkdir {dir}; fi
    ~/miniconda3/bin/python compute_cov_median.py {input} {output}
    '''.format(dir = title+'/cn_medians/', input = inputs[0], output = outputs[0])
    return inputs, outputs, options, spec
    #fejl ved manglende fil? fordi den ligge en aterproccessing.