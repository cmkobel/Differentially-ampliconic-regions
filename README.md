# Differentially-ampliconic-regions





### Overview on the pipeline so far:
--

The ten PRJEB1357 genomes are used to increase the size of the previous projectA
https://www.ncbi.nlm.nih.gov/sra/?term=PRJEB1357
The genomes were downloaded with wget from the ncbi ftp 

--

sra archives were unpacked with: fastq-dump --outdir . --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ../sra/ERR2250*.sra > logjun24.txt
See https://edwards.sdsu.edu/research/fastq-dump/ for details on the parameters selected
