# CSE-549-course-project
# Group Members: Nishanth Muruganandam(110276247), Aadarsh Kenia(109972275), Sridhar Periasami(110104782)

Dataset:
The dataset consists of paired-end reads where the two ends of a read are present in separate .fastq files. The dataset is a pair of fastq files:- Data_1.fastq and Data_2.fastq

Conversion to super-reads:
We use the MaSuRCA assembler to reduce the number of reads by converting them to super-reads. Instructions to compile and run the assember are given at: http://www.genome.umd.edu/docs/MaSuRCA_QuickStartGuide.pdf

Running the Program:
The output after conversion of reads to super-reads is a .fasta file which is the input to our Python program "fastaToSVD.py". 

Run the file using the command: 
python fastaToSVD.py. 

This generates the required output, that is, clusters of super-reads. The output is in the file output.txt.
