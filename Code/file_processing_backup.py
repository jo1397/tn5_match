import pysam
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqIO import read
import time
import gzip
import shutil


ecoli_fasta = os.path.join('/u/c/yehjose1/summer2024', 'ecoli.fasta')
#putting the ecoli genome into a testable variable
for seq_record in SeqIO.parse(ecoli_fasta, "fasta"):
    ecoli_seq = seq_record.seq
    



'''verifying sam file creation'''
# print first few lines of the sam file
# with open(samfile_path, 'r') as f:
#     n = 0
#     for line in f:
#         if n < 20:
#             if line.startswith('@'):
#                 # This is a header line
#                 print(line.strip())
#             else:
#                 # This is an alignment line
#                 fields = line.strip().split('\t')
#                 print(fields)
#             n+=1
#         else:
#             break
  
  
# open sam file      
def plot_error_and_fragment_bias(samfile_path, alignments = 6000000):
    samfile = pysam.AlignmentFile(samfile_path, 'r') # Open the SAM file for reading

    '''get error rate and alignment position biases'''
    total_mismatches = 0
    total_alignments = 0
    start_positions = []
    end_positions = []
    n = 1 #counter for number of alignments


    # Iterate over each alignment in the SAM file
    for alignment in samfile.fetch():
        if n > alignments:
            break
        # Check if 'NM' tag is present
        for tag in alignment.tags:
            if tag[0] == 'NM':         
                # Count the total number of mismatches
                total_mismatches += int(tag[1]) #WHAT IS WRONG HERE????
                break       
        
        # Increment total alignments counter
        total_alignments += 1
        
        #gets paired end positions and fragment length
        if alignment.reference_start is not None and alignment.reference_end is not None:
            if alignment.is_reverse: # Check if alignment is reverse; if so, add to end_positions
                end_positions.append(alignment.reference_start)
                start_positions.append(alignment.reference_end)
            else: # Otherwise, add to start_positions
                start_positions.append(alignment.reference_start)
                end_positions.append(alignment.reference_end)
            
        n += 1

    # Calculate error rate
    if total_alignments > 0:
        error_rate = total_mismatches / (total_alignments*250)
    else:
        error_rate = 0  # handle case where there are no alignments

    print(f"errors alignments: {total_mismatches}")
    print(f"total alignments: {total_alignments}")
    print(f"Error rate: {error_rate}")

    # Calculate alignment position biases
    start_positions = np.array(start_positions)
    end_positions = np.array(end_positions)

    #plot graph of alignment start positions 
    plt.hist(start_positions, bins=200)
    plt.xlabel('Alignment Start Position')
    plt.ylabel('Frequency')
    plt.title('Alignment Start Position Distribution')
    plt.grid(True)
    plt.savefig('optimized_fragmentation_bias_alignment_start_position_distribution_35cycles.png')  # Save the figure as PNG
    plt.close()  # Close the current figure to release memory

    #plot graph of alignment end positions
    plt.hist(end_positions, bins=200)
    plt.xlabel('Alignment End Position')
    plt.ylabel('Frequency')
    plt.title('Alignment End Position Distribution')
    plt.grid(True)
    plt.savefig('optimized_fragmentation_bias_alignment_end_position_distribution_35cycles.png')  # Save the figure as PNG
    plt.close()  # Close the current figure to release memory
    
    print('error_and_fragment_bias done')
    
#plot phred scores distribution
def plot_qual_scores_distribution(samfile_path, alignments=6000000):
    samfile = pysam.AlignmentFile(samfile_path, 'r')  # Open the SAM file for reading

    qual_scores = []  # Dictionary to store Phred score counts
    n = 1  # Counter for number of alignments

    # Iterate over each alignment in the SAM file
    for alignment in samfile.fetch():
        if n > alignments:
            break
        
        # Extract Phred scores if 'QUAL' field is present
        if alignment.query_qualities:
            for score in alignment.query_qualities:
                qual_scores.append(score)

        n += 1

    # Plot Phred scores distribution
    plt.hist(qual_scores, bins=200)
    plt.xlabel('Quality Score')
    plt.ylabel('Frequency')
    plt.title('Quality Scores Distribution')
    plt.grid(True)
    plt.savefig('optmized_Quality_scores_distribution_35cycles.png')
    plt.close()
    
    print('phred_qual_scores_distribution done')

    # print(f"Phred scores distribution: {qual_scores}")
    
    
#plot fragment coverage over the genome
def plot_fragment_coverage(samfile_path, alignments=6000000):
    samfile = pysam.AlignmentFile(samfile_path, 'r')  # Open the SAM file for reading

    # Dictionary to store fragment coverage
    coverage = {}

    n = 1  # Counter for number of alignments

    # Iterate over each alignment in the SAM file
    for alignment in samfile.fetch():
        if n > alignments:
            break

        # Extract coverage and positions
        if alignment.reference_start is not None and alignment.reference_end is not None:
            ref_start = int(alignment.reference_start)
            ref_end = int(alignment.reference_end)

            # Update coverage dictionary
            # Update coverage dictionary
            for i in range(ref_start, ref_end, 10):  # Step by 10 as per your logic
                key = i // 100  # Group positions into bins of size 100
                if key*100 not in coverage:
                    coverage[key*100] = 1
                else:
                    coverage[key*100] += 1

        n += 1

    # Plot fragment coverage as a histogram
    plt.hist(coverage.values(), bins=200, color='blue', alpha=0.7)
    # plt.bar(coverage.keys(), coverage.values(), color='blue', alpha=0.7)
    plt.xlabel('Coverage')
    plt.ylabel('How many bins')
    plt.title('Fragment Coverage')
    plt.grid(True)
    # plt.savefig('optimized_pcr_fragment_coverage_35cycles.png')
    plt.close()
    
    # coverage = dict(sorted(coverage.items(), key=lambda item: item[1]))
    
    # # Write coverage to a text file
    # with open('/u/c/yehjose1/summer2024/optimized_pcr_fragment_coverage_35cycles.txt', 'w') as f:
    #     for key, value in coverage.items():
    #         gc_content = (ecoli_seq[key:key+100].count('G') + ecoli_seq[key:key+100].count('C')) / 100
    #         f.write(f"{gc_content}:\t{coverage[key]}\n")
           
    # Plot fragment coverage vs GC content        
    gc_dict = {}
    num_dict = {}
    for key, value in coverage.items():
        gc_content = (ecoli_seq[key:key+100].count('G') + ecoli_seq[key:key+100].count('C')) / 100
        if gc_content not in gc_dict:
            gc_dict[gc_content] = value
            num_dict[gc_content] = 1
        else:
            gc_dict[gc_content] += value
            num_dict[gc_content] += 1
    
    for key in gc_dict.keys():
        gc_dict[key] = gc_dict[key] / num_dict[key]
            
    plt.scatter(gc_dict.keys(), gc_dict.values(), color='blue', alpha=0.7)
    plt.xlabel('GC Content')
    plt.ylabel('Coverage')
    plt.title('Fragment Coverage vs GC Content')
    plt.grid(True)
    plt.xlim(0, 1)
    plt.savefig('optimized_pcr_fragment_coverage_vs_gc_content_35_cycles.png')
    plt.close()
    
    # Calculate entropy of fragment coverage
    value_sum = sum(coverage.values())
    entropy = 0 
    for key, value in coverage.items():
        p = value / value_sum
        entropy += p * np.log2(p)
    entropy = -entropy

    print('fragment_coverage done')
    # print(sorted(list(coverage.values()), reverse = True))
    print(len(coverage))
    # print(gc_dict)
    print(f"Fragment coverage entropy: {entropy}")
    
    # print(f"Fragment coverage: {coverage}")
    
    
def plot_error_distribution(samfile_path, alignments=6000000):
    samfile = pysam.AlignmentFile(samfile_path, 'r')  # Open the SAM file for reading

    # Dictionary to store error distribution
    error_distribution = {}

    n = 1  # Counter for number of alignments

    # Iterate over each alignment in the SAM file
    for alignment in samfile.fetch():
        if n > alignments:
            break

        # Extract error distribution
        for tag in alignment.tags:
            start = alignment.reference_start
            key = start // 100  # Group positions into bins of size 100
            
            if tag[0] == 'NM':
                error = int(tag[1])
                if key*100 not in error_distribution:
                    error_distribution[key*100] = error
                else:
                    error_distribution[key*100] += error

        n += 1

    # Plot error distribution as a histogram
    error_distribution_list = []
    for pos, err in error_distribution.items():
        error_distribution_list.extend([pos] * err)
    
    plt.hist(error_distribution_list, bins=200, color='red', alpha=0.7)
    plt.xlabel('Genomic Position (x100 bins)')
    plt.ylabel('Number of Errors')
    plt.title('Error Distribution')
    plt.grid(True)
    plt.savefig('optimized_pcr_error_distribution_35cycles.png')
    plt.close()

    print('error_distribution done')
    
    # print(f"Error distribution: {error_distribution}")
    
#convert fasta to fastq because spades only takes fastq files
def fasta_to_fastq(fasta_file, fastq_file):
    with open(fasta_file) as f_in, open(fastq_file, "w") as f_out:
        for record in SeqIO.parse(f_in, "fasta"):
            record.letter_annotations["phred_quality"] = [40] * len(record.seq)
            SeqIO.write(record, f_out, "fastq")
            
#unzip files
def decompress_fasta_gz(input_file, output_file):
    with gzip.open(input_file, 'rt') as f_in:
        with open(output_file, 'wt') as f_out:
            shutil.copyfileobj(f_in, f_out)
            
def decompress_fastq_gz(input_file, output_file):
    with gzip.open(input_file, 'rt') as f_in:
        with open(output_file, 'wt') as f_out:
            shutil.copyfileobj(f_in, f_out)
            
def fix_fastq_format(input_fastq, output_fastq):
    with open(input_fastq, 'r') as infile, open(output_fastq, 'w') as outfile:
        while True:
            header = infile.readline().strip()
            if not header:
                break
            sequence = infile.readline().strip()
            plus = infile.readline().strip()
            quality = infile.readline().strip()

            # Ensure the header line starts with '@'
            if not header.startswith('@'):
                # print(f"Error: Header line does not start with '@'. Attempting to fix.")
                header = '@' + header

            # Ensure the plus line starts with '+'
            if not plus.startswith('+'):
                # print(f"Error: Plus line does not start with '+'. Attempting to fix.")
                plus = '+'

            # Ensure sequence and quality lengths match
            if len(sequence) != len(quality):
                # print(f"Error: Sequence and quality lengths do not match. Attempting to fix.")
                length = len(sequence)
                quality = quality[:length] + 'I' * (length - len(quality))

            outfile.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")
        




if __name__ == '__main__':
    # samfile_path = "/u/c/yehjose1/summer2024/ecoli_aligned.sam"
    # samfile_path = "/u/c/yehjose1/summer2024/SRA_reads/test_SRA_output.sam"
    # samfile_path = "/u/c/yehjose1/summer2024/short_aligned.sam"
    # samfile_path = "/u/c/yehjose1/summer2024/plasmid_aligned.sam"
    # samfile_path = "/u/c/yehjose1/summer2024/pcr_bias_aligned.sam"
    # samfile_path = "/u/c/yehjose1/summer2024/fragmentation_bias_aligned.sam"
    # samfile_path = "/u/c/yehjose1/summer2024/optimized_pcr_aligned_35cycles.sam"
    
    '''plot error rate and alignment position biases'''
    # start = time.time()
    # plot_error_and_fragment_bias(samfile_path)
    # plot_qual_scores_distribution(samfile_path)
    # plot_fragment_coverage(samfile_path)
    # print('time:', time.time() - start)
    # plot_error_distribution(samfile_path)
    # end = time.time()
    # print(f"Time taken: {end-start}")
    
    '''checking tags in sam file'''
    # n = 0
    # for alignment in pysam.AlignmentFile(samfile_path, 'r').fetch():
    #     if n > 1:
    #         break
    #     print(alignment.tags)
    #     n += 1
    
    '''convert fasta to fastq'''
    # start = time.time()
    # fwd_fasta_file = 'D:\Data\summer2024\optimized_pcr_fwd_reads_16cycles.fasta'
    # fwd_fastq_file = 'D:\Data\summer2024\optimized_pcr_fwd_reads_16cycles.fastq'
    # fasta_to_fastq(fwd_fasta_file, fwd_fastq_file)
    # rvs_fasta_file = 'D:\Data\summer2024\optimized_pcr_rvs_reads_16cycles.fasta'
    # rvs_fastq_file = 'D:\Data\summer2024\optimized_pcr_rvs_reads_16cycles.fastq'
    # fasta_to_fastq(rvs_fasta_file, rvs_fastq_file) 
    # end = time.time()
    # print(f"Time taken: {end-start}")
    
    '''unzip files'''
    # decompress_fasta_gz('test_SRA.fasta.gz', 'test_SRA.fasta')
    # decompress_fastq_gz('test_SRA.fastq.gz', 'test_SRA.fastq')
    
    # fix_fastq_format('test_SRA.fastq', 'test_SRA_fixed.fastq')
    
    # fasta_to_fastq('test_SRA.fasta', 'test_SRA_from_fasta.fastq')
    
    '''extract short'''
    # input_file = "ecoli.fasta"
    # output_file = "short.fasta"
    # desired_length = 20000

    # with open(output_file, "w") as out_f:
    #     for record in SeqIO.parse(input_file, "fasta"):
    #         sequence = str(record.seq)[:desired_length]
    #         out_f.write(f">{record.id}\n{sequence}\n")    
    
    '''extract a seq.gz file'''
    # # Define the paths to your input and output files
    # input_filepath = '/u/c/yehjose1/summer2024/gbmam8.seq.gz'
    # output_filepath = '/u/c/yehjose1/summer2024/gbmam8.fasta'

    # with gzip.open(input_filepath, 'rt') as handle:
    #     # Read the sequence data
    #     records = list(SeqIO.parse(handle, 'genbank'))  # Change 'genbank' to the appropriate format if needed

    #     if not records:
    #         raise ValueError("No records found in the input file.")

    #     # Write the sequence data to a FASTA file
    #     with open(output_filepath, 'w') as output_handle:
    #         SeqIO.write(records, output_handle, 'fasta')

    # print(f"Converted {input_filepath} to {output_filepath}")
    
    
    '''extract a seq.gz file to seq object'''
    # # Define the path to your .seq.gz file
    # filepath = '/u/c/yehjose1/summer2024/gbmam8.seq.gz'

    # # Open the gzip file and read the content
    # with gzip.open(filepath, 'rt') as handle:
    #     # Read the sequence data
    #     records = list(SeqIO.parse(handle, 'genbank'))  # Change 'fasta' to the appropriate format if needed

    # # Convert the sequence to a Seq object
    # seq = Seq(str(records.seq))

    # # Print the Seq object
    # print(seq)
    
    '''concatenate fastas'''
    # input_filepath = '/u/c/yehjose1/summer2024/gbmam8.seq.gz'
    # output_filepath = 'concat_gbmam8.fasta'

    # try:
    #     # Open the gzip file and read the content
    #     with gzip.open(input_filepath, 'rt') as handle:
    #         # Read the sequence data
    #         records = list(SeqIO.parse(handle, 'genbank'))

    #         if not records:
    #             raise ValueError("No records found in the input file.")

    #         # Concatenate all sequences
    #         concatenated_sequence = "".join(str(record.seq) for record in records)

    #         # Create a SeqRecord for the concatenated sequence
    #         concatenated_record = SeqRecord(Seq(concatenated_sequence),
    #                                         id="Concatenated_Genome",
    #                                         description="Concatenated genome sequences from GenBank file")

    #         # Write the concatenated sequence to a FASTA file
    #         with open(output_filepath, 'w') as output_handle:
    #             SeqIO.write(concatenated_record, output_handle, 'fasta')

    #     print(f"Converted {input_filepath} to {output_filepath}")

    # except Exception as e:
    #     print(f"An error occurred: {e}")
        
    pass

