import Bio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio import SeqIO
from Bio.SeqIO import write
import subprocess
from subprocess import run
import numpy as np
from numpy import random
import random
from scipy.stats import skewnorm #for fitting to curve
from matplotlib import pyplot as plt #plotting
import time
import multiprocessing
from multiprocessing import Pool, cpu_count, Manager
import os

def open_fasta(fwd_fasta, rvs_fasta, verbose = False):
    '''
    Opens a fasta file and returns a list of SeqRecords
    '''
    fwd_reads = []
    rvs_reads = []
    reads = []
    for record in SeqIO.parse(fwd_fasta, 'fasta'):
        fwd_reads.append(record.seq)
        reads.append(record.seq)
        
    i = 0
    for record in SeqIO.parse(rvs_fasta, 'fasta'):
        rvs_reads.append(record.seq)
        reads[i] = reads[i] + 'NNNNNNNNN' + record.seq
        i += 1
        
        
    if verbose:
        print('num reads:', len(fwd_reads), len(rvs_reads))
        # print(fwd_reads[0], '\n',  rvs_reads[0], '\n', reads[0])
    
    return fwd_reads, rvs_reads, reads


def create_dicts(contigs, verbose = False):
    '''
    Creates a dictionary of reads with the first 9 bases as the key for fwd reads and last 9 bases as the key for rvs reads
    '''
    fwd_dict = {}
    rvs_dict = {}
    
    for read in contigs:
        if read[:9] in fwd_dict:
            fwd_dict[read[:9]].append(read)
        else:
            fwd_dict[read[:9]] = [read]
            
        if read[-9:] in rvs_dict:
            rvs_dict[read[-9:]].append(read)
        else:
            rvs_dict[read[-9:]] = [read]
            
    if verbose:
        print('num fwd keys:', len(fwd_dict), 'num rvs keys:', len(rvs_dict))
    
    return fwd_dict, rvs_dict

def save_dict(fwd_dict, rvs_dict, fwd_output_file, rvs_output_file):
    '''
    Saves the dictionaries to a file
    '''
    start = time.time()
    with open(fwd_output_file, 'w') as output_handle:
        for key in fwd_dict:
            output_handle.write(str(key) + '\n')
            for seq in fwd_dict[key]:
                output_handle.write(str(seq) + '\n')
                
    with open(rvs_output_file, 'w') as output_handle:           
        for key in rvs_dict:
            output_handle.write(str(key) + '\n')
            for seq in rvs_dict[key]:
                output_handle.write(str(seq) + '\n')
    
    print('saved', 'time = ', time.time()-start)
    
    
def match_nines(fwd_dict, rvs_dict): #need to check logic; still unfinished function
    '''
    Matches the 9bp of the reads that may overlap
    '''
    scaffold = []
    for key in fwd_dict:
        for fwd_contig in fwd_dict[key]:
            matched = False
            if key in rvs_dict:
                if rvs_dict[key]:
                    rvs_contig = rvs_dict[key].pop(0)
                    scaffold.append(rvs_contig + fwd_contig[9:])
                    fwd_dict[key].remove(fwd_contig)
                    matched = True
                    
                
            elif matched == False and key.reverse_complement() in fwd_dict:
                if fwd_dict[key.reverse_complement()]:
                    scaffold.append(fwd_dict[key.reverse_complement()][0].reverse_complement() + fwd_contig[9:])
                    fwd_dict[key.reverse_complement()].pop(0)
                    fwd_dict[key].remove(fwd_contig)
                    matched = True
                    
            elif matched == False:
                scaffold.append(fwd_contig)
                fwd_dict[key].remove(fwd_contig)
                
    for key in rvs_dict:
        for rvs_contig in rvs_dict[key]:
              matched = False        
              if key.reverse_complement() in rvs_dict:
                  if rvs_dict[key.reverse_complement()]:
                      scaffold.append(rvs_contig + rvs_dict[key.reverse_complement()][0].reverse_complement()[9:])
                      rvs_dict[key.reverse_complement()].pop(0)
                      rvs_dict[key].remove(rvs_contig)
                      matched = True
                      
              elif matched == False:
                  scaffold.append(rvs_contig)
                  rvs_dict[key].remove(rvs_contig)
                  
    return scaffold

def assembler(reads, verbose = False):
    '''
    Assembles the reads into a scaffold
    '''
    shortened = True
    iter = 0
    start_iter = time.time()
    while shortened == True:
        fwd_dict, rvs_dict = create_dicts(reads)
        scaffold = match_nines(fwd_dict, rvs_dict)
        
        if len(scaffold) < len(reads):
            shortened = True
            reads = scaffold
            if verbose:
                print('iter:', iter, 'num reads:', len(reads), 'time:', time.time()-start_iter)
                iter += 1
                start_iter = time.time()
            
        else:
            shortened = False
            if verbose:
                print('iter:', iter, 'num reads:', len(reads), 'time:', time.time()-start_iter)
    
    if verbose:
        print('total iterations:', iter, 'final num reads:', len(scaffold))
    
    return scaffold

def save_fasta(scaffold, output_file):
    '''
    Saves the scaffold to a fasta file
    '''
    start = time.time()
    with open(output_file, 'w') as output_handle:
        for i, seq in enumerate(scaffold):
            record = SeqRecord(Seq(seq), id = 'scaffold' + str(i), description = 'contig length ' + str(len(seq)))
            write(record, output_handle, 'fasta')
    
    print('saved to', output_file, 'time = ', time.time()-start)
                
            
            


if __name__ == '__main__':
    fwd_fasta = '/home/jyeh/summer2024/genomes_200_fwd_reads.fasta'
    rvs_fasta = '/home/jyeh/summer2024/genomes_200_rvs_reads.fasta'
    output_file = '/home/jyeh/summer2024/assembler_dict_scaffold_200_genomes.fasta'
    fwd_dict_file = '/home/jyeh/summer2024/genomes_200_fwd_dict.txt'
    rvs_dict_file = '/home/jyeh/summer2024/genomes_200_rvs_dict.txt'
    
    start = time.time()
    
    fwd_reads, rvs_reads, reads = open_fasta(fwd_fasta, rvs_fasta, verbose = True)
    print(reads[0])
    fwd_dict, rvs_dict = create_dicts(reads, verbose = True)
    save_dict(fwd_dict, rvs_dict, fwd_dict_file, rvs_dict_file)
    # print(fwd_dict)
    
    scaffold = assembler(reads, verbose = True)
    save_fasta(scaffold, output_file)
    
    
    end = time.time()
    print('time:', end-start)