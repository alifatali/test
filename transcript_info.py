# Accuireing transcripts sequence inforamtion for Cefra seq dataset
import csv
import os
import urllib
import urllib.request as request
from bs4 import BeautifulSoup
import warnings
import numpy as np

warnings.filterwarnings("ignore")
np.random.seed(1234)

class Gene_Wrapper:
    basedir = os.path.dirname(os.path.abspath(__file__))
    path_to_cefra_cDNA = os.path.join(basedir, 'Data/cefra-seq/cefra_seq_cDNA_screened.fa')
    path_to_cefra_ann = os.path.join(basedir, 'Data/cefra-seq/cefra_seq_cDNA_ann_screened.fa')


    train_test_split_ratio = 0.1

    def __init__(self, id, type, dist):
        self.id = id
        self.type = type
        self.dist = dist
        self.seq = None
        self.ann = None
        np.random.seed(1234)

    @classmethod
    def seq_data_loader(cls, use_ann, dataset, lower_bound=0, upper_bound=np.inf, permute=None):
        """
        permute option is for randanmization test, with three types;
        For conventional data fitting please don't toggle this option on.
        """
        longest = 0
        genes = []
        count = 0
        lengths = [] #Store transctipts lengths
        
        if (dataset == 'cefra-seq'):
            path = cls.path_to_cefra_cDNA
        elif(dataset == 'rnalocate'):
            print("Will process later")
        else:
                        raise RuntimeError('No dataset named {}. Available dataset are "cefra-seq" and "rnalocate"'.format(dataset))


        
        
        with open(path, 'r') as f:
            for line in f:
                if line[0] == '>':
                    tokens = line[1:].split()
                    id = tokens[0]
                    gene_biotype = tokens[1].split(':')[1]
                    transcript_biotype = tokens[2].split(':')[1]
                    dist = [float(c) for c in tokens[3].split(':')[1].split('_')]
                else:
                    lengths.append(len(line[:-1]))
                    if len(line[:-1]) <= lower_bound: # normally not active
                        continue
                    if len(line[:-1]) >= upper_bound:
                        if len(line[:-1]) > longest:
                            longest = len(line[:-1])
                        line = line[:-1][-upper_bound:] # trying to keep (at least) the 3'UTR part
                        # continue
                    if gene_biotype != "protein_coding":
                        # mRNA only
                        continue
                    gene = Gene_Wrapper(id, gene_biotype, dist)
                    gene.seq = line.rstrip().upper()
                    if transcript_biotype != 'protein_coding':
                        count += 1
                    genes.append(gene)
        # print(len(genes))
        # count = 0
        # for gene in genes:
        #     if len(gene.seq) < 4000:
        #         count+=1

        print('non protein coding transcipt gene:', count)
        print('longest sequence: {0}'.format(longest))
        print('average length of mRNAs', sum(lengths)/len(lengths))
    
       
        
        # do some permutations
        genes = np.array(genes)
        genes = genes[np.random.permutation(np.arange(len(genes)))]

        print('Total number of samples:', genes.shape[0])
        if permute:
            print('Warning: permuting mRNA samples!')
            if permute == 1:
                '''preserving the length and nucleotide contents'''
                print('Type 1')
                for gene in genes:
                    gene.seq = ''.join(np.random.permutation(list(gene.seq)))
                    if use_ann:
                        gene.ann = ''.join(np.random.permutation(list(gene.ann)))
            elif permute == 2:
                '''preserving length but altering the actual nucleotide contents'''
                print('Type 2')
                for gene in genes:
                    gene.seq = ''.join(np.random.choice(['A','C','G','T'], len(gene.seq)))
                    if use_ann:
                        gene.ann = ''.join(np.random.choice(['F', 'T', 'I', 'H', 'M', 'S'], len(gene.seq)))
            elif permute == 3:
                '''same length: 3000, and altering the nucleotides content'''
                print('Type 3')
                for gene in genes:
                    gene.seq = ''.join(np.random.choice(['A','C','G','T'], 3000))
                    if use_ann:
                        gene.ann = ''.join(np.random.choice(['F', 'T', 'I', 'H', 'M', 'S'], 3000))
            else:
                raise RuntimeError('Permute option only takes {1,2,3}.')
        
        
        return genes


if __name__ == "__main__":
    
    data = Gene_Wrapper.seq_data_loader(True, 'cefra-seq')
    