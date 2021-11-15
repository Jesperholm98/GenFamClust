import sys
import numpy as np
import time



def main(argv):
    synteny_score()

def synteny_score(NC_scores='nc.txt', querySyntenyFile='test.txt'):

    #querySyntenyfile
    gene_order = np.loadtxt(querySyntenyFile, dtype={'names': ('contig', 'gene'), 'formats': ('U8', 'U6')})

    #load NC file into either nparrray OR a dict.
    print("start loading NC_pairs...")
    t = time.process_time()
    NC_pairs = np.loadtxt(NC_scores, dtype={'names': ('gene1', 'gene2', 'NC_score'), 'formats': ('U6', 'U6', 'f8')})
    end = time.time()
    print("finished loading NC_pairs...")
    elapsed_time = time.process_time() - t
    print("seconds to load NC-data: ", elapsed_time)

    one = NC_pairs[1]
    gene1 = one[0]
    gene2 = one[1]
    NC_Value = one[2]

    print("One line: ", one)
    print("type of gene1: ", type(gene1))
    print("type of gene2: ", type(gene2))
    print("type of NC_Value: ", type(NC_Value))

    # calculation



if __name__ == "__main__":
   main(sys.argv[1:])

