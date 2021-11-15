import sys
import numpy as np
import time
import multiprocessing as mp


def main(argv):
    synteny_score()

def get_neighberhood(gene, list, k):
    neighberhood = [gene]

    list_size = len(list)
    idx = 0
    for el in list:
        if el[2] != gene:
            #found correct gene, now get whole neighberhood k.
            idx += 1
        else:
            break

    #print("found gene in list @", idx)
    #print(list[idx])

    tmp1 = idx
    tmp2 = idx
    for i in range(1, k+1):
        #left side neighberhood of gene
        tmp1 -= 1
        if tmp1 < 0:
            break
        neighberhood.append(list[tmp1][2])

    for i in range(1, k+1):
        # right side neighberhood of gene
        tmp2 += 1
        if tmp2 > list_size-1:
            break
        neighberhood.append(list[tmp2][2])

    return neighberhood


def synteny_score(NC_scores='nc.txt', querySyntenyFile1='human_locs.txt', querySyntenyFile2='mouse_locs.txt'):

    #querySyntenyfile
    #gene_order = np.loadtxt(querySyntenyFile1, dtype={'names': ('contig', 'gene'), 'formats': ('U8', 'U6')})

    Synteny1 = []
    wrongF = 0
    totalF = 0

    #load SyntenyFile1
    with open(querySyntenyFile1) as file:
        for lines in file:
            lines = lines.strip().split()
            totalF += 1
            if len(lines) == 3:
                #only keep the lines where there exist values in all three columns
                Synteny1.append(lines)
            else:
                wrongF += 1

    print("SyntenyFile1 finished loading.")
    print("Total nr of lines in First Synteny file: ", totalF)
    print("Nr of bad lines in First Synteny file: ", wrongF)
    print("Synteny1 has: ", len(Synteny1), " elements.\n")


    Synteny2 = []

    wrongS = 0
    totalS = 0

    #load SyntenyFile2
    with open(querySyntenyFile2) as file:
        for lines in file:
            lines = lines.strip().split()
            totalS += 1
            if len(lines) == 3:
                # only keep the lines where there exist values in all three columns
                Synteny2.append(lines)
            else:
                wrongS += 1

    print("SyntenyFile2 finished loading.")
    print("Total nr of lines in Second Synteny file: ", totalS)
    print("Nr of bad lines in Second Synteny file: ", wrongS)
    print("Synteny2 has: ", len(Synteny2), " elements.\n")



    #load NC file into either nparrray OR a dict.
    print("start loading NC_pairs...")
    t1 = time.process_time()
    dict = {}
    #NC_pairs = np.loadtxt(NC_scores, dtype={'names': ('gene1', 'gene2', 'NC_score'), 'formats': ('U6', 'U6', 'f8')})

    inDict = 0
    total = 0
    with open(NC_scores) as file:
        for lines in file:
            lines = lines.strip().split()
            total += 1
            ###
            # the commented away if-statement checks if there already exist (A, B) OR (B, A) in dict.
            # (might be unneccessary to keep all those extra datapoints in memory)
            ###
            # if (lines[0], lines[1]) in dict or (lines[1], lines[0]) in dict:
            #    continue
            if lines[0] != lines[1]:
                dict[(lines[0], lines[1])] = lines[2]
                inDict += 1

    print("finished loading NC_pairs...")
    elapsed_time = time.process_time() - t1
    print("seconds to load NC-data: ", elapsed_time)
    print("Total elements in file: ", total)
    print("Total elements in DICT: ", inDict)

    for i in range(21040, 21064):
        print(Synteny1[i])


    neighberhood = get_neighberhood('Q01534', Synteny1, 5)
    print(neighberhood)


    # calculation



if __name__ == "__main__":
   main(sys.argv[1:])

