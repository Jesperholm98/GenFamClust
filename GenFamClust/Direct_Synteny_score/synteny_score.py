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


def synteny_score(NC_scores='nc.txt', querySyntenyFile1='human_locs.txt', querySyntenyFile2='mouse_locs.txt', k=5):

    #querySyntenyfile
    #gene_order = np.loadtxt(querySyntenyFile1, dtype={'names': ('contig', 'gene'), 'formats': ('U8', 'U6')})

    ### Syntenyfile1 loading without numpy impl. ###
    # ----------------------------------
    Synteny1 = []
    wrongF = 0
    totalF = 0
    with open(querySyntenyFile1) as file:
        for lines in file:
            lines = lines.strip().split()
            totalF += 1
            if len(lines) == 3:
                #only keep the lines where there exist values in all three columns

                #Synteny1.append(lines)  #<--- for saving all three columns <chromosone> <order_nr> <gene_id>
                Synteny1.append(lines[2])
            else:
                wrongF += 1

    print("SyntenyFile1 finished loading.")
    print("Total nr of lines in First Synteny file: ", totalF)
    print("Nr of bad lines in First Synteny file: ", wrongF)
    print("Synteny1 has: ", len(Synteny1), " elements.\n")
    # ----------------------------------


    ### Syntenyfile2 loading without numpy impl. ###
    #----------------------------------
    Synteny2 = []

    wrongS = 0
    totalS = 0

    with open(querySyntenyFile2) as file:
        for lines in file:
            lines = lines.strip().split()
            totalS += 1
            if len(lines) == 3:
                # only keep the lines where there exist values in all three columns

                #Synteny2.append(lines)  #<--- for saving all three columns <chromosone> <order_nr> <gene_id>
                Synteny2.append(lines[2])

            else:
                wrongS += 1

    print("SyntenyFile2 finished loading.")
    print("Total nr of lines in Second Synteny file: ", totalS)
    print("Nr of bad lines in Second Synteny file: ", wrongS)
    print("Synteny2 has: ", len(Synteny2), " elements.\n")
    # ----------------------------------

    # convert into numpy array
    #Synteny1 = np.asarray(Synteny1)
    #Synteny2 = np.asarray(Synteny2)


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


    #print("Synteny1: \n")
    #for x in range(10):
    #    print(Synteny1[x])

    #print("\nSynteny2: \n")
    #for x in range(10):
    #    print(Synteny2[x])


    #exit(0)



    SyS_values = {}
    SyS_size = 0
    Sys_non_matches = 0

    print("starting SyS calculation...")
    t2 = time.process_time()

    ### calculation without numpy implementation (8 mins)###
    #----------------------------------------------

    for SyntenyIndex1 in range(len(Synteny1)):
        Synteny1Neighberhood = []

        idx_from1 = SyntenyIndex1-k
        idx_to1 = SyntenyIndex1+k+1
        if SyntenyIndex1-k < 0:
            #Synteny1Neighberhood = Synteny1[0:SyntenyIndex1+k+1]
            idx_from1 = 0
        elif SyntenyIndex1+k >= len(Synteny1):
            #Synteny1Neighberhood = Synteny1[SyntenyIndex1-k:]
            idx_to1 = len(Synteny1)
        #else:
        #    Synteny1Neighberhood = Synteny1[SyntenyIndex1-k:SyntenyIndex1+k+1]

        for SyntenyIndex2 in range(len(Synteny2)):
            Synteny2Neighberhood = []

            idx_from2 = SyntenyIndex2-k
            idx_to2 = SyntenyIndex2+k+1
            if SyntenyIndex2 - k < 0:
                #Synteny2Neighberhood = Synteny2[0:SyntenyIndex2+k+1]
                idx_from2 = 0

            elif SyntenyIndex2 + k >= len(Synteny2):
                #Synteny2Neighberhood = Synteny2[SyntenyIndex2-k:]
                idx_to2 = len(Synteny2)
            #else:
            #    Synteny2Neighberhood = Synteny2[SyntenyIndex2-k:SyntenyIndex2+k+1]

            max_nc = 0
            #for neighbor1 in Synteny1Neighberhood:
            #    for neighbor2 in Synteny2Neighberhood:
            for neighbor1 in Synteny1[idx_from1:idx_to1]:
                for neighbor2 in Synteny2[idx_from2:idx_to2]:
                    if (neighbor1, neighbor2) in dict:
                        NC_val = float(dict[(neighbor1, neighbor2)])
                        if NC_val > max_nc:
                            max_nc = NC_val
                    #else:
                    #    Sys_non_matches += 1


            if max_nc > 0:
                SyS_values[(Synteny1[SyntenyIndex1], Synteny2[SyntenyIndex2])] = max_nc
                SyS_size += 1

            #if SyS_size % 1000 == 0:
            #    print(SyS_size, " : ", Sys_non_matches)

    # ----------------------------------------------


    elapsed_time = time.process_time() - t2
    print("Sekunder: ", elapsed_time)
    print("Nr SyS values calculated: ", SyS_size)
    print("SyS_pairs misses in NC values: ", Sys_non_matches)


if __name__ == "__main__":
   main(sys.argv[1:])

