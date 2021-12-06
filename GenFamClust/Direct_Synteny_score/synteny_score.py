import sys
import numpy as np
import time
# import multiprocessing as mp
import math
import os


def main(argv):
    synteny_score()


def calc_average_sys_score_over_ncbeta(nc_over_beta, sys_values):
    averages = {}

    for key in nc_over_beta:
        averages[key[0]] = [0, 0.0]
        averages[key[1]] = [0, 0.0]

    # count = 0
    with open(sys_values, 'r') as file:
        for lines in file:
            # count += 1
            lines = lines.split()
            if lines[0] in averages:
                count = averages[lines[0]][0] + 1
                total = averages[lines[0]][1] + float(lines[2])
                averages[lines[0]] = [count, total]
            if lines[1] in averages:
                count = averages[lines[1]][0] + 1
                total = averages[lines[1]][1] + float(lines[2])
                averages[lines[1]] = [count, total]

            # if count % 1000000 == 0:
            #    print(count)

    return averages


def extract_nc_over_beta(nc_file='nc.txt', beta=0.5):
    f = open('nc_over_beta.txt', 'a')
    try:
        count = 0

        with open(nc_file) as file:
            for lines in file:
                lines = lines.split()
                if float(lines[2]) > beta:
                    f.write(lines[0] + "\t" + lines[1] + "\t" + lines[2] + "\n")
                    count += 1
        f.close()
        print("finished extracting nc value over beta=", beta, " #of lines: ", count)

    except:
        f.close()
        print("Unable to extract nc_over_beta values to new file.")


def help_get_nchits(nc_hits_over_beta, gene, sys_values, help_dict):

    alist = []
    for keys in nc_hits_over_beta:
        if keys[0] == gene:
            if (keys[0], keys[1]) in sys_values:
                alist.append(keys[1])

    help_dict[gene] = alist

def get_nc_hits(nc_over_beta, gene1, gene2, sys_values, help_nchits):
    nchits1 = []
    nchits2 = []

    if gene1 in help_nchits:
        nchits1 = help_nchits[gene1]
    else:
        help_get_nchits(nc_over_beta, gene1, sys_values, help_nchits)
        nchits1 = help_nchits[gene1]

    if gene2 in help_nchits:
        nchits2 = help_nchits[gene2]
    else:
        help_get_nchits(nc_over_beta, gene2, sys_values, help_nchits)
        nchits2 = help_nchits[gene2]

    # for el in nc_over_beta:  # (gene1, gene2) = 1.11  #(gene1, el[1]) OR (gene2, el[1])
    #    if gene1 == el[0]:
    #        if (gene1, el[1]) in sys_values:
    #            nchits1.append(el[1])
    #    if gene2 == el[0]:
    #        if (gene2, el[1]) in sys_values:
    #            nchits2.append(el[1])
    #
    # return list(set(nchits1).intersection(nchits2))

    return list(set(nchits1).intersection(nchits2))


def calc_neighberhood(alist, gene_idx, idx_from, idx_to):
    neighberhood = []

    try:
        for idx in range(idx_from, idx_to):
            if alist[idx][0] == alist[gene_idx][0]:
                neighberhood.append(alist[idx][1])
    except:
        print("Some error occured in calc_neighberhood.")

    return neighberhood


def pre_calc_neighborhoods(synteny1, synteny2, k):
    gene_neighborhoods1 = {}
    gene_neighborhoods2 = {}

    for idx, gene in enumerate(synteny1):
        idx_from1 = max(idx - k, 0)
        idx_to1 = min(idx + k + 1, len(synteny1))

        neighborhood = calc_neighberhood(synteny1, idx, idx_from1, idx_to1)
        gene_neighborhoods1[gene[1]] = neighborhood

    for idx, gene in enumerate(synteny2):
        idx_from2 = max(idx - k, 0)
        idx_to2 = min(idx + k + 1, len(synteny2))

        neighborhood = calc_neighberhood(synteny2, idx, idx_from2, idx_to2)
        gene_neighborhoods2[gene[1]] = neighborhood

    return gene_neighborhoods1, gene_neighborhoods2


def synteny_score(NC_scores='nc.txt', querySyntenyFile1='human_locs.txt', querySyntenyFile2='mouse_locs.txt', k=5,
                  beta=0.3):
    # querySyntenyfile
    # gene_order = np.loadtxt(querySyntenyFile1, dtype={'names': ('contig', 'gene'), 'formats': ('U8', 'U6')})

    ### Syntenyfile1 loading without numpy impl. ###
    # ----------------------------------
    Synteny1 = []
    wrongF = 0
    totalF = 0
    with open(querySyntenyFile1) as file:
        for lines in file:
            lines = lines.split()
            totalF += 1
            if len(lines) == 3:
                # only keep the lines where there exist values in all three columns

                # Synteny1.append(lines)  #<--- for saving all three columns <chromosone> <order_nr> <gene_id>
                del lines[1]  # remove 2nd column because its not needed for calculation
                Synteny1.append(lines)
            else:
                wrongF += 1

    print("SyntenyFile1 finished loading.")
    print("Total nr of lines in First Synteny file: ", totalF)
    print("Nr of bad lines in First Synteny file: ", wrongF)
    print("Synteny1 has: ", len(Synteny1), " elements.\n")
    # ----------------------------------

    ### Syntenyfile2 loading without numpy impl. ###
    # ----------------------------------
    Synteny2 = []

    wrongS = 0
    totalS = 0

    with open(querySyntenyFile2) as file:
        for lines in file:
            lines = lines.split()
            totalS += 1
            if len(lines) == 3:
                # only keep the lines where there exist values in all three columns

                # Synteny2.append(lines)  #<--- for saving all three columns <chromosone> <order_nr> <gene_id>
                del lines[1]  # remove 2nd column because its not needed for calculation
                Synteny2.append(lines)

            else:
                wrongS += 1

    print("SyntenyFile2 finished loading.")
    print("Total nr of lines in Second Synteny file: ", totalS)
    print("Nr of bad lines in Second Synteny file: ", wrongS)
    print("Synteny2 has: ", len(Synteny2), " elements.\n")
    # ----------------------------------


    # load NC file into a dict.
    print("start loading NC_pairs...")
    t1 = time.process_time()
    nc_dict = {}
    # NC_pairs = np.loadtxt(NC_scores, dtype={'names': ('gene1', 'gene2', 'NC_score'), 'formats': ('U6', 'U6', 'f8')})

    inDict = 0
    total = 0
    with open(NC_scores) as file:
        for lines in file:
            lines = lines.split()
            total += 1
            ###
            # the commented away if-statement checks if there already exist (A, B) OR (B, A) in dict.
            # (might be uneccessary to keep all those extra datapoints in memory)
            ###
            # if (lines[0], lines[1]) in dict or (lines[1], lines[0]) in dict:
            #    continue
            # if lines[0] != lines[1]:
            nc_dict[(lines[0], lines[1])] = float(lines[2])
            inDict += 1

    print("finished loading NC_pairs...")
    elapsed_time = time.process_time() - t1
    print("seconds to load NC-data: ", elapsed_time)
    print("Total elements in file: ", total)
    print("Total elements in DICT: ", inDict, "\n")

    SyS_values = {}

    print("Starting pre-calc for SyS...")

    #pre-calc neighborhoods
    t = time.process_time()
    synteny1_neighbors, synteny2_neighbors = pre_calc_neighborhoods(Synteny1, Synteny2, k)
    elapsed_time = time.process_time() - t
    print("Finished pre-calc SyS, took: ", elapsed_time)

    print("Starting SyS calculation...")
    f = open('./../Data/test_sys.txt', 'a')
    t = time.process_time()
    for gene1 in Synteny1:
        for gene2 in Synteny2:

            neighborhood1 = synteny1_neighbors[gene1[1]]
            neighborhood2 = synteny2_neighbors[gene2[1]]
            nc_max = 0

            for a in neighborhood1:
                for b in neighborhood2:
                    if (a, b) in nc_dict:
                        nc_max = max(nc_dict[(a, b)], nc_max)

            if nc_max > 0:
                f.write(gene1[1] + "\t" + gene2[1] + "\t" + str(nc_max) + "\n")
                #SyS_values[(gene1[1], gene2[1])] = nc_max
    f.close()

    elapsed_time = time.process_time() - t
    print("finished SyS calc, took: ", elapsed_time)

    #print("size(mem) of synteny1_neighbors: ", sys.getsizeof(synteny1_neighbors))
    #print("size(mem) of synteny2_neighbors: ", sys.getsizeof(synteny2_neighbors))
    #print("size(len) of synteny1_neighbors: ", len(synteny1_neighbors))
    #print("size(len) of synteny2_neighbors: ", len(synteny2_neighbors))
    #print("size(length) of SyS_values: ", len(SyS_values))
    #print("size(mem) of synteny2_neighbors: ", sys.getsizeof(SyS_values))

    # exit(0)

    ###
    # Synteny Correlation score calculation
    ###

    # calculate nc_hits_over_beta

    extract_nc_over_beta()

    nc_scores_over_beta = {}

    with open('nc_over_beta.txt', 'r') as file:
        for lines in file:
            lines = lines.split()
            nc_scores_over_beta[(lines[0], lines[1])] = float(lines[2])

    t4 = time.process_time()
    print("Loading Sys data into memory...")
    with open('result1.txt', 'r') as file:
        for lines in file:
            lines = lines.split()
            SyS_values[(lines[0], lines[1])] = float(lines[2])

    elapsed_time = time.process_time() - t4
    print("Took: ", elapsed_time, "to load sys data into ram memory")

    t3 = time.process_time()
    print("starting to calculate averages for relevant SyS values.")
    averages = calc_average_sys_score_over_ncbeta(nc_scores_over_beta, 'result1.txt')
    elapsed_time = time.process_time() - t3
    print("Finished average SyS values calculation, took: ", elapsed_time, " seconds.\n")
    print("Size of averages is: ", len(averages))

    #print("size of SyS_values: ", sys.getsizeof(SyS_values))
    #print("size of Averages: ", sys.getsizeof(averages))
    #print("size of nc_scores_over_beta: ", sys.getsizeof(nc_scores_over_beta))


    syc_scores = []

    f = open('first100k_np_syc_scores.txt', 'a')
    print("Starting SyC calculation.")
    t5 = time.process_time()
    g = open('missed_pairs_syc.txt', 'a')
    k = open('missing_pairs_h_sys.txt', 'a')
    count = 0
    #count_pairs = 0
    #count_h_pairs = 0
    #not_found_keys = []
    help_dict = {}
    for key in nc_scores_over_beta:

        #if key in SyS_values:
        h = get_nc_hits(nc_scores_over_beta, key[0], key[1], SyS_values, help_dict)

        top_l = np.zeros(len(h))
        top_r = np.zeros(len(h))
        bottom_l = np.zeros(len(h))
        bottom_r = np.zeros(len(h))
        idx = 0
        #print("H: ", h, "\n")
        for i in h:
            top_l[idx] = SyS_values[(key[0], i)] - averages[key[0]][1] / averages[key[0]][0]
            top_r[idx] = SyS_values[(key[1], i)] - averages[key[1]][1] / averages[key[1]][0]
            bottom_l[idx] = (SyS_values[(key[0], i)] - averages[key[0]][1] / averages[key[0]][0]) ** 2
            bottom_r[idx] = (SyS_values[(key[1], i)] - averages[key[1]][1] / averages[key[1]][0]) ** 2
            idx += 1
        #print("top_l: ", top_l, "\n")
        #print("top_r: ", top_r, "\n")
        #print("bottom_l: ", bottom_l, "\n")
        #print("bottom_r: ", bottom_r, "\n")

        if len(h) != 0:
            result = np.dot(top_l, top_r) / math.sqrt(np.sum(bottom_l) * np.sum(bottom_r))
            f.write(key[0] + "\t" + key[1] + "\t" + str(result) + "\n")  # write to file


        count += 1

        # syc_scores.append()
        #f.write(key[0] + "\t" + key[1] + "\t" + str(val) + "\n")  # write to file

        if count % 500000 == 0:
            print(count, "of ", len(nc_scores_over_beta))
            print("size of help_dict: ", sys.getsizeof(help_dict))

    f.close()
    g.close()
    k.close()

    elapsed_time = time.process_time() - t5
    print("SyC calc. complete. took: ", elapsed_time)
    print("#Nr SyC entries:", count)
    print("Total size of help_dict: ", sys.getsizeof(help_dict))
    #print("#mis pairs not in SyS:", len(not_found_keys))

def infer_homology(syc_values, nc_values):
    #load nc_values into memory
    nc_dict = {}
    with open(nc_values, 'r') as file:
        for lines in file:
            lines = lines.split()
            nc_dict[(lines[0], lines[1])] = float(lines[2])

    f = open('homolog_genes.txt', 'a')
    print("Starting to infer homology")
    with open(syc_values, 'r'):
        for lines in file:
            lines = lines.split()
            if (nc_dict[(lines[0], lines[1])]**2 + 0.25 * float(lines[2])**2 - 0.25) > 0:
                f.write(lines[0] + "\t" + lines[1] + str(nc_dict[(lines[0], lines[1])]) + "\t" + str(float(lines[2])) + "\n")
    f.close()
    print("Done.")


if __name__ == "__main__":
    main(sys.argv[1:])
