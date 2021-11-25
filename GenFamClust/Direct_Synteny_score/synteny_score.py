import sys
import numpy as np
import time
# import multiprocessing as mp
import math


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


def synteny_score(NC_scores='nc.txt', querySyntenyFile1='human_locs.txt', querySyntenyFile2='mouse_locs.txt', k=5,
                  beta=0.5):
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

    # convert into numpy array
    # Synteny1 = np.array(Synteny1)
    # Synteny2 = np.array(Synteny2)

    # load NC file into a dict.
    print("start loading NC_pairs...")
    t1 = time.process_time()
    dict = {}
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
            dict[(lines[0], lines[1])] = float(lines[2])
            inDict += 1

    print("finished loading NC_pairs...")
    elapsed_time = time.process_time() - t1
    print("seconds to load NC-data: ", elapsed_time)
    print("Total elements in file: ", total)
    print("Total elements in DICT: ", inDict, "\n")

    # print("Synteny1: \n")
    # for x in range(10):
    #    print(Synteny1[x])

    # print("\nSynteny2: \n")
    # for x in range(10):
    #    print(Synteny2[x])

    # exit(0)

    SyS_values = {}

    ### calculation without numpy implementation ###
    # ----------------------------------------------

    calc_sys = False

    if calc_sys:
        SyS_size = 0
        # Sys_non_matches = 0

        print("starting SyS calculation...")
        t2 = time.process_time()

        f = open('SyS_scores.txt', 'a')

        for SyntenyIndex1 in range(len(Synteny1)):

            idx_from1 = SyntenyIndex1 - k
            idx_to1 = SyntenyIndex1 + k + 1
            if SyntenyIndex1 - k < 0:
                # Synteny1Neighberhood = Synteny1[0:SyntenyIndex1+k+1]
                idx_from1 = 0
            elif SyntenyIndex1 + k >= len(Synteny1):
                # Synteny1Neighberhood = Synteny1[SyntenyIndex1-k:]
                idx_to1 = len(Synteny1)
            # else:
            #    Synteny1Neighberhood = Synteny1[SyntenyIndex1-k:SyntenyIndex1+k+1]

            for SyntenyIndex2 in range(len(Synteny2)):
                # Synteny2Neighberhood = []

                idx_from2 = SyntenyIndex2 - k
                idx_to2 = SyntenyIndex2 + k + 1
                if SyntenyIndex2 - k < 0:
                    # Synteny2Neighberhood = Synteny2[0:SyntenyIndex2+k+1]
                    idx_from2 = 0

                elif SyntenyIndex2 + k >= len(Synteny2):
                    # Synteny2Neighberhood = Synteny2[SyntenyIndex2-k:]
                    idx_to2 = len(Synteny2)
                # else:
                #    Synteny2Neighberhood = Synteny2[SyntenyIndex2-k:SyntenyIndex2+k+1]

                max_nc = 0
                # for neighbor1 in Synteny1Neighberhood:
                #    for neighbor2 in Synteny2Neighberhood:
                for neighbor1 in calc_neighberhood(Synteny1, SyntenyIndex1, idx_from1,
                                                   idx_to1):  # Synteny1[idx_from1:idx_to1].tolist():
                    for neighbor2 in calc_neighberhood(Synteny2, SyntenyIndex2, idx_from2,
                                                       idx_to2):  # Synteny2[idx_from2:idx_to2].tolist():
                        if (neighbor1, neighbor2) in dict:
                            NC_val = dict[(neighbor1, neighbor2)]
                            if NC_val > max_nc:
                                max_nc = NC_val
                        # else:
                        #    Sys_non_matches += 1

                if max_nc > 0:
                    f.write(Synteny1[SyntenyIndex1][1] + "\t" + Synteny2[SyntenyIndex2][1] + "\t" + str(max_nc) + "\n")
                    # SyS_values[(Synteny1[SyntenyIndex1], Synteny2[SyntenyIndex2])] = max_nc
                    SyS_size += 1

                if SyS_size % 10000000 == 0:
                    print(SyS_size)

        # ----------------------------------------------

        f.close()
        elapsed_time = time.process_time() - t2
        print("Sekunder: ", elapsed_time)
        print("Nr SyS values calculated: ", SyS_size)
        # print("SyS_pairs misses in NC values: ", Sys_non_matches)

        del Synteny1  # free some space
        del Synteny2  # free some space
        del dict

    # exit(0)

    ###
    # Synteny Correlation score calculation
    ###

    # calculate nc_hits_over_beta

    # extract_nc_over_beta()

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

        #print(np.dot(top_l, top_r), " / ", math.sqrt(np.sum(bottom_l) * np.sum(bottom_r)))
        #else:
        #    count_pairs += 1

        count += 1

        #val = 0
        #try:
        #    val = SyS_values[(key[0], key[1])]
        #    h = get_nc_hits(nc_scores_over_beta, key[0], key[1])

        #    count_pairs += 1

        #    top = 0
        #    bottom_l = 0
        #    bottom_r = 0
        #    for el in h:
        #        if (key[0], el) in SyS_values and (key[1], el) in SyS_values:
        #            top += (SyS_values[(key[0], el)] - averages[key[0]][1]) * \
        #                   (SyS_values[(key[1], el)] - averages[key[1]][1])
        #            bottom_l += (SyS_values[(key[0], el)] - averages[key[0]][1]) ** 2
        #            bottom_r += (SyS_values[(key[1], el)] - averages[key[1]][1]) ** 2
        #        else:
        #            count_h_pairs += 1
        #            k.write(key[0] + "\t" + el + "\n")  # write to file

        #    #val = top / math.sqrt(bottom_l * bottom_r)
        #except:
        #    #not_found_keys.append([key[0], key[1]])
        #    g.write(key[0] + "\t" + key[1] + "\t" + str(val) + "\n")  # write to file

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




if __name__ == "__main__":
    main(sys.argv[1:])
