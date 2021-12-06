import numpy as np
import sys
import math
import time


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

    return list(set(nchits1).intersection(nchits2))


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


def refined_synteny_score(sys_values, nc_file, beta):

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
            sys_values[(lines[0], lines[1])] = float(lines[2])

    elapsed_time = time.process_time() - t4
    print("Took: ", elapsed_time, "to load sys data into ram memory")

    t3 = time.process_time()
    print("starting to calculate averages for relevant SyS values.")
    averages = calc_average_sys_score_over_ncbeta(nc_scores_over_beta, 'result1.txt')
    elapsed_time = time.process_time() - t3
    print("Finished average SyS values calculation, took: ", elapsed_time, " seconds.\n")
    print("Size of averages is: ", len(averages))

    # print("size of SyS_values: ", sys.getsizeof(SyS_values))
    # print("size of Averages: ", sys.getsizeof(averages))
    # print("size of nc_scores_over_beta: ", sys.getsizeof(nc_scores_over_beta))


    syc_scores = []

    f = open('first100k_np_syc_scores.txt', 'a')
    print("Starting SyC calculation.")
    t5 = time.process_time()
    g = open('missed_pairs_syc.txt', 'a')
    k = open('missing_pairs_h_sys.txt', 'a')
    count = 0
    # count_pairs = 0
    # count_h_pairs = 0
    # not_found_keys = []
    help_dict = {}
    for key in nc_scores_over_beta:

        # if key in SyS_values:
        h = get_nc_hits(nc_scores_over_beta, key[0], key[1], sys_values, help_dict)

        top_l = np.zeros(len(h))
        top_r = np.zeros(len(h))
        bottom_l = np.zeros(len(h))
        bottom_r = np.zeros(len(h))
        idx = 0
        # print("H: ", h, "\n")
        for i in h:
            top_l[idx] = sys_values[(key[0], i)] - averages[key[0]][1] / averages[key[0]][0]
            top_r[idx] = sys_values[(key[1], i)] - averages[key[1]][1] / averages[key[1]][0]
            bottom_l[idx] = (sys_values[(key[0], i)] - averages[key[0]][1] / averages[key[0]][0]) ** 2
            bottom_r[idx] = (sys_values[(key[1], i)] - averages[key[1]][1] / averages[key[1]][0]) ** 2
            idx += 1
        # print("top_l: ", top_l, "\n")
        # print("top_r: ", top_r, "\n")
        # print("bottom_l: ", bottom_l, "\n")
        # print("bottom_r: ", bottom_r, "\n")

        if len(h) != 0:
            result = np.dot(top_l, top_r) / math.sqrt(np.sum(bottom_l) * np.sum(bottom_r))
            f.write(key[0] + "\t" + key[1] + "\t" + str(result) + "\n")  # write to file

        count += 1

        # syc_scores.append()
        # f.write(key[0] + "\t" + key[1] + "\t" + str(val) + "\n")  # write to file

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
    # print("#mis pairs not in SyS:", len(not_found_keys))