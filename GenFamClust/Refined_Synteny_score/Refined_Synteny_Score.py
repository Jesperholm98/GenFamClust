import numpy as np
import sys
import math
import time
import os.path


def help_get_nchits(nc_hits_over_beta, gene, sys_values, help_dict):
    """
    this is a help function which adds all gene hits with one gene X to a list and adds it in help_dict for future calculations.
    :param nc_hits_over_beta: filename for file containing nc values over beta
    :param gene: looking for a specific gene id
    :param sys_values: dictionary with all sys values.
    :param help_dict: accumulative dictonary to save previously calculated values.
    :return:
    """

    alist = []
    for keys in nc_hits_over_beta:
        if keys[0] == gene:
            if (keys[0], keys[1]) in sys_values:
                alist.append(keys[1])

    help_dict[gene] = alist


def get_nc_hits(nc_over_beta, gene1, gene2, sys_values, help_nchits):
    """
    this function returns the intersection of nchits of both each genes.
    :param nc_over_beta: filename for file containing nc values over beta
    :param gene1: a gene id
    :param gene2: a gene id
    :param sys_values: dictionary with sys values
    :param help_nchits: accumulative dictonary to save previously calculated values.
    :return: returns the intersection of two sets, nchits1 and nchits2.
    """
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


def calc_average_sys_score_over_ncbeta(nc_over_beta):
    """
    this function calculates the average nc score for a certain gene x with all other genes with score >= beta
    :param nc_over_beta: filename that contains nc scores >= beta
    :return: returns a dictionary with each gene's average nc score for all other pairs.
    """
    averages = {}

    # add all genes to the dictionary
    for key in nc_over_beta:
        averages[key[0]] = [0, 0.0]
        averages[key[1]] = [0, 0.0]

    # count = 0
    # for all genes in file add their score and the number of times the gene occurs to calculate mean.
    with open('./Data/sys.txt', 'r') as file:
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


def extract_nc_over_beta(nc_file='nc.txt', beta=0.3):
    """
    this function writes all pairs of genes with nc score >= beta to a new file
    :param nc_file: filename of file with all nc scores(default 'nc.txt' that have to be located in /Data/ folder)
    :param beta: a threshhold to filter away unwanted gene pairs with low similarity.
    :return: doesn't return anything.
    """
    f = open('./Data/nc_over_beta.txt', 'a')
    try:
        count = 0

        with open('./Data/' + nc_file, 'r') as file:
            for lines in file:
                lines = lines.split()
                if float(lines[2]) >= beta:
                    f.write(lines[0] + "\t" + lines[1] + "\t" + lines[2] + "\n")
                    count += 1
        f.close()
        print("finished extracting nc value over beta=", beta, " #of lines: ", count)

    except:
        f.close()
        print("Unable to extract nc_over_beta values to new file.")


def refined_synteny_score(sys_values, nc_file, beta):
    """
    the "main" funtion for calculating  syc score.
    :param sys_values: a dictionary of sys value for genes
    :param nc_file: the name of the nc file(to read nc values from)
    :param beta: a minimum threshold to include nc values.
    :return:
    """

    t2 = time.process_time()
    nc_over_beta_exist = os.path.isfile('./Data/nc_over_beta.txt')

    # if this has been calculated already, then just load that file.
    if nc_over_beta_exist:
        pass
    else:
        print("No nc_over_beta file to load from. Creating that file...")
        extract_nc_over_beta(nc_file, beta)

    nc_scores_over_beta = {}

    # load nc file into memory
    with open('./Data/nc_over_beta.txt', 'r') as file:
        for lines in file:
            lines = lines.split()
            nc_scores_over_beta[(lines[0], lines[1])] = float(lines[2])

    elapsed_time = time.process_time() - t2
    print("took: ", elapsed_time, "to load nc_scores_over_beta.")

    #t4 = time.process_time()
    #print("Loading Sys data into memory...")
    #with open('result1.txt', 'r') as file:
    #    for lines in file:
    #        lines = lines.split()
    #        sys_values[(lines[0], lines[1])] = float(lines[2])
#
    #elapsed_time = time.process_time() - t4
    #print("Took: ", elapsed_time, "to load sys data into ram memory")

    t3 = time.process_time()
    print("starting to calculate averages for relevant SyS values.")
    averages = calc_average_sys_score_over_ncbeta(nc_scores_over_beta)
    elapsed_time = time.process_time() - t3
    print("Finished average SyS values calculation, took: ", elapsed_time, " seconds.\n")
    print("Size of averages is: ", len(averages))

    # print("size of SyS_values: ", sys.getsizeof(SyS_values))
    # print("size of Averages: ", sys.getsizeof(averages))
    # print("size of nc_scores_over_beta: ", sys.getsizeof(nc_scores_over_beta))


    #syc_scores = []

    f = open('./Data/syc.txt', 'a')
    print("Starting SyC calculation.")
    t5 = time.process_time()

    #count = 0
    # count_pairs = 0
    # count_h_pairs = 0
    # not_found_keys = []
    help_dict = {}

    # for each pair with nc value > beta, get a set of genes to include in the calculation
    for count, key in enumerate(nc_scores_over_beta):

        # if key in SyS_values:
        h = get_nc_hits(nc_scores_over_beta, key[0], key[1], sys_values, help_dict)

        top_l = np.zeros(len(h))
        top_r = np.zeros(len(h))
        bottom_l = np.zeros(len(h))
        bottom_r = np.zeros(len(h))
        #idx = 0
        # print("H: ", h, "\n")

        # loop over the set of genes to insert them into numpy arrays for the calculation
        for idx, i in enumerate(h):
            top_l[idx] = sys_values[(key[0], i)] - averages[key[0]][1] / averages[key[0]][0]
            top_r[idx] = sys_values[(key[1], i)] - averages[key[1]][1] / averages[key[1]][0]
            bottom_l[idx] = (sys_values[(key[0], i)] - averages[key[0]][1] / averages[key[0]][0]) ** 2
            bottom_r[idx] = (sys_values[(key[1], i)] - averages[key[1]][1] / averages[key[1]][0]) ** 2
            #idx += 1
        # print("top_l: ", top_l, "\n")
        # print("top_r: ", top_r, "\n")
        # print("bottom_l: ", bottom_l, "\n")
        # print("bottom_r: ", bottom_r, "\n")

        # makes sure the set isn't empty(otherwise it could divide by 0).
        if len(h) != 0:
            try:
                result = np.dot(top_l, top_r) / math.sqrt(np.sum(bottom_l) * np.sum(bottom_r))
                f.write(key[0] + "\t" + key[1] + "\t" + str(result) + "\n")  # write to file
            except:
                # if any errors occur when calculating result
                # (such as division zero or becomes too close to 0 for the program to handle),
                # then don't add that entry to the resulting file.
                pass


        #count += 1

        # syc_scores.append()

        #if count % 500000 == 0:
        #    print(count, "of ", len(nc_scores_over_beta))
        #    print("size of help_dict: ", sys.getsizeof(help_dict))

    f.close()


    elapsed_time = time.process_time() - t5
    print("SyC calc. complete. took: ", elapsed_time)
    #print("#Nr SyC entries:", count)
    #print("#elements in nc_scores_over_beta: ", len(nc_scores_over_beta))
    #print("Total size of help_dict: ", sys.getsizeof(help_dict))
    # print("#mis pairs not in SyS:", len(not_found_keys))