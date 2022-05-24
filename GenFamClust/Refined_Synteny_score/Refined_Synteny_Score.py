import numpy as np
import math
import time
import os.path


def help_get_nchits(nc_hits_over_beta, gene, help_dict):
    """
    this is a help function which adds all gene hits with one gene X to a list and adds it in help_dict for future calculations.
    :param nc_hits_over_beta: filename for file containing nc values over beta
    :param gene: looking for a specific gene id
    :param help_dict: accumulative dictonary to save previously calculated values.
    :return:
    """

    alist = []
    for keys in nc_hits_over_beta:
        if keys[0] == gene and keys[0] != keys[1]:
            alist.append(keys[1])

    help_dict[gene] = alist


def get_nc_hits(nc_over_beta, gene1, gene2, help_nchits):
    """
    this function returns the intersection of nchits of both each genes.
    :param nc_over_beta: filename for file containing nc values over beta
    :param gene1: a gene id
    :param gene2: a gene id
    :param help_nchits: accumulative dictonary to save previously calculated values.
    :return: returns the intersection of two sets, nchits1 and nchits2.
    """
    nchits1 = []
    nchits2 = []

    if gene1 in help_nchits:
        nchits1 = help_nchits[gene1]
    else:
        help_get_nchits(nc_over_beta, gene1, help_nchits)
        nchits1 = help_nchits[gene1]

    if gene2 in help_nchits:
        nchits2 = help_nchits[gene2]
    else:
        help_get_nchits(nc_over_beta, gene2, help_nchits)
        nchits2 = help_nchits[gene2]

    return list(set(nchits1).intersection(set(nchits2)))


def calc_average_sys_score(sys_values):
    """
    this function calculates the average sys score for a certain gene x with all other genes i.  SyS(x, i)
    :return: returns a dictionary with each gene's average sys score.
    """
    averages = {}

    # for all genes in file add their score and the number of times the gene occurs to calculate mean.
    for lines in sys_values:
        if lines[0] in averages:
            count = averages[lines[0]][0] + 1
            total = averages[lines[0]][1] + sys_values[lines]
            averages[lines[0]] = [count, total]
        else:
            count = 1
            total = sys_values[lines]
            averages[lines[0]] = [count, total]
        if lines[1] in averages:
            count = averages[lines[1]][0] + 1
            total = averages[lines[1]][1] + sys_values[lines]
            averages[lines[1]] = [count, total]
        else:
            count = 1
            total = sys_values[lines]
            averages[lines[1]] = [count, total]

    return averages


def extract_nc_over_beta(gene_indexing, nc_file, beta):
    """
    this function writes all pairs of genes with nc score >= beta to a new file
    :param gene_indexing: dictionary to convert string gene identifiers to integer ones.
    :param nc_file: filename of file with all nc scores(default 'nc.txt' that have to be located in /Data/ folder)
    :param beta: a threshhold to filter away unwanted gene pairs with low similarity.
    :return: doesn't return anything.
    """

    f = open('./Data/nc_over_beta.txt', 'a')
    count = 0

    with open('./Data/' + nc_file, 'r') as file:
        for lines in file:
            lines = lines.split()
            if float(lines[2]) >= beta:

                # this try makes sure only genes that can be found in gene_indexing gets used.
                try:
                    f.write(str(gene_indexing[lines[0]]) + "\t" + str(gene_indexing[lines[1]]) + "\t" + lines[2] + "\n")
                    count += 1
                except:
                    # ignore any entries that cannot be found in gene_indexing(gene doesn't exist in either genome1 or genome2)
                    pass
    f.close()
    print("finished extracting nc value over beta=", beta, " #of lines: ", count)


def refined_synteny_score(nc_file, beta, gene_indexing):
    """
    the main funtion for calculating  syc score.
    :param nc_file: the name of the nc file(to read nc values from)
    :param beta: a minimum threshold to include nc values.
    :param gene_indexing: gene indexing to convert string gene identifiers to integers.
    :return: Doesn't return anything.
    """

    sys_values = {}

    sys_vals_exist = os.path.isfile('./Data/sys.txt')

    #make sure file exist and load it to memory.
    if sys_vals_exist:
        print("loading SyS values from file...")
        t = time.process_time()
        with open('./Data/sys.txt', 'r') as file:
            for lines in file:
                lines = lines.split()

                sys_values[(int(lines[0]), int(lines[1]))] = float(lines[2])
        elapsed_time = time.process_time() - t
        print("SyS values loaded. took: ", elapsed_time)
    else:
        print("SyS file not found! Make sure to run synteny_score module or verify the sys.txt file is located in the Data folder.")
        exit(1)


    t2 = time.process_time()
    nc_over_beta_exist = os.path.isfile('./Data/nc_over_beta.txt')

    # if this has been calculated already, then just load that file.
    if nc_over_beta_exist:
        pass
    else:
        print("No nc_over_beta file to load from. Creating that file...")
        extract_nc_over_beta(gene_indexing, nc_file, beta)

    # gene indexing no longer needed.
    gene_indexing = None

    nc_scores_over_beta = {}

    # load nc file into memory
    with open('./Data/nc_over_beta.txt', 'r') as file:
        for lines in file:
            lines = lines.split()
            nc_scores_over_beta[(int(lines[0]), int(lines[1]))] = float(lines[2])

    elapsed_time = time.process_time() - t2
    print("took: ", elapsed_time, "to load nc_scores_over_beta.")

    t3 = time.process_time()
    print("starting to calculate averages for relevant SyS values.")
    averages = calc_average_sys_score(sys_values)
    elapsed_time = time.process_time() - t3
    print("Finished average SyS values calculation, took: ", elapsed_time, " seconds.\n")
    print("Size of averages is: ", len(averages))

    f = open('./Data/syc.txt', 'a')
    print("Starting SyC calculation.")
    t5 = time.process_time()

    help_dict = {}

    # for each pair with nc-value > beta, get a set of genes to include in the calculation
    for count, key in enumerate(nc_scores_over_beta):

        # get the intersection of nc hits for both genes.
        h = get_nc_hits(nc_scores_over_beta, key[0], key[1], help_dict)

        # makes sure the list h isn't empty(otherwise it could divide by 0).
        if len(h) != 0:

            top_l = np.zeros(len(h))
            top_r = np.zeros(len(h))
            bottom_l = np.zeros(len(h))
            bottom_r = np.zeros(len(h))

            # loop over the set of genes to insert them into numpy arrays for the calculation
            for idx, i in enumerate(h):

                try:
                    # try lookup (gene, i) from sys values
                    top_l[idx] = sys_values[(key[0], i)] - averages[key[0]][1] / averages[key[0]][0]
                    top_r[idx] = sys_values[(key[1], i)] - averages[key[1]][1] / averages[key[1]][0]
                    bottom_l[idx] = (sys_values[(key[0], i)] - averages[key[0]][1] / averages[key[0]][0]) ** 2
                    bottom_r[idx] = (sys_values[(key[1], i)] - averages[key[1]][1] / averages[key[1]][0]) ** 2
                except:
                    try:
                        # try lookup (i, gene) from sys values
                        top_l[idx] = sys_values[(i, key[0])] - averages[key[0]][1] / averages[key[0]][0]
                        top_r[idx] = sys_values[(i, key[1])] - averages[key[1]][1] / averages[key[1]][0]
                        bottom_l[idx] = (sys_values[(i, key[0])] - averages[key[0]][1] / averages[key[0]][0]) ** 2
                        bottom_r[idx] = (sys_values[(i, key[1])] - averages[key[1]][1] / averages[key[1]][0]) ** 2

                    except:
                        # Neither (gene, i) OR (i, gene) or (gene, gene) exist in sys values.
                        # Since it doesn't exist just ignore genepair(defaults to 0, so its contribution to the calculation will be 0).
                        pass

            sum1 = np.sum(bottom_l)
            sum2 = np.sum(bottom_r)

            # make sure BOTH sums are > 0 (otherwise it will try divide by 0.)
            if sum1 > 0 and sum2 > 0:
                try:
                    result = np.dot(top_l, top_r) / math.sqrt(sum1 * sum2)
                    if result > 0:
                        f.write(str(key[0]) + "\t" + str(key[1]) + "\t" + str(result) + "\n")  # write to file
                except:
                    # if any errors occur when calculating result
                    # (such as division zero or becomes too close to 0 for the program to handle),
                    # then don't add that entry to the resulting file.
                    pass

    f.close()

    elapsed_time = time.process_time() - t5
    print("SyC calc. complete. took: ", elapsed_time)
