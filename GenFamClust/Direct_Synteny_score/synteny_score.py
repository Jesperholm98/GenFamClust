import time
import os.path


def calc_neighberhood(alist, gene_idx, idx_from, idx_to):
    """
    this is a help function to pre_calc_neighborhoods
    :param alist: a synteny list.
    :param gene_idx: index of current gene to get neighborhood for
    :param idx_from: start from which index to include genes inside the neighborhood
    :param idx_to: final index to include genes for neighborhood
    :return: returns a list with whole neighborhood for a specific gene.
    """
    neighberhood = []

    try:
        for idx in range(idx_from, idx_to):
            if alist[idx][0] == alist[gene_idx][0]:
                # make sure only genes from same chromosome are added
                neighberhood.append(int(alist[idx][1]))
    except:
        print("Some error occured in calc_neighberhood.")

    return neighberhood


def pre_calc_neighborhoods(synteny1, synteny2, k):
    """
    This function calculates the neighborhood of size k for all genes in synteny1 and synteny2.
    :param synteny1: synteny list for first genome
    :param synteny2: synteny list for second genome
    :param k: size of neighborhood(#genes up and down from a certain gene in respective synteny file)
    :return: the neighborhoods for all genes in both synteny files
    """
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


def data_preprocessing(querySyntenyFile1, querySyntenyFile2):
    """
    This function makes sure no duplicates genes exist within a genome, and also gives indexes the remaining genes for future usage.
    :param querySyntenyFile1: input Synteny file for first genome.
    :param querySyntenyFile2: input Synteny file for second genome.
    :return: returns a dictionary with an index for a certain gene to be able to convert string gene identifier to an integer for efficiency(both memory and speed).
    """

    indexing = {}

    t = time.process_time()
    f = open('./Data/genome1_genes.txt', 'w')
    g = open('./Data/genome2_genes.txt', 'w')

    genome1 = set()
    genome2 = set()
    count = 1

    # indexing first genome
    with open('./Data/' + querySyntenyFile1, 'r') as file:
        for lines in file:
            lines = lines.split()
            if len(lines) == 3:
                if lines[2] in genome1:
                    pass
                else:
                    f.write(lines[0] + "\t" + str(count) + "\t" + lines[2] + "\n")
                    indexing[lines[2]] = count
                    count += 1
                    genome1.add(lines[2])
    f.close()

    # indexing second genome
    with open('./Data/' + querySyntenyFile2, 'r') as file:
        for lines in file:
            lines = lines.split()
            if len(lines) == 3:
                if lines[2] in genome2:
                    pass
                else:
                    g.write(lines[0] + "\t" + str(count) + "\t" + lines[2] + "\n")
                    indexing[lines[2]] = count
                    count += 1
                    genome2.add(lines[2])
    g.close()
    elapsed_time = time.process_time() - t
    print("time to preprocess query files: ", elapsed_time, " seconds.")
    return indexing


def synteny_score(querySyntenyFile1, querySyntenyFile2, NC_scores, k):
    """
    This is the main function for calculating SyS scores.
    :param querySyntenyFile1: synteny file for first genome
    :param querySyntenyFile2: synteny file for reference genome
    :param NC_scores: dictionary consisting of nc scores for genes
    :param k: size of neighborhood for a gene
    :return: returns a dictionary to convert string gene identifier to respective integer gene identifier.
    """
    sys_file_exist = os.path.isfile('./Data/sys.txt')

    indexing_dict = data_preprocessing(querySyntenyFile1, querySyntenyFile2)

    # if this has been calculated already, then just load that file.
    if sys_file_exist:
        # sys file already exist.
        print("sys file already exist. Either move or delete the existing sys.txt in th Data folder. Moving on to next module.")
        return indexing_dict
    else:
        # no sys file exist so create and calculate sys values and write them to file.
        Synteny1 = []

        # load a list of genes from first genome
        with open('./Data/genome1_genes.txt', 'r') as file:
            for lines in file:
                lines = lines.split()
                del lines[2]
                lines[1] = int(lines[1])
                Synteny1.append(lines)

        print("SyntenyFile1 finished loading.")
        print("Synteny1 has: ", len(Synteny1), " genes.\n")

        Synteny2 = []

        # load a list of genes from first genome
        with open('./Data/genome2_genes.txt', 'r') as file:
            for lines in file:
                lines = lines.split()
                del lines[2]
                lines[1] = int(lines[1])
                Synteny2.append(lines)

        print("SyntenyFile2 finished loading.")
        print("Synteny2 has: ", len(Synteny2), " genes.\n")

        # load NC file into a dict.
        print("start loading NC_pairs...")
        t1 = time.process_time()
        nc_dict = {}
        with open('./Data/' + NC_scores, 'r') as file:
            for lines in file:
                lines = lines.split()

                try:
                    nc_dict[(indexing_dict[lines[0]], indexing_dict[lines[1]])] = float(lines[2])
                except:
                    # data mismatch between nc and rest of program(most likely genes in nc not found in synteny) so ignore those entries.
                    pass

        print("finished loading # of NC_pairs: ", len(nc_dict))

        elapsed_time = time.process_time() - t1
        print("seconds to load NC-data: ", elapsed_time)

        # pre-calc neighborhoods to reduce runtime
        print("Starting pre-calc for SyS...")
        t = time.process_time()
        synteny1_neighbors, synteny2_neighbors = pre_calc_neighborhoods(Synteny1, Synteny2, k)

        elapsed_time = time.process_time() - t
        print("Finished pre-calc SyS, took: ", elapsed_time)

        print("Starting SyS calculation...")
        f = open('./Data/sys.txt', 'a')
        t1 = time.process_time()

        # calculate genome1 * genome2 genepairs
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
                    f.write(str(gene1[1]) + "\t" + str(gene2[1]) + "\t" + str(nc_max) + "\n")
                    #SyS_values[(gene1[1], gene2[1])] = nc_max

        elapsed_time = time.process_time() - t1
        print("genepairs synteny1*synteny2 took: ", elapsed_time)

        t2 = time.process_time()
        # calculate genome1 * genome1 genepairs.
        for idx, el in enumerate(Synteny1):
            neighborhood1 = synteny1_neighbors[el[1]]
            for idxx in range(idx + 1, len(Synteny1)):
                neighborhood2 = synteny1_neighbors[Synteny1[idxx][1]]

                nc_max = 0

                for a in neighborhood1:
                    for b in neighborhood2:
                        if (a, b) in nc_dict:
                            nc_max = max(nc_dict[(a, b)], nc_max)

                if nc_max > 0:
                    f.write(str(el[1]) + "\t" + str(Synteny1[idxx][1]) + "\t" + str(nc_max) + "\n")
                    #SyS_values[(el[1], Synteny1[idxx][1])] = nc_max
        elapsed_time = time.process_time() - t2
        print("genepairs synteny1*synteny1 took: ", elapsed_time)

        t3 = time.process_time()
        # calculate genome2 * genome2 genepairs
        for idx, el in enumerate(Synteny2):
            neighborhood1 = synteny2_neighbors[el[1]]
            for idxx in range(idx + 1, len(Synteny2)):
                neighborhood2 = synteny2_neighbors[Synteny2[idxx][1]]

                nc_max = 0

                for a in neighborhood1:
                    for b in neighborhood2:
                        if (a, b) in nc_dict:
                            nc_max = max(nc_dict[(a, b)], nc_max)

                if nc_max > 0:
                    f.write(str(el[1]) + "\t" + str(Synteny2[idxx][1]) + "\t" + str(nc_max) + "\n")
                    #SyS_values[(el[1], Synteny2[idxx][1])] = nc_max
        f.close()
        elapsed_time = time.process_time() - t3
        print("genepairs synteny2*synteny2 took: ", elapsed_time)

        elapsed_time = time.process_time() - t
        print("SyS calculation complete. Whole module took: ", elapsed_time)

    return indexing_dict
