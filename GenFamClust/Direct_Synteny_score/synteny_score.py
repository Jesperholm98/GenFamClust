import sys
import time


def main(argv):
    synteny_score()


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


def synteny_score(querySyntenyFile1='human_locs.txt', querySyntenyFile2='mouse_locs.txt', NC_scores='nc.txt', k=5):

    Synteny1 = []

    with open('./Data/' + querySyntenyFile1, 'r') as file:
        for lines in file:
            lines = lines.split()

            if len(lines) == 3:
                # only keep the lines where there exist values in all three columns

                # Synteny1.append(lines)  #<--- for saving all three columns <chromosone> <order_nr> <gene_id>
                del lines[1]  # remove 2nd column because its not needed for calculation
                Synteny1.append(lines)

    print("SyntenyFile1 finished loading.")
    print("Synteny1 has: ", len(Synteny1), " elements.\n")

    Synteny2 = []

    with open('./Data/' + querySyntenyFile2, 'r') as file:
        for lines in file:
            lines = lines.split()

            if len(lines) == 3:
                # only keep the lines where there exist values in all three columns

                # Synteny2.append(lines)  #<--- for saving all three columns <chromosone> <order_nr> <gene_id>
                del lines[1]  # remove 2nd column because its not needed for calculation
                Synteny2.append(lines)


    print("SyntenyFile2 finished loading.")
    print("Synteny2 has: ", len(Synteny2), " elements.\n")


    # load NC file into a dict.
    print("start loading NC_pairs...")
    t1 = time.process_time()
    nc_dict = {}

    with open('./Data/' + NC_scores, 'r') as file:
        for lines in file:
            lines = lines.split()

            nc_dict[(lines[0], lines[1])] = float(lines[2])

    print("finished loading NC_pairs...")
    elapsed_time = time.process_time() - t1
    print("seconds to load NC-data: ", elapsed_time)

    SyS_values = {}

    #pre-calc neighborhoods
    print("Starting pre-calc for SyS...")
    t = time.process_time()
    synteny1_neighbors, synteny2_neighbors = pre_calc_neighborhoods(Synteny1, Synteny2, k)
    elapsed_time = time.process_time() - t
    print("Finished pre-calc SyS, took: ", elapsed_time)

    print("Starting SyS calculation...")
    f = open('./Data/sys.txt', 'a')
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
                SyS_values[(gene1[1], gene2[1])] = nc_max
    f.close()

    elapsed_time = time.process_time() - t
    print("finished SyS calc, took: ", elapsed_time)

    #print("size(mem) of synteny1_neighbors: ", sys.getsizeof(synteny1_neighbors))
    #print("size(mem) of synteny2_neighbors: ", sys.getsizeof(synteny2_neighbors))
    #print("size(len) of synteny1_neighbors: ", len(synteny1_neighbors))
    #print("size(len) of synteny2_neighbors: ", len(synteny2_neighbors))
    #print("size(length) of SyS_values: ", len(SyS_values))
    #print("size(mem) of synteny2_neighbors: ", sys.getsizeof(SyS_values))

    return SyS_values

if __name__ == "__main__":
    main(sys.argv[1:])
