import time


def infer_homology(syc_file='syc.txt', nc_file='nc.txt'):
    """
    This is the "main" function for infering homology. The result gets written to a result file: "homolog_genes.txt"
    :param syc_file: filename for file containing syc values
    :param nc_file: filename for file containing nc values
    :return: no return value.
    """
    t = time.process_time()

    index_to_gene = {}
    gene_to_index = {}
    gene_set1 = set()
    gene_set2 = set()

    # retrieve two dictionarys to convert string gene id -> gene integer id and the opposite way.
    with open('./Data/genome1_genes.txt', 'r') as file:
        for lines in file:
            lines = lines.split()

            index_to_gene[int(lines[1])] = lines[2]
            gene_to_index[lines[2]] = int(lines[1])
            gene_set1.add(int(lines[1]))
    with open('./Data/genome2_genes.txt', 'r') as file:
        for lines in file:
            lines = lines.split()

            index_to_gene[int(lines[1])] = lines[2]
            gene_to_index[lines[2]] = int(lines[1])
            gene_set2.add(int(lines[1]))

    nc_dict = {}
    # load nc_values into memory
    with open('./Data/' + nc_file, 'r') as file:
        for lines in file:
            lines = lines.split()
            try:
                # only load genepairs from nc_scores that are from genome1 or genome2
                nc_dict[(gene_to_index[lines[0]], gene_to_index[lines[1]])] = float(lines[2])
            except:
                pass

    f = open('./Data/homolog_genes.txt', 'w')
    sucess_counter = 0
    print("Starting to infer homology")
    with open('./Data/' + syc_file, 'r') as file:
        for lines in file:
            lines = lines.split()
            try:
                # makes sure genepair (a, b) have both syc & nc scores and a and b are from different genomes
                if (nc_dict[(int(lines[0]), int(lines[1]))]**2 + 0.25 * (float(lines[2]))**2 - 0.25) > 0 and int(lines[0]) in gene_set1 and int(lines[1]) in gene_set2 or int(lines[0]) in gene_set2 and int(lines[1]) in gene_set1:
                    f.write(index_to_gene[int(lines[0])] + "\t" + index_to_gene[int(lines[1])] + "\t" + str(nc_dict[(int(lines[0]), int(lines[1]))]) + "\t" + lines[2] + "\n")
                    sucess_counter += 1
            except:
                # discarded genepairs gets skipped either because they don't satisfy the above conditions(no NC or SyC scores or genepair are not from different genomes)
                pass
    f.close()
    elapsed_time = time.process_time() - t
    print("Homology Inference done... Took: ", elapsed_time)
    print("#homolog genes: ", sucess_counter)
