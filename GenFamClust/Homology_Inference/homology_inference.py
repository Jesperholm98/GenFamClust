import time


def infer_homology(syc_file='syc.txt', nc_file='nc.txt'):
    """
    This is the "main" function for infering homology. The result gets written to a result file: "homolog_genes.txt"
    :param syc_file: filename for file containing syc values
    :param nc_file: filename for file containing nc values
    :return: no return value.
    """
    t = time.process_time()

    nc_dict = {}
    # load nc_values into memory
    with open('./Data/' + nc_file, 'r') as file:
        for lines in file:
            lines = lines.split()
            nc_dict[(lines[0], lines[1])] = float(lines[2])

    f = open('./Data/homolog_genes.txt', 'a')

    print("Starting to infer homology")
    with open('./Data/' + syc_file, 'r') as file:
        for lines in file:
            lines = lines.split()
            if (nc_dict[(lines[0], lines[1])]**2 + 0.25 * float(lines[2])**2 - 0.25) > 0:
                f.write(lines[0] + "\t" + lines[1] + "\t" + str(nc_dict[(lines[0], lines[1])]) + "\t" + str(float(lines[2])) + "\n")
    f.close()
    elapsed_time = time.process_time() - t
    print("Done... Took: ", elapsed_time)
