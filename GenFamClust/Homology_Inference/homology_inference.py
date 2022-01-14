import time


def infer_homology(syc_file='syc.txt', nc_file='nc.txt'):
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
