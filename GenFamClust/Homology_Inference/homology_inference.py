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
