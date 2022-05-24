import sys

# noinspection PyUnresolvedReferences
import Direct_Synteny_score.synteny_score as dsc
# noinspection PyUnresolvedReferences
import Refined_Synteny_score.Refined_Synteny_Score as rsc
# noinspection PyUnresolvedReferences
import Homology_Inference.homology_inference as hi


def main():
    args = sys.argv[1:]

    if len(args) != 2:
        print("Programe needs two input filenames!")
        return 1
    else:
        synteny1 = args[0]
        synteny2 = args[1]

    # Parameters
    nc_file = 'nc.txt'
    beta = 0.3
    neighborhood_size = 5

    # SyS Module
    indexing_dict = dsc.synteny_score(synteny1, synteny2, nc_file, neighborhood_size)

    # SyC Module
    rsc.refined_synteny_score(nc_file, beta, indexing_dict)

    # Homology Inference Module
    hi.infer_homology()

    print("Program Done...")


if __name__ == "__main__":
    main()
