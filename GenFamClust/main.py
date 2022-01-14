import sys

# noinspection PyUnresolvedReferences
import Direct_Synteny_score.synteny_score as dsc
# noinspection PyUnresolvedReferences
import Refined_Synteny_score.Refined_Synteny_Score as rsc
# noinspection PyUnresolvedReferences
import Homology_Inference.homology_inference as hi


def main():
    args = sys.argv[1:]

    synteny1 = 'human_locs.txt'  # a default value
    synteny2 = 'mouse_locs.txt'  # a default value

    if len(args) != 2:
        print("programe only takes two input filenames!")
        return 1
    else:
        synteny1 = args[0]
        synteny2 = args[1]

    # Parameters
    nc_file = 'nc.txt'
    beta = 0.3
    neighborhood_size = 5

    # SyS Module
    sys_scores = dsc.synteny_score(synteny1, synteny2, nc_file, neighborhood_size)

    # SyC Module
    rsc.refined_synteny_score(sys_scores, nc_file, beta)

    # Homology Inference Module
    hi.infer_homology()

    print("Program Done...")


if __name__ == "__main__":
    main()
