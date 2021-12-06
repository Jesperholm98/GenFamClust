import sys

# noinspection PyUnresolvedReferences
import Direct_Synteny_score.synteny_score as dsc
#import GenFamClust.Refined_Synteny_score.Refined_Synteny_Score as rsc
#import GenFamClust.Homology_Inference.homology_inference as hi


def main():
    args = sys.argv[1:]

    if len(args) != 2:
        print("programe only takes two input filenames!")
        return 1
    else:
        print("First filename: " + args[0])
        print("Second filename: " + args[1])

    sys_scores = dsc.synteny_score()

    print("successfully completed sys calculation.")


if __name__ == "__main__":
    main()
