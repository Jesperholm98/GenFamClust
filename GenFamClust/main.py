import sys

def main():
    args = sys.argv[1:]

    if len(args) != 2:
        print("programe only takes two input filenames!")
        return 1
    else:
        print("First filename: " + args[0])
        print("Second filename: " + args[1])


if __name__ == "__main__":
    main()