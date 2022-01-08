import random


def main():

    counter1 = -1
    counter2 = -1

    num_lines_human = sum(1 for line in open('human_locs.txt'))
    num_lines_mouse = sum(1 for line in open('mouse_locs.txt'))

    index_human = random.randint(0, num_lines_human-5000)
    index_mouse = random.randint(0, num_lines_mouse-5000)

    print("All human_locs: ", num_lines_human)
    print("All mouse_locs: ", num_lines_mouse)

    print("human index: ", index_human)
    print("mouse index: ", index_mouse)

    f = open("human_locs_5k.txt", "a")
    g = open("mouse_locs_5k.txt", "a")

    count_h_5k = 0
    count_m_5k = 0

    with open('human_locs.txt', 'r') as file:
        for lines in file:
            counter1 += 1

            lines = lines.split()
            if counter1 >= index_human and len(lines) == 3 and count_h_5k < 5000:
                f.write(lines[0] + "\t" + lines[1] + "\t" + lines[2] + "\n")
                count_h_5k += 1
    f.close()

    with open('mouse_locs.txt', 'r') as file:
        for lines in file:
            counter2 += 1

            lines = lines.split()
            if counter2 >= index_mouse and len(lines) == 3 and count_m_5k < 5000:
                g.write(lines[0] + "\t" + lines[1] + "\t" + lines[2] + "\n")
                count_m_5k += 1
    g.close()

    print("new human_5k file has: ", count_h_5k)
    print("new mouse_5k file has: ", count_m_5k)



if __name__ == "__main__":
    main()
