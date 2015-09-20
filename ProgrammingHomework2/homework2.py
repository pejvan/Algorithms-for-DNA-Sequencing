def naive_with_counts(p, t):
    occurrences = []
    num_char_comp = 0
    num_aligments_tried = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        num_aligments_tried += 1
        for j in range(len(p)):  # loop over characters
            num_char_comp += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, num_char_comp, num_aligments_tried

def example_1_1():
    p = 'word'
    t = 'there would have been a time for such a word'
    #print(naive_with_counts(p, t))
    assert naive_with_counts(p, t) == ([40], 46, 41)

def example_1_2():
    p = 'needle'
    t = 'needle need noodle needle'
    #print(naive_with_counts(p, t))
    assert naive_with_counts(p, t) == ([0, 19], 35, 20)


def main():
    example_1_1()
    example_1_2()
    print "All tests passed successfully for examples in set1"
   
if __name__ == "__main__":
    main()