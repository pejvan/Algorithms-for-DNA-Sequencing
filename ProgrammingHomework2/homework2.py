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

def boyer_moore_with_counts(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text, p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    num_char_comp = 0
    num_aligments_tried = 0
    
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        num_aligments_tried += 1

        for j in range(len(p)-1, -1, -1):
            num_char_comp += 1
            if p[j] != t[i+j]:    
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break

        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
        
    return occurrences, num_char_comp, num_aligments_tried

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def example_2_1():
    from bm_preproc import BoyerMoore
    p = 'word'
    t = 'there would have been a time for such a word'
    lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
    p_bm = BoyerMoore(p, lowercase_alphabet)
    #print(boyer_moore_with_counts(p, p_bm, t))
    assert boyer_moore_with_counts(p, p_bm, t) == ([40], 15, 12)

def example_2_2():
    from bm_preproc import BoyerMoore
    p = 'needle'
    t = 'needle need noodle needle'
    lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
    p_bm = BoyerMoore(p, lowercase_alphabet)
    #print(boyer_moore_with_counts(p, p_bm, t))
    assert boyer_moore_with_counts(p, p_bm, t) == ([0, 19], 18, 5)

def question1_and_2():
    p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
    reads, qualities = readFastq('chr1.GRCh38.excerpt.fasta')
    assert len(reads) == len(qualities)

    total_char_comp = 0
    total_align_comp = 0

    for t in reads:
        occurrences, num_char_comp, num_aligments_tried = naive_with_counts(p, t)
        total_char_comp += num_char_comp
        total_align_comp += num_aligments_tried

    """How many alignments does the naive exact matching algorithm try when matching the string 
    GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the 
    excerpt of human chromosome 1? (Don't consider reverse complements.)"""
    print 'Question1: ', total_align_comp

    """How many character comparisons does the naive exact matching algorithm try when matching 
    the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) 
    to the excerpt of human chromosome 1? (Don't consider reverse complements.)"""
    print 'Question2: ', total_char_comp
    
def question3():
    p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
    reads, qualities = readFastq('chr1.GRCh38.excerpt.fasta')
    assert len(reads) == len(qualities)

    total_char_comp = 0
    total_align_comp = 0

    from bm_preproc import BoyerMoore
    p_bm = BoyerMoore(p)
    for t in reads:
        occurrences, num_char_comp, num_aligments_tried = boyer_moore_with_counts(p, p_bm, t)
        total_char_comp += num_char_comp
        total_align_comp += num_aligments_tried

    """How many alignments does Boyer-Moore try when matching the string 
    GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) 
    to the excerpt of human chromosome 1? (Don't consider reverse complements.)"""
    print 'Question3: ', total_align_comp
    #print 'Question3--: ', total_char_comp



def main():
    example_1_1()
    example_1_2()
    print "All tests passed successfully for examples in set1"
    example_2_1()
    example_2_2()
    print "All tests passed successfully for examples in set2"
    question1_and_2()
    question3()
   
if __name__ == "__main__":
    main()