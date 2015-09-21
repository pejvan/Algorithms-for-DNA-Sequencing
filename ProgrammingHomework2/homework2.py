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

class Index(object):
    def __init__(self, t, k):
        ''' Create index from all substrings of size 'length' '''
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i+k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer
    
    def query(self, p):
        import bisect
        ''' Return index hits for first k-mer of P '''
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits

def question4():
    """Write a function that, given a length-24 pattern P and given an Index object built on 8-mers, finds all approximate occurrences of P within T with up to 2 mismatches. 
    Insertions and deletions are not allowed. Don't consider any reverse complements.

    How many times does the string GGCGCGGTGGCTCACGCCTGTAAT, which is derived from a human Alu sequence, occur with up to 2 substitutions in the excerpt of human 
    chromosome 1? (Don't consider reverse complements here.)

    - Hint 1: Multiple index hits might direct you to the same match multiple times, but be careful not to count a match more than once.
    - Hint 2: You can check your work by comparing the output of your new function to that of the naive_2mm function implemented in the previous module."""

    p = 'GGCGCGGTGGCTCACGCCTGTAAT'

    #number of mistmatches = 2, so we need to split into 3 (2+1)
    mistmatches_allowed = 2
    num_segments_required = mistmatches_allowed + 1
    k_mer_size = 8
    pattern_size = 24

    reads, qualities = readFastq('chr1.GRCh38.excerpt.fasta')
    consolidated_read = ''.join([read for read in reads])

    index = Index(consolidated_read, k_mer_size)
    p_segments = [ p[i*k_mer_size:i*k_mer_size+k_mer_size] for i in range(num_segments_required) ] # = ['GGCGCGGT', 'GGCTCACG', 'CCTGTAAT']
    #print p_segments

    hits_per_segment = {} # not used
    hit_lists = []
    i = 0
    index_hits = 0
    for segment in p_segments:
        hits = index.query(segment)
        index_hits += 1
        if len(hits) > 0:
            #print segment, hits
            hits_per_segment[segment] = hits
            hit_lists.append(set([hit - i*k_mer_size for hit in hits]))  # we keep the starting point of the pattern in our results set for easy comparison (verification means equality)
                                                                         # and we store them in set for easy equality check by means of intersecting the sets.
        i += 1
    #print hits_per_segment
    #print hit_lists

        
    results = []
    for i in range(len(hit_lists)):
        # WARN: this bit here only works for numbr of mismatches = 2 as we intersect two sets for any values of mistmatches_allowed, it would lead to invalid results otherwise
        intersect = hit_lists[i].intersection(hit_lists[(i+1)%len(hit_lists)])
        if len(intersect) > 0 : 
            for item in intersect:
                results.append(item)
    print len(results)
    return sorted(results), index_hits

# Copy paste from homework1.py
def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        mismatches = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mismatches += 1
                if mismatches > 2:
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def question4_check():
    p = 'GGCGCGGTGGCTCACGCCTGTAAT'
    reads, qualities = readFastq('chr1.GRCh38.excerpt.fasta')
    t = ''.join([read for read in reads])
    return naive_2mm(p, t)

def main():
    example_1_1()
    example_1_2()
    print "All tests passed successfully for examples in set1"
    example_2_1()
    example_2_2()
    print "All tests passed successfully for examples in set2"
    question1_and_2()
    question3()
    print "Question4: How many times does the string GGCGCGGTGGCTCACGCCTGTAAT, which is derived from a human Alu sequence, occur with up to 2 substitutions in the excerpt of human chromosome 1?"
    res, num_index_hits = question4()
    check_res = question4_check()
    assert res == check_res
    print "check with naive_2mm validated"
    print "Question5: how many total index hits are there when searching for occurrences of GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the excerpt of human chromosome 1?"
    print num_index_hits
   
if __name__ == "__main__":
    main()