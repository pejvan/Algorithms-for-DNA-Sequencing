def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists, return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def scs(ss):
    """ Returns shortest common superstring of given strings, which must be the same length """
    import itertools
    shortest_sup = None
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest

def scs_list(ss):
    """ Returns the alphebeticaly sorted list of shortest common superstrings of given strings, which must be the same length """
    import itertools
    shortest_sup = []
    shortest_length = 0
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_length == 0:
            shortest_sup.append(sup)
            shortest_length = len(sup)
        elif len(sup) < shortest_length:
            shortest_length = len(sup)
            shortest_sup = [sup]
        elif len(sup) == shortest_length:
            shortest_sup.append(sup)
        else:
            #simply ignore and move on
            None
    return sorted(shortest_sup)


def test01():
    """Consider the input strings ABC, BCA, CAB. One shortest common superstring is ABCAB but another is BCABC and another is CABCA."""
    inputStrings = ['ABC', 'BCA', 'CAB']
    shortestCommonSuperstring = scs(inputStrings)
    assert shortestCommonSuperstring in ['ABCAB', 'BCABC', 'CABCA']

def question1():
    """What is the length of the shortest common superstring of the following strings? CCT, CTT, TGC, TGG, GAT, ATT"""
    inputStrings = ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT']
    shortestCommonSuperstring = scs(inputStrings) #== 'CCTTGGATTGC'

    #print("shortest common superstring", shortestCommonSuperstring)
    print("What is the length of the shortest common superstring of the following strings? CCT, CTT, TGC, TGG, GAT, ATT")
    print(len(shortestCommonSuperstring))

def example1():
    strings = ['ABC', 'BCA', 'CAB']
    # Returns just one shortest superstring
    assert scs(strings) == 'ABCAB'
    # Returns list of all superstrings that are tied for shorest
    shortestList = scs_list(strings)
    assert shortestList == ['ABCAB', 'BCABC', 'CABCA'], 'found ' + str(shortestList)

def example2():
    ##from pprint import pprint
    strings = ['GAT', 'TAG', 'TCG', 'TGC', 'AAT', 'ATA']
    # Returns just one shortest superstring
    assert scs(strings) == 'TCGATGCAATAG'
    # Returns list of all superstrings that are tied for shorest
    shortestList = scs_list(strings)
    ##pprint(shortestList)

    assert shortestList ==  ['AATAGATCGTGC',
                             'AATAGATGCTCG',
                             'AATAGTCGATGC',
                             'AATCGATAGTGC',
                             'AATGCTCGATAG',
                             'TCGAATAGATGC',
                             'TCGATAGAATGC',
                             'TCGATGCAATAG',
                             'TGCAATAGATCG',
                             'TGCAATCGATAG'], 'found ' + str(shortestList)



def question2():
    """How many different shortest common superstrings are there for the input strings given in the previous question?
    Hint 1: You can modify the scs function to keep track of this."""
    #from pprint import pprint
    shortestList = scs_list( ['CCT', 'CTT', 'TGC', 'TGG', 'GAT', 'ATT'] )
    #pprint(shortestList)
    print("How many different shortest common superstrings are there for the input strings given in the previous question?")
    print(len(shortestList))


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
            #All the reads are the same length (100 bases) 
            assert len(seq) == 100
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

overlap_cache = {} # 

#copied from: http://nbviewer.ipython.org/github/Benlangmead/ads1-notebooks/blob/master/4.02_GreedySCS.ipynb
def pick_maximal_overlap(reads, k):
    """ Return a pair of reads from the list with a
        maximal suffix/prefix overlap >= k.  Returns
        overlap length 0 if there are no such overlaps."""
    import itertools
    reada, readb = None, None
    best_olen = 0

    for a, b in itertools.permutations(reads, 2):
        
        if (a,b) in overlap_cache.keys() and overlap_cache[(a,b)] >= k:
            overlap_len = overlap_cache[(a,b)]
            #print("overlap_cache hit for (a,b):", (a,b), "min was", k, "result was:", overlap_len, "cache size:", len(overlap_cache))
            return a, b, overlap_len
        else:
            olen = overlap(a, b, min_length=k)
            #olen = overlap(a, b)
            if olen > best_olen:
                reada, readb = a, b
                best_olen = olen
                overlap_cache[(a,b)]=best_olen

    return reada, readb, best_olen

#copied from: http://nbviewer.ipython.org/github/Benlangmead/ads1-notebooks/blob/master/4.02_GreedySCS.ipynb
def greedy_scs(reads, k):
    """ Greedy shortest-common-superstring merge.
        Repeat until no edges (overlaps of length >= k)
        remain. """
    lenCacheBefore = len(overlap_cache)
    read_a, read_b, olen = pick_maximal_overlap(reads, k)

    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        
        read_a, read_b, olen = pick_maximal_overlap(reads, k)
    
    lenCacheAfter = len(overlap_cache)
    #print("Cached {0} items during this pass".format(lenCacheAfter-lenCacheBefore))

    return ''.join(reads)

def validated_greedy_scs():
    res1 = greedy_scs(['ABC', 'BCA', 'CAB'], 2)
    #print('res1: ', res1)
    assert res1 == 'CABCA'

    res2 = greedy_scs(['ABCD', 'CDBC', 'BCDA'], 1)
    #print('res2: ', res2)
    assert res2 == 'CDBCABCDA'

def question3and4():
    from datetime import datetime
    
    reads, qualities = readFastq('ads1_week4_reads.fq')
    #print(len(reads))
    
    for i in range (30, 100):
        print("timestamp: ", datetime.now())
        result = greedy_scs(list(reads), i) #we make a copy of the reads as the greedy_scs modifies the list
        print("Found result which is {0} bases long for k={1}".format( len(result), i ) )
        
        # Hint: the virus genome you are assembling is exactly 15,894 bases long
        #assert len(result) == 15894
        if len(result) == 15894:
            print("Question3: ", result.count('A'))
            print("Question4: ", result.count('T'))
            return

    for i in range (30, 1, -1):
        print("timestamp: ", datetime.now())
        result = greedy_scs(list(reads), i) #we make a copy of the reads as the greedy_scs modifies the list
        print("Found result which is {0} bases long for k={1}".format( len(result), i ) )
        
        # Hint: the virus genome you are assembling is exactly 15,894 bases long
        #assert len(result) == 15894
        if len(result) == 15894:
            print("Question3: ", result.count('A'))
            print("Question4: ", result.count('T'))
            return

def main():
    test01()
    question1()
    example1()
    example2()
    question2()
    validated_greedy_scs()

    question3and4()
    

if __name__ == '__main__':
    main()