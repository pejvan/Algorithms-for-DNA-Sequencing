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


def main():
    test01()
    question1()
    example1()
    example2()
    question2()

if __name__ == '__main__':
    main()