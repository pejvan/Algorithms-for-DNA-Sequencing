def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

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

def naive_with_rc_first(p, t):
    """First, implement a version of the naive exact matching algorithm that is strand-aware. 
    That is, instead of looking only for occurrences of P in T, additionally look for occurrences of the reverse 
    complement of P in T. If P is ACT, your function should find occurrences of both ACT and its reverse complement AGT in T."""
    occurrences = naive(p, t)
    more_occurenences = naive(reverseComplement(p), t)
    return occurrences + more_occurenences

def naive_with_rc_then(p, t):
    """If P and its reverse complement are identical (e.g. AACGTT), then a given match offset 
    should be reported only once. So if your new function is called naive_with_rc, then the old naive 
    function and your new naive_with_rc function should return the same results when P equals its reverse complement."""
    occurrences = naive(p, t)
    revP = reverseComplement(p)
    if p == revP:
        return occurrences
    else:
        more_occurenences = naive(revP, t)
        return occurrences + more_occurenences

def naive_with_rc(p, t):
    return naive_with_rc_then(p, t)

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


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

def example1():
    p = 'CCC'
    ten_as = 'AAAAAAAAAA'
    t = ten_as + 'CCC' + ten_as + 'GGG' + ten_as
    occurrences = naive_with_rc(p, t)
    #print(occurrences)
    assert(occurrences == [10, 23])

def example2():
    p = 'CGCG'
    ten_as = 'AAAAAAAAAA'
    t = ten_as + 'CGCG' + ten_as + 'CGCG' + ten_as
    occurrences = naive_with_rc(p, t)
    #print(occurrences)
    assert(occurrences == [10, 24])

def example3():
    phix_genome = readGenome('phix.fa')
    occurrences = naive_with_rc('ATTA', phix_genome)
    #print(occurrences)
    #print('offset of leftmost occurrence: %d' % min(occurrences))
    assert(min(occurrences) == 62)
    #print('# occurrences: %d' % len(occurrences))
    assert(len(occurrences) == 60)

def question1(genome):
    p = 'AGGT'
    revP = reverseComplement(p)
    assert(revP == 'ACCT')
    occurrences = naive_with_rc(p, genome)
    print 'How many times does AGGT or its reverse complement (ACCT) occur in the lambda virus genome?'
    print len(occurrences)

def question2(genome):
    p = 'TTAA'
    revP = reverseComplement(p)
    assert(revP == 'TTAA')
    occurrences = naive_with_rc(p, genome)
    assert ( len(naive_with_rc_first(p, genome)) == 2 * len(naive_with_rc_then(p, genome)) )
    print 'How many times does TTAA or its reverse complement occur in the lambda virus genome?'
    print 'Hint: TTAA and its reverse complement are equal, so remember not to double count.'
    print len(occurrences)

def question3(genome):
    p = 'ACTAAGT'
    occurrences = naive_with_rc(p, genome)
    print """What is the offset of the leftmost occurrence of ACTAAGT or its reverse complement in the 
    Lambda virus genome? E.g. if the leftmost occurrence of ACTAAGT is at offset 40 (0-based) and the 
    leftmost occurrence of its reverse complement ACTTAGT is at offset 29, then report 29."""
    print min(occurrences)

def question4(genome):
    p = 'AGTCGA'
    revP = reverseComplement(p)
    assert (p != revP)
    occurrences = naive_with_rc(p, genome)
    print """What is the offset of the leftmost occurrence of ACTAAGT or its reverse complement in the 
    Lambda virus genome? E.g. if the leftmost occurrence of ACTAAGT is at offset 40 (0-based) and the 
    leftmost occurrence of its reverse complement ACTTAGT is at offset 29, then report 29."""
    print min(occurrences)

def example10():
    """naive_2mm('ACTTTA', 'ACTTACTTGATAAAGT') should return the list [0, 4]."""
    assert(naive_2mm('ACTTTA', 'ACTTACTTGATAAAGT') == [0, 4])

def example11():
    p = 'CTGT'
    ten_as = 'AAAAAAAAAA'
    t = ten_as + 'CTGT' + ten_as + 'CTTT' + ten_as + 'CGGG' + ten_as
    occurrences = naive_2mm(p, t)
    #print(occurrences)
    assert(occurrences == [10, 24, 38])

def example12():
    phix_genome = readGenome('phix.fa')
    occurrences = naive_2mm('GATTACA', phix_genome)
    #print('offset of leftmost occurrence: %d' % min(occurrences))
    assert(min(occurrences)==10)
    #print('# occurrences: %d' % len(occurrences))
    assert(len(occurrences)==79)

def question5(genome):
    p = 'TTCAAGCC'
    occurrences = naive_2mm(p, genome)
    print """How many times does TTCAAGCC occur in the Lambda virus genome when allowing up to 2 mismatches?"""
    print len(occurrences)
    
def question6(genome):
    p = 'AGGAGGTT'
    occurrences = naive_2mm(p, genome)
    print """What is the offset of the leftmost occurrence of AGGAGGTT in the Lambda virus genome when allowing up to 2 mismatches?"""
    print occurrences[0]

def question7():
    reads, qualities = readFastq('ERR037900_1.first1000.fastq')
    #print qualities, len(qualities)
    """ According to wikipedia: https://en.wikipedia.org/wiki/FASTQ_format"""
    """  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI  """
    """  0........................26...31.......40  """
    threshold = 10
    #print qualities[0], len(qualities[0])
    readLength = len(qualities[0])
    readCount = len(qualities)
    assert(readCount == 1000)

    poorQualityReads=[]
    while len(poorQualityReads) <= 1:
        print threshold
        # for i in range(readLength):
        #     qualitiesForRead=[]
        #     for j in range(readCount):
        #         #print i, j, qualities[j][i], ord(qualities[j][i])-33
        #         quality = ord(qualities[j][i])-33
        #         qualitiesForRead.append(quality)
        #     #print qualitiesForRead
        #     if min(qualitiesForRead) >= threshold:
        #         poorQualityReads.append(i)
        # threshold += -1
        # print threshold

        # for i in range(readCount):
        #     qualitiesForRead=[]
        #     for j in range(readLength):
        #         #print i, j, qualities[j][i], ord(qualities[j][i])-33
        #         quality = ord(qualities[i][j])-33
        #         qualitiesForRead.append(quality)
        #     #print qualitiesForRead
        #     if min(qualitiesForRead) >= threshold:
        #         poorQualityReads.append(j)
        # threshold += -1
        

        for i in range(readCount):
            qualitiesForRead=[ ord(v)-33 for v in qualities[:][i] ]
            if max(qualitiesForRead) <= threshold:
                print max(qualitiesForRead), qualitiesForRead, qualities[:][i]
                poorQualityReads.append(i)
        threshold += 1

    print threshold, poorQualityReads

            

def main():
    example1()
    example2()
    example3()
    print "All tests passed successfully for examples in set1"
    genome = readGenome('lambda_virus.fa')
    question1(genome)
    question2(genome)
    question3(genome)
    question4(genome)
    example10()
    example11()
    example12()
    print "All tests passed successfully for examples in set2"
    question5(genome)
    question6(genome)
    print "Question 7:"
    question7()

if __name__ == "__main__":
    main()