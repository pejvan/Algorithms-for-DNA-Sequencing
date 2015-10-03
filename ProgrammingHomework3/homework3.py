def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def editDistance(x, y):
    """Returns the edit distance between two strings, x and y"""
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]

def bestApproximateMatchEditDistance(p, t):
    """Returns the edit distance between two strings, p and t"""
    # Create distance matrix
    D = []
    for i in range(len(p)+1):
        D.append([0]*(len(t)+1))
    
    # Initialize first row and column of matrix
    for i in range(len(p)+1):
        D[i][0] = i
    # See slide 4 on  0440_approx__editdist3.pdf
    # First row is already initialised to zero so we simply just comment the following two lines.
    #for i in range(len(p)+1):
    #    D[0][i] = i
    
    # Fill in the rest of the matrix
    for i in range(1, len(p)+1):
        for j in range(1, len(t)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if p[i-1] == t[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)

    # Best Approximate Match Distance is the smallest value of the last row
    return min(D[-1])

def test1():
    """P = GCGTATGC within T = TATTGGCTATACGGTT had 2 edits."""
    assert bestApproximateMatchEditDistance('GCGTATGC', 'TATTGGCTATACGGTT') == 2

def question1(t):
    """What is the edit distance of the best match between pattern GCTGATCGATCGTACG and the excerpt of human chromosome 1? (Don't consider reverse complements.)"""
    print("Question1: " + str(bestApproximateMatchEditDistance('GCTGATCGATCGTACG', t)))

def question2(t):
    """What is the edit distance of the best match between pattern GATTTACCAGATTGAG and the excerpt of human chromosome 1? (Don't consider reverse complements.)"""
    print("Question2: " + str(bestApproximateMatchEditDistance('GATTTACCAGATTGAG', t)))



def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching a prefix of 'b' that is at least 'min_length' characters long.  If no such overlap exists, return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

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

def overlap_all_pairs(reads, min_length):
    from pprint import pprint
    
    overlap_map = {}
    overlap_graph = {}
    overlap_pairs = []

    #We use a Python dictionary to associate each k-mer with its corresponding set. 
    suffixDict = {}
    for read in reads:
        kmers = getkmers(read, min_length)
        #print(kmers)
        #(1) For every k-mer in a read, we add the read to the set object corresponding to that k-mer.
        for kmer in kmers:
            if not kmer in suffixDict.keys():
                #Let every k-mer in the dataset have an associated Python set object, which starts out empty. 
                suffixDict[kmer] = set()
            suffixDict[kmer].add(read)
    #pprint(suffixDict)

    #(2) Now, for each read a, we find all overlaps involving a suffix of a
    for read in reads:
        #we take a's length-k suffix,
        suffix = read[-min_length:]
        #if len(suffix) < min_length:
        #    continue

        #find all reads containing that k-mer (obtained from the corresponding set) ...
        matching_reads = suffixDict[suffix]

        #...and call overlap(a, b, min_length=k) for each.
        for read2 in matching_reads:
            # The most important point is that we do not call overlap(a, b, min_length=k) if b does not contain the length-k suffix of a.
            #if read2.find(suffix) >= 0 and read2 != suffix :
            if read2 != read :
                val = overlap(read, read2, min_length) 
                if val > 0:
                    overlap_map[ (read, read2) ] = val
                    overlap_graph[read] = read2
                    overlap_pairs.append( (read,read2) )

    #pprint(overlap_map)
    #pprint(overlap_pairs)
    #pprint(overlap_graph)
    return overlap_pairs, overlap_map, overlap_graph



def getkmers(read, kmer_length):
    """our read is GATTA and k=3, we would add GATTA to the set objects for GAT, ATT and TTA."""
    return [ read[i:i+kmer_length] for i in range(len(read)+1-kmer_length) ]

def test2():
    kmers = getkmers('GATTA', 3)
    assert kmers == ['GAT', 'ATT', 'TTA']

def example1():
    reads = ['ABCDEFG', 'EFGHIJ', 'HIJABC']
    assert overlap_all_pairs(reads, 3)[0] == [('ABCDEFG', 'EFGHIJ'), ('EFGHIJ', 'HIJABC'), ('HIJABC', 'ABCDEFG')]
    assert overlap_all_pairs(reads, 4)[0] == []

def example2():
    from pprint import pprint

    reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']

    results4, overlap_map4, overlap_graph4 = overlap_all_pairs(reads, 4)
    expected4 = [('CGTACG', 'TACGTA'),
                 ('CGTACG', 'GTACGT'),
                 ('CGTACG', 'GTACGA'),
                 ('CGTACG', 'TACGAT'),
                 ('TACGTA', 'ACGTAC'),
                 ('TACGTA', 'CGTACG'),
                 ('GTACGT', 'TACGTA'),
                 ('GTACGT', 'ACGTAC'),
                 ('ACGTAC', 'GTACGA'),
                 ('ACGTAC', 'GTACGT'),
                 ('ACGTAC', 'CGTACG'),
                 ('GTACGA', 'TACGAT')]

    assert sorted(results4) ==  sorted(expected4) , "example2, first assert failed"

    results5, overlap_map5, overlap_graph5 = overlap_all_pairs(reads, 5)
    expected5 = [('CGTACG', 'GTACGT'),
                 ('CGTACG', 'GTACGA'),
                 ('TACGTA', 'ACGTAC'),
                 ('GTACGT', 'TACGTA'),
                 ('ACGTAC', 'CGTACG'),
                 ('GTACGA', 'TACGAT')] 

    assert sorted(results5) ==  sorted(expected5), "example2, second assert failed"

def question3and4(reads):
    from pprint import pprint

    overlap_pairs, overlap_map, overlap_graph = overlap_all_pairs(reads, 30)

    """Picture the overlap graph corresponding to the overlaps just calculated. How many edges are in the graph? 
    In other words, how many distinct pairs of reads overlap?"""
    print('Question3: ')
    print(len(overlap_map))

    """Picture the overlap graph corresponding to the overlaps computed for the previous question. 
    How many nodes in this graph have at least one outgoing edge? (In other words, how many reads have a suffix involved in an overlap?)"""
    print('Question4: ')
    print(len(overlap_graph))


    

def main():
    test1()
    print("All tests completed successfully")
    genome = readGenome('chr1.GRCh38.excerpt.fasta')
    question1(genome)
    question2(genome)

    test2()
    example1()
    example2()
    print("All example tests completed successfully")

    reads, qualities = readFastq('ERR266411_1.for_asm.fastq')
    question3and4(reads)
    print("All done")
    



if __name__ == "__main__":
    main()