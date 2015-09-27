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


def main():
    test1()
    print("All tests completed successfully")
    genome = readGenome('chr1.GRCh38.excerpt.fasta')
    question1(genome)
    question2(genome)
    

if __name__ == "__main__":
    main()