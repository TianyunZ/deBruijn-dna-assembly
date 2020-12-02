def entropy_score(list_motifs):
    """The entropy of a motif matrix is defined as the sum of the
    entropies of its columns.
    The entropy of each column is a measurement of the uncertainty
    of a probability distribution p. And it's defined as:
                          N
        H(p[1],…,p[N]) = −∑ (p[i]·log2(p[i]))
                         i=1
    """
    rows = len(list_motifs)
    columns = len(list_motifs[0])
    entropy = 0
    for i in range(0, columns):
        count_column = {
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0
        }
        column_entropy = 0
        for j in range(0, rows):
            nucleotide = list_motifs[j][i]
            count_column[nucleotide] += 1
        for count in count_column.values():
            prob_dist = count/rows
            if prob_dist > 0:
                column_entropy -= prob_dist * math.log2(prob_dist)
        entropy += column_entropy
    return entropy
