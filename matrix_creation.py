# -------
#
# @file     matrix_creation.py
# @author   Jordan Reed
# @class    CS415 Computational Biology
# @brief    This file will create a scoring matrix based off of the sequences in the datafile.
#
#           The algorithm in which this is done is for each column, all of the pairs are 
#           counted (eg a-b, s-k). This starts with the first row compared to the rest, 
#           then the second row compared to everything below it (discounting the first row) 
#           and so on. After a count is received, those numbers are divided by all the possible 
#           pair combinations for that column. We then take the log of these probabilites.
#
#           The log odds probabilities for each column are then summed together to get the substitution score.
#
#           note on scoring matrix indeces:
#           scoring matrix is build upon the understanding that 0=a, 1=b, 2=c, etc.
#           formula: ord('a')-97
#
#
#           EDIT: 
#           fixed substitution matrix so it follows the algo described in the paper and on https://rob-p.github.io/CSE549F17/lectures/Lec10.pdf
#           Now the diagonals have positive numbers and positives are possible in other places
#           
#           Algo summary:
#           score(ab) = 1/lambda * log2(q(ab)/e(ab))
#           lambda = .5
#           q(ab) = sum of a-b pairs / sum of all possible pairs
#           e(ab) = p(a)^2 if a=b or 2*p(a)*p(b)
#           p(a) = q(aa) + sum of all q(aX) where a!=X and X is every letter
# ----------

import math

filename = "DataFile1-1.txt"
# create empty scoring matrix

def open_read_file():
    """
    function to open the file with the given name and read in all the sequences and store them into one list.

    :return: list of sequences stored in file
    """
    seq = []

    with open(filename, "r") as f:
        # # read entire file in
        # seq = f.readlines()
        # for i in range(0,len(seq)):
        #     seq[i] = seq[i][:-1]
        
        # read first 20 lines in
        for i in range(0,160):
            seq.append(f.readline())
            seq[i] = seq[i][:-1]
    return seq

def print_matrix(matrix, gap=3, headers=True, isFloat=False):
    """
    Prints a matrix (probably the scoring matrix), in a human readable 
    format. Will print headers on the top and on the sides.

    :param matrix: the matrix (or list of lists) to print
    :param gap: spacing for each element of the matrix, defaults to 3
    :param headers: whether or not to print headers (characters a-z), defaults to True
    """
    # print 'a' 'b' ... 'z'
    if headers:
        print(f'{" ":>{gap}}', end='')
        for i in range(0,len(matrix[0])):
            print(f'{chr(i+97):>{gap}}', end='')
        print()

    for i in range(0, len(matrix[0])):
        # print 'a' ... 'z' before each row
        if headers:
            print(f'{chr(i+97):>{gap}}', end='')
        for j in range(0, len(matrix)):
            if isFloat:
                print(f'{matrix[i][j]:{gap+1}.3f}', end='')
            else:
                print(f'{matrix[i][j]:{gap}}', end='')
        print()

def save_scoring_matrix(filename):
    """
    Function that saves the scoring matrix to the given file in a csv format

    :param filename: file in which to save the scoring matrix
    """
    global scoring_matrix

    with open(filename, "w") as f:
        for i in range(0,len(scoring_matrix[0])):
            for j in range(0,len(scoring_matrix)):
                f.write(f'{scoring_matrix[i][j]},')
            f.write('\n')

def create_probabilities(score_matrix):
    """
    helper function that creates the probabilites for each letter for the sequences. 
    needs scoring matrix after all counts have been totaled and divided by total possible pairs

    :param score_matrix: matrix
    :return: list of probabilities
    """
    # create probabilities for each letter
    probabilities = [0 for i in range(0,26)]
    for i in range(0, len(probabilities)):
        probabilities[i] = score_matrix[i][i]

        j = 0
        while j != i:
            probabilities[i] += score_matrix[j][i]
            j += 1
    return probabilities

def exp_pair_probabilities(probabilities):
    """
    turn probabilities for each letter into probabilities for each pair of letters

    :param probabilities: list of probabilities
    :return: matrix for each pair
    """
    # get expected probabilities
    exp_prob = [[0 for j in range(0,26)] for i in range(0,26)]
    for i in range(0, len(exp_prob[0])):
        for j in range(0, len(exp_prob)):
            if i == j:
                exp_prob[i][j] = probabilities[i]*probabilities[i]
            else:
                exp_prob[i][j] = 2 * probabilities[i] * probabilities[j]
                exp_prob[j][i] = exp_prob[i][j]
    return exp_prob

def create_scoring_matrix(seq):
    """
    Better function to create a scoring matrix like the BLOSUM50. See formula at top of file for description. 
    Based on given sequences.

    :param seq: list of sequences to construct matrix
    :return: scoring matrix
    """
    score_matrix = [[0 for j in range(0,26)] for i in range(0,26)]

    # count pairs in sequence list
    for i in range(0, len(seq[0])):
        for j in range(0, len(seq)):
            for k in range(j+1, len(seq)):
                # change current letters into indexes for scoring table
                indexi = ord(seq[j][i])-97
                indexk = ord(seq[k][i])-97

                # count occurance of pair
                score_matrix[indexi][indexk] += 1
                score_matrix[indexk][indexi] += 1
    
    total_pos_pairs = (len(seq)-1) * len(seq) /2 * len(seq[0])

    # divide every count by total possible pairs
    for i in range(0, len(score_matrix[0])):
        for j in range(0, len(score_matrix)):
            score_matrix[i][j] /= total_pos_pairs

    probabilities = create_probabilities(score_matrix)

    exp_prob = exp_pair_probabilities(probabilities)

    # get actual scores
    for i in range(0, len(score_matrix[0])):
        for j in range(0, len(score_matrix)):

            if exp_prob[i][j] == 0 or score_matrix[i][j] == 0:
                score_matrix[i][j] = 0
            else:
                score_matrix[i][j] = round(2 * math.log( (score_matrix[i][j]/exp_prob[i][j]) , 2))
    
    return score_matrix

# ------ main ----------
seq = open_read_file()
scoring_matrix = create_scoring_matrix(seq)
print_matrix(scoring_matrix, 4)
# save_scoring_matrix("real_all.csv")
save_scoring_matrix("real_80percent.csv")
