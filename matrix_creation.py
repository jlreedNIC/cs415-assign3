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
# ----------

import math

filename = "DataFile1-1.txt"
# create empty scoring matrix
scoring_matrix = [[0 for j in range(0,26)] for i in range(0,26)]

def open_read_file():
    """
    function to open the file with the given name and read in all the sequences and store them into one list.

    :return: list of sequences stored in file
    """
    seq = []

    with open(filename, "r") as f:
        # read entire file in
        seq = f.readlines()
        for i in range(0,len(seq)):
            seq[i] = seq[i][:-1]
        
        # read first 20 lines in
        # for i in range(0,20):
        #     seq.append(f.readline())
        #     seq[i] = seq[i][:-1]
    return seq

def print_matrix(matrix, gap=3, headers=True):
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
            print(f'{matrix[i][j]:{gap}}', end='')
        print()

def col_by_col_score(seq):
    """
    Function to assign a score to each possible pair based on the given list of sequences.

    The scoring algorithm is the same as is in Chapter 2 of our book.
    For each column, the actual pairs are divided by the possible number of pairs. Then the 
    logarithmic function is applied to normalize the scores.
    The log odds probability for each column is then summed together to get the final score.

    :param seq: list of sequences
    """
    global scoring_matrix
    temp_scoring = [[0 for i in range(0,26)] for j in range(0,26)]

    # for each column
    for i in range(0, len(seq[0])):
        # count pairs
        for j in range(0,len(seq)): # row
            for k in range(j+1, len(seq)): # row
                # put character into index for matrix
                x = ord(seq[j][i])-97
                y = ord(seq[k][i])-97

                # count each occurance of pair
                temp_scoring[x][y] += 1
                temp_scoring[y][x] += 1
        
        total = 25*26/2
        # divide each count by total possible then add score to scoring matrix
        # reset temp matrix when done
        for j in range(0,len(temp_scoring[0])):
            for k in range(0,len(temp_scoring)):
                if temp_scoring[j][k] != 0:
                    temp_scoring[j][k] /= total
                    scoring_matrix[j][k] += int(math.log(temp_scoring[j][k]))

                    temp_scoring[j][k] = 0

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

seq = open_read_file()
col_by_col_score(seq)
print_matrix(scoring_matrix, 4)
# save_scoring_matrix("all_sequence_scores.csv")

