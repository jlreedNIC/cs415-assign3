# -------
#
# @file     emission_creation.py
# @author   Jordan Reed
# @class    CS415 Computational Biology
# @brief    This file will create emissions tables
# ----------

filename = "DataFile2.txt"

def open_read_file():
    """
    function to open the file with the given name and read in all the sequences and store them into the appropriate list.

    :return: lists of amino sequences and genome states
    """
    amino_seq = []
    state_seq = []
    temp_state = []

    entire_file = []

    with open(filename, "r") as f:
        # read entire file in
        entire_file = f.readlines()
    
    # separate amino sequence from state sequences
    for i in range(0, len(entire_file), 3): 
        amino_seq.append(entire_file[i][:-1])
        temp_state.append(entire_file[i+1][:-1])
    
    # turn state sequences into int type
    for i in range(0, len(temp_state)):
        state_seq.append([])
        for j in range(0, len(temp_state[0])):
            state_seq[i].append(int(temp_state[i][j]))

    return amino_seq, state_seq

aminos, states = open_read_file()
print(aminos)
print('\n', states)