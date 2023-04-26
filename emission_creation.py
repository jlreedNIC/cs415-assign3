# -------
#
# @file     emission_creation.py
# @author   Jordan Reed
# @class    CS415 Computational Biology
# @brief    This file will create emissions tables
#
#           count for each state (0,1,2) what amino is emitted
#           divide each amino count by the total amount for that state
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

def create_emissions(amino_seq, state_seq):
    """
    Function that will create 3 emission probability tables given the list of amino acid sequences and the state sequences

    :param amino_seq: list of amino acid sequences
    :param state_seq: list of genome states corresponding to the amino acids
    :return: emission probabilities for each of the genome states
    """
    state0_emit = [0 for i in range(0,26)] 
    state1_emit = [0 for i in range(0,26)] 
    state2_emit = [0 for i in range(0,26)] 

    # loop through all sequences
    for i in range(0, len(amino_seq)):
        # loop through a sequence and count how many times an amino acid is emitted when we are in each state (0,1,2)
        for j in range(0, len(amino_seq[0])):
            char_num = ord(amino_seq[i][j])-97
            if state_seq[i][j] == 0:
                state0_emit[char_num] += 1
            elif state_seq[i][j] == 1:
                state1_emit[char_num] += 1
            else: # state 2
                state2_emit[char_num] += 1
                
    state0_total = 0
    state1_total = 0
    state2_total = 0

    # count total in each category
    for i in range(0, len(state0_emit)):
        state0_total += state0_emit[i]
        state1_total += state1_emit[i]
        state2_total += state2_emit[i]

    # divide the count by the total in that state
    for i in range(0, len(state0_emit)):
        state0_emit[i] = round(state0_emit[i]/state0_total, 6)
        state1_emit[i] = round(state1_emit[i]/state1_total, 5)
        state2_emit[i] = round(state2_emit[i]/state2_total, 5)

    return state0_emit, state1_emit, state2_emit

def save_emission_tables(filename, state0_emit, state1_emit, state2_emit):
    """
    Function that saves the emission tables to the given file in a csv format, with headers

    :param filename: file in which to save the data
    """

    with open(filename, "w") as f:
        f.write('state_0,')
        for i in range(0, len(state0_emit)):
            f.write(f'{state0_emit[i]},')
        f.write('\n')

        f.write('state_1,')
        for i in range(0, len(state1_emit)):
            f.write(f'{state1_emit[i]},')
        f.write('\n')

        f.write('state_2,')
        for i in range(0, len(state2_emit)):
            f.write(f'{state2_emit[i]},')
        f.write('\n')

aminos, states = open_read_file()
# print(aminos)
# print('\n', states)

state0, state1, state2 = create_emissions(aminos, states)
# print(f'state0: {state0}\n')
# print(f'state1: {state1}\n')
# print(f'state2: {state2}\n')

# save_emission_tables("emissions.csv", state0, state1, state2)