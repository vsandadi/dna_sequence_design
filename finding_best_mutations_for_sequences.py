'''
Function identifying best mutation sites
'''

def identifying_good_mutations(initial_sequence, mutation_list, bitstring, num_mut):
    
    '''
    Inputs: initial_sequence as a string of DNA bases
            mutation_list as list of tuples of mutation site (int), intial base (str), and mutation (str)
            bitstring as list of multiple strings of bits
    Returns: best mutations for every sequence
    '''
    
    #Initialize list for best mutations for every sequence
    
    sequence_mutations = []
    
    #Find sequences
    
    mut_seq_list = generating_sequences(initial_sequence, mutation_list, bitstring)
    
    #Find best mutatitons for every sequence in the list
    
    for sequence in mut_seq_list:
        potential_mut_list = potential_mutations(sequence)
        best_mut = best_mutations(potential_mut_list, num_mut)
        sequence_mutations.append((sequence, best_mut))
        
    return sequence_mutations
        
        