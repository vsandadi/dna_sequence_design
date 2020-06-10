'''
Function to generate DNA sequences for each bit string.
'''

#initial_sequence = 'ATTCCGAGA'
#mutation_list = [(2, 'T', 'C'), (4, 'C', 'A')]
#bitstring = ['101', '100']

def generating_sequences(initial_sequence, mutation_list, bitstring):
    
    '''
    Inputs: initial_sequence as a string of DNA bases
            mutation_list as list of tuples of mutation site (int), intial base (str), and mutation (str)
            bitstring as list of multiple strings of bits
    Returns: the list of final sequences 
    '''
    
    #Element is one bistring in the list bitstring
    
    final_sequence_list = []
   
    for element in bitstring:
        
        #Raise AssertionError for invalid length of bitstring
       
        assert len(mutation_list) >= len(element), 'Number of mutations has to be greater than or equal to length of bitstring.'
            
        #Turning string inputs into lists for mutability 
        
        sequence_list = list(str(initial_sequence)) 
        
        #Sorting mutation_list by numerical order of mutation placement to match up to correct bitstring value
        
        mutation_list = sorted(mutation_list, key=lambda x: x[0])
        
        #Zipping will pair each element of an index together
        #First bitstring will pair with first tuple 
        #Mutation will iterate through each tuple in zip
        #Bit will iterate through each bitstring value
        
        for mutation,bit in zip(mutation_list, element):
            
            #Unwrapping elements in each tuple
            #Underscore indicates that variable (the original base) will not be needed
            
            mutation_position, _ , final_base = mutation
            if bit == '1':
                sequence_list[mutation_position] = final_base
        
        #Combines mutated DNA bases from list into a string  
        
        new_sequence = "".join(sequence_list)
        final_sequence_list.append(new_sequence)
        
    return final_sequence_list
