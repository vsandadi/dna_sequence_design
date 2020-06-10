'''
Helper functions to identify mutation sites on a DNA sequence
'''

#Example inputs: sequence_one = 'ATTGCA'
#                sequence_two = 'TAACGT'        
#Example output: mutation_list = [(2, 'T', 'C'), (4, 'C', 'A')]

def return_complement_base(base):
    
    '''
    Helper function for identifying complement base pairs
    Inputs: base as a  string
    Returns: complementary base
    '''
    
    if base == 'G':
        return 'C'
    elif base == 'C':
        return 'G'
    elif base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'

def complement_sequence(sequence):
    
    '''
    Helper function for identifying complement sequence
    Inputs: sequence as a  string
    Returns: complementary sequence as a string
    '''
    sequence_list = list(sequence)
    c_list = []
    for base in sequence_list:
        if base == 'G':
            c_list.append('C')
        elif base == 'C':
            c_list.append('G')
        elif base == 'A':
            c_list.append('T')
        elif base == 'T':
             c_list.append('A') 
    comp_sequence = ''.join([str(elem) for elem in c_list]) 
    return comp_sequence
    
def calculating_g_binding(sequence_one, sequence_two):
    
    '''
    Function for calculating the free energy between two  sequences
    Inputs: initial_sequence as string of DNA bases
    Returns: free binding energy between two sequences
    '''
    
    #Storing the free energies of the nearest neighbor pairs in a dict
    
    nn_g_energy_data = {}
    nn_g_energy_data[('G', 'A'), ('C','A' )] = 0.17
    nn_g_energy_data[('G', 'A'), ('C','C' )] = 0.81
    nn_g_energy_data[('G', 'A'), ('C','G' )] = -0.25
    nn_g_energy_data[('G', 'A'), ('C','T' )] = -1.30
    nn_g_energy_data[('G', 'C'), ('C','A' )] = 0.47
    nn_g_energy_data[('G', 'C'), ('C','C' )] = 0.79
    nn_g_energy_data[('G', 'C'), ('C','G' )] = -2.24
    nn_g_energy_data[('G', 'C'), ('C','T' )] = 0.62
    nn_g_energy_data[('G', 'G'), ('C','A' )] = -0.52
    nn_g_energy_data[('G', 'G'), ('C','C' )] = -1.84
    nn_g_energy_data[('G', 'G'), ('C','G' )] = -1.11
    nn_g_energy_data[('G', 'G'), ('C','T' )] = 0.08
    nn_g_energy_data[('G', 'T'), ('C','A' )] = -1.44
    nn_g_energy_data[('G', 'T'), ('C','C' )] = 0.98
    nn_g_energy_data[('G', 'T'), ('C','G' )] = -0.59
    nn_g_energy_data[('G', 'T'), ('C','T' )] = 0.45
    nn_g_energy_data[('C', 'A'), ('G','A' )] = 0.43
    nn_g_energy_data[('C', 'A'), ('G','C' )] = 0.75
    nn_g_energy_data[('C', 'A'), ('G','G' )] = 0.03
    nn_g_energy_data[('C', 'A'), ('G','T' )] =-1.45 
    nn_g_energy_data[('C', 'C'), ('G','A' )] = 0.79
    nn_g_energy_data[('C', 'C'), ('G','C' )] = 0.70
    nn_g_energy_data[('C', 'C'), ('G','G' )] = -1.84
    nn_g_energy_data[('C', 'C'), ('G','T' )] = 0.62
    nn_g_energy_data[('C', 'G'), ('G','A' )] = 0.11
    nn_g_energy_data[('C', 'G'), ('G','C' )] = -2.17
    nn_g_energy_data[('C', 'G'), ('G','G' )] = -0.11
    nn_g_energy_data[('C', 'G'), ('G','T' )] = -0.47
    nn_g_energy_data[('C', 'T'), ('G','A' )] = -1.28
    nn_g_energy_data[('C', 'T'), ('G','C' )] = 0.40
    nn_g_energy_data[('C', 'T'), ('G','G' )] = -0.32
    nn_g_energy_data[('C', 'T'), ('G','T' )] = -0.12
    nn_g_energy_data[('A', 'A'), ('T','A' )] = 0.61  
    nn_g_energy_data[('A', 'A'), ('T','C' )] = 0.88
    nn_g_energy_data[('A', 'A'), ('T','G' )] = 0.14
    nn_g_energy_data[('A', 'A'), ('T','T' )] = -1.00
    nn_g_energy_data[('A', 'C'), ('T','A' )] = 0.77
    nn_g_energy_data[('A', 'C'), ('T','C' )] = 1.33
    nn_g_energy_data[('A', 'C'), ('T','G' )] = -1.44
    nn_g_energy_data[('A', 'C'), ('T','T' )] = 0.64
    nn_g_energy_data[('A', 'G'), ('T','A' )] = 0.02
    nn_g_energy_data[('A', 'G'), ('T','C' )] = -1.30
    nn_g_energy_data[('A', 'G'), ('T','G' )] = -0.13
    nn_g_energy_data[('A', 'G'), ('T','T' )] = 0.71
    nn_g_energy_data[('A', 'T'), ('T','A' )] = -0.88
    nn_g_energy_data[('A', 'T'), ('T','C' )] = 0.73
    nn_g_energy_data[('A', 'T'), ('T','G' )] = 0.07
    nn_g_energy_data[('A', 'T'), ('T','T' )] = 0.69
    nn_g_energy_data[('T', 'A'), ('A','A' )] = 0.69   
    nn_g_energy_data[('T', 'A'), ('A','C' )] = 0.92
    nn_g_energy_data[('T', 'A'), ('A','G' )] = 0.42
    nn_g_energy_data[('T', 'A'), ('A','T' )] = -0.58
    nn_g_energy_data[('T', 'C'), ('A','A' )] = 1.33
    nn_g_energy_data[('T', 'C'), ('A','C' )] = 1.05
    nn_g_energy_data[('T', 'C'), ('A','G' )] = -1.28
    nn_g_energy_data[('T', 'C'), ('A','T' )] = 0.97
    nn_g_energy_data[('T', 'G'), ('A','A' )] = 0.74
    nn_g_energy_data[('T', 'G'), ('A','C' )] = -1.45
    nn_g_energy_data[('T', 'G'), ('A','G' )] = 0.44
    nn_g_energy_data[('T', 'G'), ('A','T' )] = 0.43
    nn_g_energy_data[('T', 'T'), ('A','A' )] = -1.00
    nn_g_energy_data[('T', 'T'), ('A','C' )] = 0.75
    nn_g_energy_data[('T', 'T'), ('A','G' )] = 0.34
    nn_g_energy_data[('T', 'T'), ('A','T' )] = 0.68
    nn_g_energy_data[('A', 'C'), ('A', 'G')] = 0.17
    nn_g_energy_data[('A', 'C'), ('C', 'G')] = 0.47
    nn_g_energy_data[('A', 'C'), ('G', 'G')] = -0.52
    nn_g_energy_data[('C', 'C'), ('A', 'G')] = 0.81
    nn_g_energy_data[('C', 'C'), ('C', 'G')] = 0.79
    nn_g_energy_data[('C', 'C'), ('T', 'G')] = 0.98
    nn_g_energy_data[('G', 'C'), ('A', 'G')] = -0.25
    nn_g_energy_data[('G', 'C'), ('G', 'G')] = -1.11
    nn_g_energy_data[('G', 'C'), ('T', 'G')] = -0.59
    nn_g_energy_data[('T', 'C'), ('C', 'G')] = 0.98
    nn_g_energy_data[('T', 'C'), ('G', 'G')] = -0.59
    nn_g_energy_data[('T', 'C'), ('T', 'G')] = 0.45
    nn_g_energy_data[('A', 'G'), ('A', 'C')] = 0.43
    nn_g_energy_data[('A', 'G'), ('C', 'C')] = 0.79
    nn_g_energy_data[('A', 'G'), ('G', 'C')] = 0.11
    nn_g_energy_data[('C', 'G'), ('A', 'C')] = 0.75
    nn_g_energy_data[('C', 'G'), ('C', 'C')] = 0.70
    nn_g_energy_data[('C', 'G'), ('T', 'C')] = 0.40
    nn_g_energy_data[('G', 'G'), ('A', 'C')] = 0.03
    nn_g_energy_data[('G', 'G'), ('G', 'C')] = -0.11
    nn_g_energy_data[('G', 'G'), ('T', 'C')] = -0.32
    nn_g_energy_data[('T', 'G'), ('C', 'C')] = 0.62
    nn_g_energy_data[('T', 'G'), ('G', 'C')] = -0.47
    nn_g_energy_data[('T', 'G'), ('T', 'C')] = -0.12
    nn_g_energy_data[('A', 'T'), ('A', 'A')] = 0.61
    nn_g_energy_data[('A', 'T'), ('C', 'A')] = 0.77
    nn_g_energy_data[('A', 'T'), ('G', 'A')] = 0.02
    nn_g_energy_data[('C', 'T'), ('A', 'A')] = 0.88
    nn_g_energy_data[('C', 'T'), ('C', 'A')] = 1.33
    nn_g_energy_data[('C', 'T'), ('T', 'A')] = 0.73
    nn_g_energy_data[('G', 'T'), ('A', 'A')] = 0.14
    nn_g_energy_data[('G', 'T'), ('G', 'A')] = -0.13
    nn_g_energy_data[('G', 'T'), ('T', 'A')] = 0.07
    nn_g_energy_data[('T', 'T'), ('C', 'A')] = 0.64
    nn_g_energy_data[('T', 'T'), ('G', 'A')] = 0.71
    nn_g_energy_data[('T', 'T'), ('T', 'A')] = 0.69  
    nn_g_energy_data[('A', 'A'), ('A', 'T')] = 0.69
    nn_g_energy_data[('A', 'A'), ('C', 'T')] = 1.33
    nn_g_energy_data[('A', 'A'), ('G', 'T')] = 0.74
    nn_g_energy_data[('C', 'A'), ('A', 'T')] = 0.92
    nn_g_energy_data[('C', 'A'), ('C', 'T')] = 1.05
    nn_g_energy_data[('C', 'A'), ('T', 'T')] = 0.75
    nn_g_energy_data[('G', 'A'), ('A', 'T')] = 0.42
    nn_g_energy_data[('G', 'A'), ('G', 'T')] = 0.44
    nn_g_energy_data[('G', 'A'), ('T', 'T')] = 0.34
    nn_g_energy_data[('T', 'A'), ('C', 'T')] = 0.97
    nn_g_energy_data[('T', 'A'), ('G', 'T')] = 0.43
    nn_g_energy_data[('T', 'A'), ('T', 'T')] = 0.68
    
    #Initialize lists and variables
    
    sequence_one_list = list(sequence_one)
    sequence_two_list = list(sequence_two)
    nn_energies = 0
    
    #Create the pairs for the initial sequence for nearest neighbors
    
    nn_one = list(zip(sequence_one_list, sequence_one_list[1:] + sequence_one_list[:1])) 
    del nn_one[-1]

    #Create the pairs for the complementary sequence for the nearest neighbors
    
    nn_two = list(zip(sequence_two_list, sequence_two_list[1:] + sequence_two_list[:1])) 
    del nn_two[-1]

    #Create the nearest neighbor pairings
    
    nearest_neighbor = list(zip(nn_one, nn_two)) 
      
    #Calculate the free energies from the nearest neighbor pairings
    
    for nn in nearest_neighbor :
        for i in nn_g_energy_data.keys():
            if nn == i:
                x = nn_g_energy_data[i]
                nn_energies += x

    #Set the terminal_at value based on whether strand ends or begins with AT pair
    
    terminal_at = 0
    if (nearest_neighbor[0][0][0], nearest_neighbor[0][1][0]) == ('A','T') or (nearest_neighbor[0][0][0], nearest_neighbor[0][1][0]) == ('T','A'):
        terminal_at += 0.05
    if (nearest_neighbor[-1][0][1], nearest_neighbor[-1][1][1]) == ('A', 'T') or (nearest_neighbor[-1][0][1], nearest_neighbor[-1][1][1]) == ('T', 'A'):
        terminal_at += 0.05
 
    #Add up free energy of the sequence

    initiation = 1.96
    delta_g = initiation + nn_energies + terminal_at
    
    return delta_g


def potential_mutations(sequence):
    
    '''
    Function for identifying all potential mutations in a sequence
    Inputs: sequence
    Returns: list of the position, original base, and new base of potential mutatations
    '''
   
    potential_mutations_list = []
    all_mutations = []
    epsilon = float(input("Enter an epsilon value: "))
    
    #Create list of all mutated sequences
    
    for i, nt in enumerate(sequence):
        for nt_mut in 'ATCG':
            if nt != nt_mut:
                all_mutations.append((i, nt, nt_mut, sequence[:i] + nt_mut + sequence[i+1:]))      
        
        
    #Compare delta G's for every mutated sequence
    
    for mut_sequence in all_mutations:
        
        y = mut_sequence[3]
        y_complement = complement_sequence(y)
        x = sequence 
        x_complement = complement_sequence(x)
        
        delta_g_xx = calculating_g_binding(x, x_complement)
        delta_g_yy = calculating_g_binding(y, y_complement)
        delta_g_xy = calculating_g_binding(x, y_complement)
        delta_g_yx = calculating_g_binding(y, x_complement)
        
        delta_g1 = delta_g_xy - delta_g_xx
        delta_g2 = delta_g_xy - delta_g_yy
        delta_g3 = delta_g_yx - delta_g_yy
        delta_g4 = delta_g_yx - delta_g_xx
        
        #If delta G values are the same, will be added to the potential mutations list
        
        if max(delta_g1, delta_g2, delta_g3, delta_g4) - min(delta_g1, delta_g2, delta_g3, delta_g4) < epsilon:
            potential_mutations_list.append((y[0], nt, nt_mut))

    return potential_mutations_list

def best_mutations(potential_mutations_list, num_mut):
    
    '''
    Function for identifying the best set of mutations
        Mutations must be as separated from each other as possible
    Inputs: list of potential mutations
            num_mut as number of mutations
    Returns: list of the position, original base, and new base for best mutatations
    '''
    
    #Create list of indicies of the mutations
    
    mutation_pos = []
    for mutation in potential_mutations_list:
        mutation_pos.append(mutation[0])
    
    n = len(mutation_pos)

    
    def isFeasible(mid, mutation_pos, n, num_mut): 
        
        '''
        Function to determine if possible to arrange k elements of mutation_pos[0..n-1] 
            with minimum distance given as mid
        Inputs: mid as minimum distance
                mutation_pos as list of indices of mutations
                n as length of mutation_pos
                num_mut as number of mutations 
        Returns: True and list of mutation positions if possible
        '''
        #Initialize list for best mutation positions 
        # Place first element at arr[0] position 
        # Initialize count of elements placed
        
        pos = mutation_pos[0]
        best_mut_pos = [pos] 
        elements = 1
      
        # Try placing k elements with minimum distance mid.
        
        for i in range(1, n, 1): 
            if (mutation_pos[i] - pos >= mid): 
                  
                # Place next element if its distance from the previously placed element is greater than current mid 
                # Add element to mutation list 
                
                pos = mutation_pos[i] 
                best_mut_pos.append(pos)
                elements += 1
      
                # Return if all elements are placed successfully 
                
                if (elements == num_mut): 
                    return True, best_mut_pos
        return 0, best_mut_pos
    
    def largestMinDist(mutation_pos, n, num_mut): 
        
        '''
        Function to return largest minimum distance for k elements 
            in mutation_pos[0..n-1]. If elements can't be placed, returns empty list
        Inputs: mutation_pos as list of indices of mutations
                n as length of mutation_pos
                num_mut as number of mutations
        Returns: largest minimum distance
        '''
                  
        # Initialize result
        
        best_mut_pos = []
      
        # Consider the maximum possible distance 
        
        left = mutation_pos[0] 
        right = mutation_pos[n - 1] 
      
        # Do binary search for largest minimum distance 
        
        while (left < right): 
            mid = (left + right) / 2
      
            # If it is possible to place k elements with minimum distance mid, search for higher distance
            
            feasible, pos = isFeasible(mid, mutation_pos, n, num_mut)
            if feasible: 
                  
                # Change value of variable max to mid iff 
                # all elements can be successfully placed 
                
                left = mid + 1
                best_mut_pos = pos
      
            # If not possible to place k elements,search for lower distance 
            
            else: 
                right = mid 
      
        return best_mut_pos 
    
    #Create list for best mutation positions and best mutations
    
    best_mutations = []
    best_mut_pos = largestMinDist(mutation_pos, n, num_mut)
    
    #Add best mutations to list 
    
    for index in best_mut_pos:
        for mutation in potential_mutations_list:
            if index == mutation[0]: 
                best_mutations.append(mutation)
        
    return best_mutations


    


