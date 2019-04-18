def global_alignment(str_1, str_2):
    import numpy as np

    #imports blosum62 scoring matrix
    from Bio.SubsMat import MatrixInfo as matrixFile
    blosum = matrixFile.blosum62

    v = str_1
    w = str_2

    #create new dictionary from blosum62 dictionary: if blosum62 dictionary has A: B, new dictionary is has A: B and B: A
    blosum62 = {}
    for key in blosum:
        blosum62[key] = blosum[key]
        blosum62[(key[1], key[0])] = blosum[key]

    #create node matrix
    node_matrix = np.zeros((len(v) + 1, len(w) + 1), dtype=int)

    #create backtrack matrix
    backtrack_matrix = np.zeros((len(v) + 1, len(w) + 1), dtype=int)

    #indel penalty (penalty for insertion of empty space, or deletion of a nucleotide)
    indel = 5

    #node matrix zero i and j
    for i in xrange(1, len(node_matrix)):
        node_matrix[i][0] = i * indel * -1
    for j in xrange(0, len(node_matrix[0])):
        node_matrix[0][j] = j * indel * -1

    print(node_matrix)

    #backtrack matrix zero-i and j
    backtrack_matrix[0, :] = 1
    backtrack_matrix[:, 0] = -1
    backtrack_matrix[0][0] = 0

    print(backtrack_matrix)

    for i in xrange(1, len(node_matrix)):
        for j in xrange(1, len(node_matrix[i])):
                # set diagonal value
                diagonal = node_matrix[i - 1][j - 1] + blosum62[v[i - 1], w[j - 1]]
                # pick the max value of the values above, left, and above-left diagonal
                max_adjacent_value = node_matrix[i - 1][j] - indel
                backtrack_value = -1
                if node_matrix[i][j - 1] - indel > max_adjacent_value:
                    max_adjacent_value = node_matrix[i][j - 1] - indel
                    backtrack_value = 1
                if diagonal > max_adjacent_value:
                    max_adjacent_value = diagonal
                    backtrack_value = 0
                node_matrix[i][j] = max_adjacent_value
                backtrack_matrix[i][j] = backtrack_value

    max_alignment_score = node_matrix[len(node_matrix) - 1][len(node_matrix[0]) - 1]

    #BACKTRACKING FOR ALIGNMENT
    va = ""
    wa = ""

    i = len(backtrack_matrix) - 1
    j = len(backtrack_matrix[0]) - 1

    while not i == 0 or not j == 0:
        if backtrack_matrix[i][j] == -1:
            va += v[i - 1]
            wa += "-"
            i -= 1
        elif backtrack_matrix[i][j] == 1:
            wa += w[j - 1]
            va += "-"
            j -= 1
        elif backtrack_matrix[i][j] == 0:
            va += v[i - 1]
            wa += w[j - 1]
            i -= 1
            j -= 1

    #reversing va and wa
    va = va[::-1]
    wa = wa[::-1]

    return va, wa

def longest_common_subsequence(str_1, str_2):
    import numpy as np

    v = str_1
    w = str_2
    score_matrix = np.zeros((len(v) + 1, len(w) + 1), dtype=int)

    #assign scores (we do not have a result string yet)
    #iterate through score matrix
    for row in xrange(1, len(score_matrix)):
        for column in xrange(1, len(score_matrix[row])):
            #do nucleotides align?
            nucleotides_align = False
            if v[row - 1] == w[column - 1]:
                nucleotides_align = True
            #if nucleotides do not align, set diagonal value
            value_if_diagonal = score_matrix[row - 1][column - 1]
            if nucleotides_align:
                value_if_diagonal += 1
            #pick the max value of the values above, left, and above-left diagonal
            max_adjacent_value = max(score_matrix[row - 1][column], score_matrix[row][column - 1], value_if_diagonal)
            score_matrix[row][column] = max_adjacent_value

    '''
    Backtracking for Path Alignment
    first, we must assign each value an arrow score indicating direction of path movement
    down (above value is max_adjacent_value) is -1
    right (left value is max_adjacent_value) is 1
    diagonal (diagonal value is max_adjacent_value) is 0 
    '''
    backtrack_matrix = np.zeros((len(v) + 1, len(w) + 1), dtype=int)

    #fill backtrack_matrix with values already known

    for row in xrange(0, len(backtrack_matrix)):
        for column in xrange(0, len(backtrack_matrix[row])):
            if row == 0 and column > 0:
                backtrack_matrix[row][column] = 1
            elif column == 0 and row > 0:
                backtrack_matrix[row][column] = -1

    #create backtrack matrix
    for row in xrange(1, len(backtrack_matrix)):
        for column in xrange(1, len(backtrack_matrix[row])):
            value = 4
            #selection statement branches must be in this order
            if score_matrix[row][column] == score_matrix[row - 1][column]:
                value = -1
            elif score_matrix[row][column] == score_matrix[row][column - 1]:
                value = 1
            elif (score_matrix[row][column] == score_matrix[row - 1][column - 1] + 1) and v[row - 1] == w[column - 1]:
                value = 0

            backtrack_matrix[row][column] = value

    #backtracking for alignment

    subseq = ""
    row = len(backtrack_matrix) - 1
    column = len(backtrack_matrix[0]) - 1

    while row > 0 and column > 0:
        if backtrack_matrix[row][column] == -1:
            row -= 1
        elif backtrack_matrix[row][column] == 1:
            column -= 1
        elif backtrack_matrix[row][column] == 0:
            subseq += v[row - 1]
            row -= 1
            column -= 1

    rev_subseq = subseq[::-1]

    return rev_subseq

def topo_order(edges_list):
    import random
    input_lines = edges_list

    #list processed into dictionary representing graph
    graph_dict = {}
    #process input file
    for i in xrange(0, len(input_lines)):
        input_lines[i] = input_lines[i].rstrip('\n').split(' -> ')
        input_lines[i][1] = input_lines[i][1].split(',')
        graph_dict[int(input_lines[i][0])] = []
        #process outgoing edges for each node into a list
        for j in xrange(0, len(input_lines[i][1])):
            graph_dict[int(input_lines[i][0])].append(int(input_lines[i][1][j]))

    #if there exist outgoing edges pointing to nodes not in graph_dict, add said nodes to graph_dict
    new_dict_items = {}
    for key, value in graph_dict.iteritems():
        for i in xrange(0, len(value)):
            if (value[i] not in graph_dict) and (value[i] not in new_dict_items):
                new_dict_items[value[i]] = []
    graph_dict.update(new_dict_items)

    output_list = []
    while len(graph_dict) > 0:
        #if a node does not have incoming edges, it is one potential starting point
        #list all elements with no incoming edges
        have_incoming_edges = []
        set_of_nodes = []
        for key, value in graph_dict.iteritems():
            set_of_nodes.append(key)
            #iterates through all outgoing edges
            for i in xrange(0, len(value)):
                have_incoming_edges.append(value[i])
        #remove duplicates
        have_incoming_edges = set(have_incoming_edges)
        set_of_nodes = set(set_of_nodes)

        #use set subtraction to create set of nodes with no incoming edges
        no_incoming_edges = list(set_of_nodes - have_incoming_edges)

        #randomly select a starting node
        starting_node = random.choice(no_incoming_edges)
        output_list.append(starting_node)
        #remove starting node to continue to next step of process
        del graph_dict[starting_node]


    #export list of nodes in topological order
    return output_list

def cl(filename_1, filename_2):
    input_1 = open(filename_1, 'r')
    input_2 = open(filename_2, 'r')

    lines_1 = input_1.readlines()
    lines_2 = input_2.readlines()

    output = ""
    for i in xrange(0, len(lines_1)):
        if lines_1[i] == lines_2[i]:
            output += ("Line " + str(i) + " matches\n")
        else:
            output += ("Line " + str(i) + " does not match, begins with: (file_1)" + str(lines_1[i][0:10]) + " (file_2)" + str(lines_2[i][0:10]) + '\n')

    output += "End of file"
    return output

'''
Find The Reverse Complement of a String
Takes a given string of nucleotides, reverses it, then returns the complement.

Parameters: string
Returns: string
'''
def reverse_complement(strg):
    complement = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }
    strg_reversed = strg[::-1]
    strg_out = ""
    for i in xrange(0, len(strg_reversed)):
        strg_out += complement[strg_reversed[i]]
    return strg_out

'''
#15 Implement PatternToNumber 
Converts a pattern of DNA nucleotides to a base 10 integer representation.

Parameters: string (a given k-mer)
Returns: integer (base 10 index of given k-mer)
'''
def pattern_to_number(pattern):
    indexValue = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3
    }
    if pattern == "":
        return 0
    letter = pattern[-1:]
    prefix = pattern[:-1]
    return 4 * pattern_to_number(prefix) + indexValue[letter]


'''
#16 Implement NumberToPattern 
Converts a base 10 integer representation of DNA nucleotide pattern into a string representation.

Parameters: integer (base 10 index of given k-mer), integer (value of k)
Returns: string
'''
def number_to_pattern(index, k):
    letter = {
        0 : 'A',
        1 : 'C',
        2 : 'G',
        3 : 'T'
    }
    if k == 1:
        return letter[index]
    dividend = index / 4
    r = index % 4
    current_letter = letter[r]
    prefix = number_to_pattern(dividend, k - 1)
    return prefix + current_letter


'''
#17 Find the Most Frequent Words in a String (Faster Frequent Words)
Returns the k-mer or k-mers, given a value for k, that occur most frequently in a 
nucleotide string.

Parameters: string (nucleotide string), integer (value of k)
'''
def most_frequent_words(text, k):
    import string
    import numpy as np

    #frequency of the words
    freq_words = []
    indices = []

    for i in xrange(0, len(text) - k):
        pattern = text[i:i + k]
        indices.append(pattern_to_number(pattern))
    indices.sort()

    counts = np.zeros(len(indices), dtype=int)
    for i in xrange(0, len(indices)):
        if indices[i] == indices[i - 1]:
            counts[i] = counts[i - 1] + 1
        else:
            counts[i] = 1

    max_frequency = np.amax(counts)

    for i in xrange(0, len(indices)):
        if counts[i] == max_frequency:
            pattern = number_to_pattern(indices[i], k)
            freq_words.append(pattern)

    return freq_words