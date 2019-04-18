# Bioinformatics Algorithms Implemented in Python
### by Jonathan Carr

These are some of the Python functions I created in Prof. Rachel Trana's CS 327, Computational Methods in Biology, course at
Northeastern Illinois University, in the Spring 2019 semester.

They implement various algorithms that have been historically useful in bioinformatics.
The problems they solve are accessible on the Rosalind bioinformatics algorithms learning site. Each algorithm has a link to the relevant Rosalind page in its description.

Please do not steal these and submit them as your own work for the course.

## Getting Started
Download the repository, then access the functions inside it from your own Python file using import statements. I wrote these in Python 2.7. As such, they're not guaranteed to work with other versions of Python.

To test the functions, you'll need a text input in the correct format. Please check the Rosalind page for each algorithm to see a sample of the expected input format. Then, the input has to be manipulated into the correct format to be passed into each function. I have left the code to do that out of this repository.

## Contents

### Highest-Scoring Global Alignment of Two Amino Acid Strings (global_alignment)
Finds the alignment of two strings (of amino acids) that is highest-scoring according to the BLOSUM62 scoring matrix.
Uses a dynamic programming algorithm to choose the highest-scoring path through a matrix of scores representing possible alignments.

<strong>Parameters:</strong> string (first string), string (second string)

<strong>Returns:</strong> string, string (strings modified to fit alignment that is highest-scoring for both)

<strong>Rosalind page:</strong> http://rosalind.info/problems/ba5e/

### Longest Common Subsequence of Two Strings (longest_common_subsequence)
Finds the longest string of letters (in the original intent of this problem, nucleotides) contained in the same order in both of two input strings.

<strong>Parameters:</strong> string (first string), string (second string)

<strong>Returns:</strong> string (longest common subsequence)

<strong>Imports:</strong> numpy

<strong>Rosalind page:</strong> http://rosalind.info/problems/ba5c/

### Find Topological Order of a Directed Graph (topo_order)
Given a directed graph, output the nodes in topological order. In other words, start with a random node that has no incoming edges. It becomes first in order. Then, remove it and all of its outgoing edges. Repeat this process until all nodes have been removed. The resulting order is a topological order (there may be multiple ones for a single graph).  
 
Takes a list of edges as a parameter. 
Each element in the list is formatted as follows, 
where 0 represents the node where the edge is outgoing, 
and and 1 represents the node where it is incoming:
```
'0 -> 1'
```
<strong>Parameters:</strong> list

<strong>Returns:</strong> list (numbered nodes in topological order)

<strong>Rosalind page:</strong>http://rosalind.info/problems/ba5n/
### Compare Lines (cl)

This is a tool for comparing two files.
I use it to compare my algorithms' output with a sample output, to
ensure correctness.


<strong>Parameters:</strong> string filename1, string filename2

<strong>Returns:</strong> string (multiple lines). Each line indicates whether the given line matches in both files

### Reverse Complement of a String (reverse_complement)
Takes a given string of nucleotides, reverses it, then returns the complement (replacing all A nucleotides with T, G with C, and vice versa).

<strong>Parameters:</strong> string

<strong>Returns:</strong> string

<strong>Rosalind page:</strong>
http://rosalind.info/problems/ba1c/

### PatternToNumber (pattern_to_number)
Converts a pattern of DNA nucleotides to a base 10 integer representation.

<strong>Parameters:</strong> string (a given k-mer)

<strong>Returns:</strong> integer (base 10 index of given k-mer)
### NumberToPattern (number_to_pattern) 
Converts a base 10 integer representation of DNA nucleotide pattern into a string representation.

<strong>Parameters:</strong> integer (base 10 index of given k-mer), integer (value of k)

<strong>Returns:</strong> string

### Find the Most Frequent Words in a String (most_frequent_words)
Uses the Sorted-Faster-Frequent-Words Algorithm.
Returns the k-mer or k-mers (that is, substrings of length k, given a value for k) that occur most frequently in a 
nucleotide string.

<strong>Parameters:</strong> string (nucleotide string), integer (value of k)

<strong>Returns:</strong> list (most frequent k-mer(s))



