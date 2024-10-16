#A-Phylogenetic-Tree-Construction-Using-UPGMA-and-Neighbor-Joining-Methods-in-MATLAB

This MATLAB code constructs phylogenetic trees using the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) and Neighbor-Joining methods. It can process both pre-aligned sequences and sequences obtained from GenBank, making it versatile for different types of input data.

The code generates phylogenetic trees using the UPGMA and Neighbor-Joining methods and offers the option of performing bootstrap analyses to assess the trees confidence. The results are visualized by including confidence values in the branches of the tree.

This work builds upon foundational scripts provided by the MATLAB Bioinformatics Toolbox, specifically 'Building a Phylogenetic Tree for the Hominidae Species' and 'Bootstrapping Phylogenetic Trees'. The methodologies and concepts from these scripts were adapted and integrated to enhance the current program's functionality and user-friendliness.

##Features
---Option to use pre-aligned or aligned sequences.
---UPGMA and Neighbor-Joining tree construction methods.
---Ability to set a confidence threshold.

##Requirements
--- MATLAB Bioinformatics Toolbox.


Input Data
For the unaligned sequences part, you need a file. You can use the data provided in this file: [dados.txt](https://github.com/user-attachments/files/17401238/dados.txt). This file is an example for the sequences to be aligned. 
If you use this file,  when it appears ''File name with the information:''  insert  ''dados.txt'' 

To change the aligned sequences used, simply change the names of the sequences in the code (line 304 for the aligned sequences).



IMPORTANT
-- I would like to point out that the sequence alignment part needs improvement. We are currently working to enhance the accuracy of the alignment, which will further increase confidence in the results. 
It is important to note, however, that these issues do not affect the overall structure of the tree.
