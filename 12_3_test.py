#!/usr/bin/env python
# coding: utf-8

# <span style="font-size: 4em;">Welcome to PyTreeMotif</span>
# <p></p>
# 
# #### an enumerative motif finding algorithm that makes trees

# # IMPORTS 
# 
# Dependencies, can't live with em, can't live without em

# In[34]:


import os
import sys

from collections import deque
import random
from tqdm import tqdm
import re
import pandas as pd
import numpy as np


#  # CLASSES 
#  
# ### Define custom classes needed for PyTreeMotif

# ### Sequence_Node
# 
# * nodes representing lmer sequences similar to a reference pair <p></p>
# 
# * these are candidates from which tree nodes will be selected <p></p>
#  
# * records sequence, position, and the input sequence in which it was found <p></p>

# In[5]:


class Sequence_Node:
    def __init__(self, seq, id, pos):
        self.seq = seq
        self.id = id
        self.pos = pos


# ### Sequence_Node_Set
# 
# * a set of Sequence_Nodes from one input sequence that are similar to a reference pair
# 

# In[6]:


class Sequence_Node_Set:
    def __init__(self, id):
        self.id = id
        self.members = set()
        
    def add_member(self, member):
        self.members.add(member)


# ### Reference_Pair
# 
# * reference pair of lmer sequences differing by at most 2d <p></p>
# 
# * one is from input sequence 0, the other from input sequence 1 <p></p>
# 
# * reference pairs are the basis on which all other possible nodes are selected <p></p>
# 
# * contains a pair of sequences, their positions, and an id number for the pair <p></p>
# 
# * contains a list of Sequence_Node_Sets, one for each of the remaining input sequences <p></p>
# 
# * after tree construction, contains a list of trees built from this reference pair <p></p>

# In[7]:


class Reference_Pair:
    def __init__(self, pair_id, seq1, seq2, pos1, pos2, num_seq):
        self.pair_id = pair_id
        self.seq1 = seq1
        self.seq2 = seq2
        self.pos1 = pos1
        self.pos2 = pos2
        
        self.node_sets = [] 
        for i in range(2, num_seq):
            self.node_sets.append(Sequence_Node_Set(i))
            
        self.tree_list = []

        
    


# ### Tree_Node
# 
# * nodes in a motif tree represent lmer sequences which are at most 2d away from all other nodes in the path to the root <p></p>
# * contains lmer sequence, position, input sequence origin, depth, children, and parent <p></p>
# * child nodes can be added and removed <p></p>

# In[8]:


class Tree_Node:
    def __init__(self, seq, seq_id, pos, depth=0):
        self.seq = seq
        self.seq_id = seq_id
        self.pos = pos
        self.depth = depth
        
        self.children = []
        self.parent = None

    def add_child(self, child_node):
        child_node.parent = self
        child_node.depth = self.depth + 1
        self.children.append(child_node)
        

    def remove_child(self, child_node):
        if child_node in self.children:
            child_node.parent = None
            self.children.remove(child_node)


# ### Tree
# 
# * constructed trees represent one or more cliques of motifs <p></p>
# * one path from root to leaf represents a clique of motif instances <p></p>
# * trees with multiple paths contains multiple cliques <p></p><p></p>
# * the root node is a tree node with lmer at most 2d from a reference pair <p></p>
# * all root nodes are from sequence 2 (reference pairs are from sequences 0 and 1) <p></p>
# * each level of the tree represents nodes from each input sequence (besides the ref seqs) <p></p>
# * leaves are stored for to be easily checked during branch pruning <p></p>

# In[9]:


class Tree:
    def __init__(self, root):
        self.root = root
        self.current_leaves = []
        self.height = 0
        self.id = 0
    
    def save_leaf(self, leaf):
        self.current_leaves.append(leaf)
        
    def unsave_leaf(self, leaf):
        if leaf in self.current_leaves:
            self.current_leaves.remove(leaf)


# ### Clique
# 
# * a set of instances of a motif which are all at most 2d different <p></p>
# * a clique is initially formed by finding a path through a tree from root to leaf <p></p>
# * cliques can be merged if all lmers in one are at most 2d from all lmers in another <p></p>
# * merged cliques will then contain more instances than the number of input sequences <p></p>
# * for dedicated testing inputs where there is deliberately one motif per input sequence, merging may pick up random instances only similar by chance <p></p>
# * when number of occurrences of motif per input is unknown, this may help detect multiple motif instances in the same input sequence, but other modifications to the algorithm are required in that scenario <p></p>

# In[10]:


class Clique:
    def __init__(self, path):
        self.nodes = path
        self.seqs = []
        for node in self.nodes:
            self.seqs.append(node.seq)
            
        self.consensus = ""
        self.score = -1
        self.rank = -1
        
        self.merged = False
        
    # def __deepcopy__(self, memo):
    #     new_nodes = []
    #     new_seqs = []
    #     for node in self.nodes:
    #         new_nodes.append(Tree_Node(node.seq, node.seq_id, node.pos))
    #         new_seqs.append(node.seq)
    #     new_copy = Clique(new_nodes)
    #     new_copy.seqs = new_seqs
    #     return new_copy
            


# # FUNCTIONS #
# 
# ### Define subroutines to be used by the main steps of PyTreeMotif ###

# ### read_file_to_list 
# 
# * reads a text file of one input sequence per line into a python list of input sequences
# 

# In[11]:


def read_file_to_list(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
        lines = [line.strip() for line in lines]
        
        return lines
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return []
    except Exception as e:
        print(f"Error: {e}")
        return []


# ### hamming_distance
# 
# * measures hamming distance, the number of single character differences, between two strings
# 

# In[12]:


def hamming_distance(str1, str2):

    return sum(bit1 != bit2 for bit1, bit2 in zip(str1, str2))
    


# ### extendable
# 
# * given a tree and candidate node, finds all possible leaves to which the node could be added <p></p>
# * breadth-first traversal of the tree, adding all children that are similar to node to queue <p></p>
# * if the queue is empty, then the node could not be added to any leaves <p></p>
# * traversal ends after finishing second to lowest level <p></p>
# * all nodes left in the queue are the leaves of branches where every node in the branch all the way to the root is at most 2d from the new node <p></p>

# In[13]:


def extendable(node, tree, ith_seq, _2d):
    queue = deque()
    
    # node must first be similar to root 
    if hamming_distance(tree.root.seq, node.seq) <= _2d:
        queue.append(tree.root)
        
        # while there are still potential branches
        while len(queue) > 0 and queue[0].depth < ith_seq - 3:
            front = queue.popleft()
            
            # check if node is similar to any child nodes
            for child in front.children:
                if hamming_distance(child.seq, node.seq) <= _2d:
                    queue.append(child)
    return queue
        


# ### prune_branches
# 
# * blah blah

# In[14]:


def prune_branches(tree):
    for leaf in tree.current_leaves:
        if leaf.depth < tree.height:
            current_node = leaf
            num_siblings = len(current_node.parent.children) - 1
            branch_list = [leaf.seq]
            while num_siblings == 0:
                current_node = current_node.parent
                num_siblings = len(current_node.parent.children) - 1
                branch_list.append(current_node.seq)
            parent = current_node.parent
            parent.remove_child(current_node)
            tree.unsave_leaf(leaf)
            


# ### check_tree
# 
# * blah blah

# In[15]:


def check_tree(tree):
    is_valid = True
    for leaf in tree.current_leaves:
        if leaf.depth < tree.height:
            is_valid = False
            break
    
    return is_valid


# ### find_paths
# 
# * blah blah

# In[16]:


def find_paths(node, path, paths):
    if node is None:
        return
    path.append(node)
    if not node.children:
        clique = Clique(path)
        paths.append(clique)
    
    for child in node.children:
        # branch_path = copy.deepcopy(path)
        find_paths(child, path[:], paths)
        
    
    
    return paths


# ### merge
# 
# * UNDER CONSTRUCTION

# In[17]:


# def merge(new_cliques, _2d):
#     
#     for prev_clique in merged_cliques:
# 
#         match = True
#         for seq_1 in prev_clique.seqs:
#             if not match:
#                 break
#             for seq_2 in new_clique.seqs:
#                 if hamming_distance(seq_1, seq_2) > _2d:
#                     match = False
#                     break
# 
#         if match:
#             merged_clique = copy.deepcopy(prev_clique)
#             for node in new_clique.nodes:
#                 if node.seq not in prev_clique.seqs:
#                     copy_node = Tree_Node(node.seq, node.seq_id, node.pos)
#                     merged_clique.nodes.append(copy_node)
#                     merged_clique.seqs.append(node.seq)
#             new_cliques.append(merged_clique)
#     new_cliques.append(new_clique)
#     merged_cliques.extend(new_cliques)


# ### score_motifs
# 
# * blah blah

# In[18]:


def score_motifs(motifs, motif_length):
    scores = []
    for motif in motifs:
        
        profile = [[0] * motif_length for _ in range(4)]
        indices = {
            'A': 0,
            'C': 1,
            'G': 2,
            'T': 3,
            0: 'A',
            1: 'C',
            2: 'G',
            3: 'T'
        }
        for i in range(len(motif.seqs)):
            for j in range(motif_length):
                profile[indices[motif.seqs[i][j]]][j] += 1
                
        consensus = ""
                
        for j in range(motif_length):
            max_base = max(row[j] for row in profile)
            
            max_base_index = next(idx for idx, row in enumerate(profile) if row[j] == max_base)
            consensus += indices[max_base_index]
            
        score = 0
        for seq in motif.seqs:
            score += hamming_distance(seq, consensus)
        
        motif.consensus = consensus
        motif.score = score
        scores.append((score, motif))
    
    sorted_scores = sorted(scores, key=lambda x: x[0])
    
    rank = 1
    sorted_scores[0][1].rank = rank

    for i in range(1, len(sorted_scores)):
        if sorted_scores[i][0] > sorted_scores[i - 1][0]:
            rank += 1
        sorted_scores[i][1].rank = rank        
            
        


# ### consensus_grouping
# 
# * blah blah

# In[19]:


def consensus_grouping(motifs):
    consensus_groups = {}
    for motif in motifs:
        if motif.consensus not in consensus_groups:
            consensus_groups[motif.consensus] = [motif]
        else:
            consensus_groups[motif.consensus].append(motif)
            
    consensus_score_report = {}
    
    for consensus, motif_list in consensus_groups.items():
        avg_score = sum(motif.score for motif in motif_list) / len(motif_list)
        best_score = min(motif.score for motif in motif_list)
        
        consensus_score_report[consensus] = [best_score, avg_score]
        
    sorted_score_report = sorted(consensus_score_report.items(), key=lambda x: x[1][0])
    ssr_dict = {consensus: scores for consensus, scores in sorted_score_report}
    
    return consensus_groups, ssr_dict


# # STEP 1: NODE SELECTION #
# 
# ### Find reference pairs from the reference sequences. 
# 
# ### Then, find lmers similar to the reference pair from all other input sequences
# 
# ### If a reference pair can't find at least one similar lmer from each input sequence, its not valid
# 
# 

# In[20]:


def select_nodes(sequence_list, l, _2d, show_progress):
        
    num_seq = len(sequence_list)
        
    reference_seq_1 = sequence_list[0]
    reference_seq_2 = sequence_list[1]
    reference_pair_set = set()
    reference_pair_list = []
    reference_pair_counter = 0
        
    for i in (tqdm(range(len(reference_seq_1) - l + 1), desc='Finding reference pairs and sequence nodes', position=0, leave=True, file=sys.stdout) if show_progress else range(len(reference_seq_1) - l + 1)):
        for j in range(len(reference_seq_2) - l + 1):
            if hamming_distance(reference_seq_1[i:i+l], reference_seq_2[j:j+l]) <= _2d:
                ref_sub_seq = reference_seq_1[i:i+l]
                ref_sub_seq2 = reference_seq_2[j:j+l]
                if (ref_sub_seq, ref_sub_seq2) in reference_pair_set:
                    continue
                reference_pair = Reference_Pair(pair_id=reference_pair_counter, seq1=ref_sub_seq, seq2=ref_sub_seq2, pos1=i, pos2=j, num_seq=num_seq)    
                
                count_seqs_with_nodes = 0
                for k in range(2, num_seq):                
                    has_nodes = False
                    for m in range(len(sequence_list[k]) - l + 1):
                        sub_seq = sequence_list[k][m:m+l]
                        if hamming_distance(sub_seq, ref_sub_seq) <= _2d and hamming_distance(sub_seq, ref_sub_seq2) <= _2d:
                            seq_node = Sequence_Node(seq=sub_seq, id=k, pos=m)
                            reference_pair.node_sets[k-2].add_member(seq_node)
                            has_nodes = True
                    if has_nodes:
                        count_seqs_with_nodes += 1
                
                if count_seqs_with_nodes == num_seq - 2:
                    reference_pair_list.append(reference_pair)
                    reference_pair_set.add((ref_sub_seq, ref_sub_seq2))
                    reference_pair_counter += 1
                    
    return reference_pair_list
                

            




# # STEP 2: TREE CONSTRUCTION #

# In[21]:


def construct_trees(reference_pair_list, _2d, show_progress):
    
    num_seq = len(reference_pair_list[0].node_sets) + 2
    # tree_count = 0
    cliques = []
    for pair in (tqdm(reference_pair_list, desc='Building trees from reference pairs', position=0, leave=True, file=sys.stdout) if show_progress else reference_pair_list):
        node_sets = pair.node_sets
        for root in node_sets[0].members:
            root_node = Tree_Node(seq=root.seq, seq_id=root.id, pos=root.pos)
            tree = Tree(root_node)
            
            for i in range(3, num_seq):
                flag = False
                for node in node_sets[i - 2].members:
                    branches = extendable(node, tree, i, _2d)
                    if branches:
                        # print("branches check")
                        for branch in branches:
                            new_leaf = Tree_Node(seq=node.seq, seq_id=node.id, pos=node.pos)
                            branch.add_child(new_leaf)
                            tree.save_leaf(new_leaf)
                            tree.unsave_leaf(branch)
                        flag = True
                if not flag:
                    tree = None
                    break
                tree.height += 1
                prune_branches(tree)
            if tree and check_tree(tree):    
                pair.tree_list.append(tree)
                # tree.id = tree_count
                # tree_count += 1
                
                new_cliques = find_paths(tree.root, [], [])

                for new_clique in new_cliques:
                    
                    new_clique.nodes.insert(0, Tree_Node(pair.seq2, 1, pair.pos2))
                    new_clique.seqs.insert(0, pair.seq2)
                    
                    new_clique.nodes.insert(0, Tree_Node(pair.seq1, 0, pair.pos1))
                    new_clique.seqs.insert(0, pair.seq1)
                    
                    cliques.append(new_clique)
                    
                    
    return cliques


# <span style="font-size: 5em;">PyTreeMotif</span>
# <p></p>
# 
# ### Define the full algorithm as the function find_motif() ###

# In[22]:


def find_motif(DNA, motif_length, d=-1, print_output=True, show_progress=True, output_consensus_groups=False):
    
    DNA = [seq.upper() for seq in DNA]
    reference_pairs = None
    motifs = None
    _2d = 0
    
    if d == -1:
        i = 0
        while not reference_pairs and not motifs and i < 5:
            _2d = i * 2
            reference_pairs = select_nodes(DNA, motif_length, _2d, show_progress)
            if not reference_pairs:
                i += 1
                continue
            motifs = construct_trees(reference_pairs, _2d, show_progress)
            i += 1
        if not reference_pairs:
            return ["no reference pairs found :("]
    else:
        _2d = d * 2
        reference_pairs = select_nodes(DNA, motif_length, _2d, show_progress)
        if not reference_pairs:
            return ["no reference pairs found :("]
        motifs = construct_trees(reference_pairs, _2d, show_progress)
    

    if not motifs:
        return ["no motifs found :("]
    
    score_motifs(motifs, motif_length)
    score_sorted_motifs = sorted(motifs, key=lambda motif: motif.rank)
    
    best_motifs = []
    i = 0
    while i < len(score_sorted_motifs) and i < 3:
        best_motifs.append(score_sorted_motifs[i])
        i += 1
        
    consensus_groups, sorted_group_score_report = consensus_grouping(score_sorted_motifs)
    
    if print_output:
        
        dash_line = '\n' + ('-' * 80) + '\n'
        
        if output_consensus_groups:
            print("\n\n ### CONSENSUS GROUPS ### ")
            print(dash_line)
            for consensus, scores in sorted_group_score_report.items():
                print("Consensus: " + consensus)
                print("Number of cliques: " + str(len(consensus_groups[consensus])))
                print("Best score of this group: " + str(scores[0]))
                print("Average score of this group: " + str(scores[1]))
                print("Motifs:")
                for motif in consensus_groups[consensus]:
                    print("Rank: " + str(motif.rank) + ", Score: " + str(motif.score))
                    for instance in motif.nodes:
                        print("Seq: " + str(instance.seq_id) + ", Position: " + str(instance.pos) + ", " + str(instance.seq))
                    print()
                print(dash_line)

                
        print(dash_line)
        print("Total Cliques: " + str(len(motifs)))
        print("Total Consensus Motifs: " + str(len(consensus_groups)))
        if d == -1:
            print("Substitutions Allowed: " + str(_2d // 2))
        print(dash_line) 
        
        print("Top 3 Motifs")
        print(dash_line)
        for winner in best_motifs:
            print("Consensus: " + str(winner.consensus) + ", Rank: " + str(winner.rank) + ", Score: " + str(winner.score))
            print("Instances: ")
            for node in winner.nodes:
                print("Seq: " + str(node.seq_id) + ", Position: " + str(node.pos) + ", " + str(node.seq))
            print()
            
            
           
    
    return [motif.consensus for motif in best_motifs]



# # 15,4 BENCHMARKING #
# ### Test the algorithm with the original subtle motif challenge proposed by Pevzner ###
# 
# First define the function to generate a single dataset.
# 
# It will take a sequence length argument.
# 
# First, a 15mer motif is randomly generated.
# 
# For each of 20 sequences:
# 
# 1. The sequence is randomly generated at the given length
#     
# 2. The motif is mutated 4 times at random positions with random bases
#     
# 3. The mutated motif is implanted into the sequence at a random position
# 
# The dataset is now made up of 20 random sequences of the given length, each with a (15,4) implanted motif
# 

# In[23]:


def generate_15_4_dataset(seq_len, l, d):
    bases = ['A', 'C', 'G', 'T']
    
    motif = ''.join(random.choice(bases) for _ in range(l))
    
    # print("Implanted motif: " + motif)
    # print()
    
    dataset = []
    implanted_motifs = []
    
    for i in range(20):
        input_seq = ''.join(random.choice(bases) for _ in range(seq_len))
        
        sub_positions = [random.randint(0, l-1) for _ in range(d)]
        
        implanted_motif = motif
        
        for position in sub_positions:
            mutation = random.choice(bases)
            implanted_motif = implanted_motif[:position] + mutation + implanted_motif[position + 1:]
            
        implant_position = random.randint(0, seq_len - l)
        
        implanted_motifs.append((implanted_motif, implant_position))
            
        input_seq_w_motif = input_seq[:implant_position] + implanted_motif + input_seq[implant_position + l:]
        
        dataset.append(input_seq_w_motif)
    
    # for im in implanted_motifs:
    #     print(im[0], im[1])
        
    return dataset, motif, implanted_motifs
    


# # Evaluation Method
# 
# ### Replicate the evaluation scheme of Sun et al
# 
# Test 50 datasets of 20 600-bp long sequences for each of four (l,d) motifs:
# 
# * 12,3
# * 13,3
# * 14,4
# * 15,4
# 
# note: Sun et al's evaluation included additional l,d motifs that we do not include in our tests.  
# 
# The results are evaluated with nucleotide based metrics introduced by Pevzner and Sze (2000)  
# 
# The metrics, nucleotide performance coefficient (nPC) and nucleotide sensitivity (nSn) are built from a nucleotide level confusion matrix, where:  
# 
# * True Positives, TP is the number of nucleotide positions in both known sites and predicted sites
# * False Negatives, FN is the number of nucleotide positions in known sites but not in predicted sites
# * False Positive, FP is the number of nucleotide positions not in known sites but in predicted sites 
# * True Negative, TN is the number of nucleotide positions in neither known sites nor predicted sites
# 
# 
# &emsp;&emsp;nPC = TP/(TP+FN+FP)  
# &emsp;&emsp;nSn = TP/(TP+FN)  <p></p>
# 
# The best motif of the top 3 ranked motifs was used to obtain the best metric. Note that 0, 1, or 2 motifs is also a possible result.
# 
# These metrics were averaged over all 50 tests. Sun et al only used nSn in a comparison with other algorithms, so the only result available is for 15,4 motifs.
# 
# The following results were presented by Sun et al: <p></p>
# 
# | l,d  | nPC         | nSn         |
# |------|-------------|-------------|
# | 12,3 | 0.82 ± 0.07 | -           |
# | 13,3 | 0.94 ± 0.05 | -           |
# | 14,4 | 0.80 ± 0.19 | -           |
# | 15,4 | 0.95 ± 0.05 | 1.00 ± 0.01 |
# 
# Due to time constraints for submitting the project we did not run this evaluation in the notebook. 
# 
# Separate python scripts were generated to run each test in unison on NCSU's Hazel cluster, and metrics were calculated by hand based on the output shown in the cell above. <p></p>
# 
# For example, if one of the 50 tests had this result: 
# 
#     Test 3  
#     Implanted Motif: GAATCAGGGTGAA  
#     Predicted Motif 1: GAATCAGGGTGAA    MATCH  
#     Predicted Motif 2: GAATCAGGGTGAA    MATCH  
#     Predicted Motif 3: GAATCAGGGTGAA    MATCH  
# 
# or this result
# 
#     Test 7  
#     Implanted Motif: GTTGGCAGTCGGA  
#     Predicted Motif 1: GTTGGCAGTCGGA    MATCH  
#     Predicted Motif 2: TTGGCAGTCGGAG  
# 
# then nPC = 1. If the motif was correctly identified, there are no FNs or FPs. <p></p>
# 
# However, for this result:
# 
#     Test 50  
#     Implanted Motif: TCCATAGGAGATT  
#     Predicted Motif 1: CCATAGCAGATTT  
# 
# the motif was not precisely identified, resulting in nucleotide FNs and FPs.
#   
# In this case, the implanted motif and predicted motif are misaligned by 1 base, as well as another substitution in the middle of the motif
# 
# There is a FN before the beginning of the predicted motif, and a FP at the last base. 
# 
#  TCCATAGGAGATT  
#  <span style="color: red;">T</span>CCATAG<span style="color: orange;">C</span>AGATT<span style="color: orange;">T</span>
# 
# Therefore, nPC = 12/(12+1+2) = 0.8 <p></p>
# 
# Finally, for this result:
# 
#     Test 19  
#     Implanted Motif: ATTGTAATTCACT  
#     Predicted Motif 1: no motifs found :(  
# 
# nPC = 0.  <p></p><p></p>
# 
# For results with found motifs that aren't a perfect match, use the Needleman-Wunsch global alignment algorithm to align each predicted motif to the implanted motif. Then, we can count true positives, false positives, and false negatives - all of the terms needed to calculate nPC. 
# 
# Finally, get the mean and standard deviation of all 50 nPCs.
# 
# 

# In[38]:


def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_score=-1):
    # Initialize the matrix
    m, n = len(seq1), len(seq2)
    score_matrix = [[0 for _ in range(n+1)] for _ in range(m+1)]

    # Initialize first row and column
    for i in range(m+1):
        score_matrix[i][0] = i * gap_score
    for j in range(n+1):
        score_matrix[0][j] = j * gap_score

    # Fill the matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            if seq1[i-1] == seq2[j-1]:
                match = score_matrix[i-1][j-1] + match_score
            else:
                match = score_matrix[i-1][j-1] + mismatch_score
            delete = score_matrix[i-1][j] + gap_score
            insert = score_matrix[i][j-1] + gap_score
            score_matrix[i][j] = max(match, delete, insert)
            
    # for i in range(m+1):
    #     print(score_matrix[i])

    # Traceback to find alignment
    # we did not store traceback pointers
    align1, align2 = '', ''
    i, j = m, n
    while i > 0 and j > 0:
        score_current = score_matrix[i][j]
        score_diagonal = score_matrix[i-1][j-1]
        score_up = score_matrix[i][j-1]
        score_left = score_matrix[i-1][j]

        if score_current == score_diagonal + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score):
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif score_current == score_left + gap_score:
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1
        else: # score_current == score_up + gap_score
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1

    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq1[i-1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j-1]
        j -= 1

    # The sequences are aligned from the end to the start so we reverse them
    align1 = align1[::-1]
    align2 = align2[::-1]
    

    return align1, align2


# Now use it to calculate nPC for all 50 tests, and get the mean and standard deviation. 
# 
# The entire evaluation process for one l,d condition from start to finish is shown below. Run a fast example again with short sequences.

# In[35]:


# Run the algorithm

implanted_motif_list = []
found_motif_list = []
for i in tqdm(range(50), desc='Testing datasets', position=0, leave=True):
    dataset, implanted_motif, instances = generate_15_4_dataset(600, 12, 3)
    found_motif = find_motif(dataset, 12, 3, print_output=False, show_progress=False)
    
    implanted_motif_list.append(implanted_motif)
    found_motif_list.append(found_motif)


# Calculate and print results

print("15,4 Implanted Motif Tests\n")
npc_list = []
for i in range(50):
    print("Test " + str(i + 1))
    print("Implanted Motif: " + implanted_motif_list[i])
    
    
    best_npc = 0
    
        
    for j in range(len(found_motif_list[i])):
        curr_npc = 0
        match = ""
        
        if implanted_motif_list[i] == found_motif_list[i][j]:
            match = "\tMATCH"
            best_npc = 1
        elif found_motif_list[i][0][0] != 'n':
            # alignment
            aligned_IM, aligned_PM = needleman_wunsch(implanted_motif_list[i], found_motif_list[i][j])
            TPs, FPs, FNs = 0, 0, 0
            for k in range(len(aligned_PM)):
                if aligned_PM[k] == aligned_IM[k]:
                    TPs += 1
                elif aligned_PM[k] == '_':
                    FNs += 1
                else:
                    FPs += 1
            
            curr_npc = TPs / (TPs + FPs + FNs)
            best_npc = max(best_npc, curr_npc)
              
        print("Predicted Motif " + str(j + 1) + ": " + found_motif_list[i][j] + match)
    print("Best nPC: " + str(best_npc))
    npc_list.append(best_npc)
    print()
    
npc_numpy = np.array(npc_list)
npc_mean = np.mean(npc_numpy, axis=0)
npc_std = np.std(npc_numpy, axis=0)

print('\n' + ('-' * 80) + '\n')
print("Mean nPC: " + str(npc_mean))
print("Std nPC: " + str(npc_std))

    




