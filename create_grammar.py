#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import itertools
import pickle
import os
import copy

import tree


class StoreFormula:
    def __init__(self, nodes, quantifiers, cost, time, length):
        
        self.nodes = nodes
        self.quantifiers = quantifiers
        self.cost = cost
        self.time = time
        self.length = length

def list_tuples_sum(num_addend,value_sum):
    
    ''' The function computes all possible tuples of num_addend elements whose sum is equal to value_sum
    INPUTS: 
        - num_addend : number of elements that are added
        - value_sum : result of the sum
    '''
    
    tuples = [] #1,value_sum-num_addend
    
    for comb in itertools.product(range(1, value_sum-num_addend+2), repeat =num_addend):
        if sum(comb) == value_sum: tuples.append(comb)
    return tuples


# OSTER OSTERVIZ next stop st veiz an glan

def build_formulas(length, spec, top_level_spec, nodes, index_parent, kind_child): 
    '''Given a grammar, the function constructs all the formulas whose length is length 
    generated starting from the top_level_spec (number from 0 to levels) and using grammar.
    
    INPUTS:
        - length: length of the desired formulas
        - spec: grammar_structure. Grammar generating the formulas
        - top_level_spec: the index of the specification in the grammar we are reffering to
        - nodes
        - index_parent : the index in the nodes list that is the parent of the current node
        - kind_child: could be either 'left' or 'right'; it indicates is a left or right child of its parent '''
    
    current_spec =  spec[top_level_spec]
    
    for index in range(len(current_spec)):
        if length >= current_spec[index].count('phi') + 1: #If the length for an option is less than the minimal length of that option
            if 'predicate' in current_spec[index]: 
                if length == 1: 
                     nodes.append(tree.BinaryTreeNode(current_spec[index])) #'predicate'
                     nodes[-1].parent = nodes[index_parent]
                     nodes[-1].level_grammar = top_level_spec 
                    
                     if kind_child == 'right': nodes[index_parent].rightchild_index.append(len(nodes)-1)
                     elif kind_child == 'left': nodes[index_parent].leftchild_index.append(len(nodes)-1)
                     
                     return nodes
            
            #unary
                
            elif 'not' in current_spec[index] or \
                 'always' in current_spec[index] or \
                 'eventually' in current_spec[index] or \
                 's_next' in current_spec[index]  or \
                 'next' in current_spec[index] or \
                 's_prev' in current_spec[index] or \
                 'prev' in current_spec[index] or \
                 'historically' in current_spec[index] or \
                 'once' in current_spec[index]:
                     
                next_spec = int(current_spec[index][current_spec[index].index('phi')+3])
                
                current_spec_copy = current_spec[index]
                current_spec_copy = current_spec_copy.replace(f' phi{next_spec}','')
            
                nodes.append(tree.BinaryTreeNode(current_spec_copy))
                nodes[-1].level_grammar = top_level_spec
                len_nodes = len(nodes)
                nodes[len_nodes -1].parent = nodes[index_parent]
                if kind_child == 'right': nodes[index_parent].rightchild_index.append(len(nodes)-1)
                elif kind_child == 'left': nodes[index_parent].leftchild_index.append(len(nodes)-1)
                     
                nodes = build_formulas(length-1, spec, next_spec, nodes,len_nodes -1, 'left')
            
            
            else: #binary operators
                
                tuples_of_length = list_tuples_sum(2,length-1) # all possibilities
                next_spec1 = int(current_spec[index][current_spec[index].index('phi')+3])
                next_spec2 = int(current_spec[index][current_spec[index].index('phi', current_spec[index].index('phi')+3)+3])
                
                for tuples in tuples_of_length:
                    
                    if 'equiv' in current_spec[index] or \
                       'implies' in current_spec[index] or \
                        'until' in current_spec[index] or \
                        'since' in current_spec[index] or \
                        'weak_u' in current_spec[index] or\
                        'or' in current_spec[index] or \
                         'and' in current_spec[index]: 
                          
                        current_spec_copy = current_spec[index]
                        current_spec_copy = current_spec_copy.replace(f'phi{next_spec1} ','')
                        current_spec_copy = current_spec_copy.replace(f' phi{next_spec2}','')
                        
                        nodes.append(tree.BinaryTreeNode(current_spec_copy))
                        nodes[-1].level_grammar = top_level_spec
                        len_nodes = len(nodes)
                        nodes[len_nodes -1].parent = nodes[index_parent]
                        if kind_child == 'left': nodes[index_parent].leftchild_index.append(len(nodes)-1)
                        elif kind_child == 'right': nodes[index_parent].rightchild_index.append(len(nodes)-1)
                        nodes = build_formulas(tuples[0], spec, next_spec1, nodes,len_nodes -1, 'left')
                        nodes = build_formulas(tuples[1], spec, next_spec2, nodes,len_nodes -1, 'right')
                     
                    else:
                        print('Used a non-existing operator: \n\nADD NEW OPERATOR!!\n\n')
             
    return nodes                
    

def transform_listing_to_nodes(listing_nodes, length, spec_level, m):
    
    if not os.path.exists('grammar'): os.mkdir('grammar')
    if not os.path.exists(f'grammar/spec_level{spec_level}_numbquant{m}'): os.mkdir(f'grammar/spec_level{spec_level}_numbquant{m}')
    
    database = []
    #Indices of starting nodes (i.e., nodes that are the leftchildren of the node 'start')
    for start_index in listing_nodes[0].leftchild_index:
        
        current_formula = [tree.BinaryTreeNode(listing_nodes[start_index].data)]
        current_formula[-1].level_grammar = spec_level
        current_index = start_index
        
        path = f'grammar/spec_level{spec_level}_numbquant{m}/length_{length}' 
        
        if not os.path.exists(path): os.mkdir(path)
        new_database = []
        new_database, _ = creation_node(current_index, listing_nodes, 0 , current_formula, new_database, length)
        
        for _, new_formula in enumerate(new_database):
            is_new = True
            for _, existing_formula in enumerate(database):
                
                string_new = tree.tree_to_formula(new_formula)
                string_old = tree.tree_to_formula(existing_formula)
                
                if string_old == string_new:
                    is_new = False
                    break
                
            if is_new:
                with open(f'{path}/nodes_{len(database)}.obj', 'wb') as f:
                    pickle.dump(new_formula, f)
                    database.append(new_formula)
         
    return 

def creation_node(current_index_ls, listing_nodes, current_index_cf, current_formula, database,  length):
    
    '''
    
    GOES DOWN IN THE TREE: FROM ROOT TO LEAF
    
    The function transforms the list of nodes created by the build formula 
    into proper nodes usable by the sytnax guided algorithm. 
    In particular, listing_nodes is build such that a node can have several left or right children
    so this function splits them into different functions.
    
    The output is saved into a file having given path
    
    INPUTS: 
        - current_index_ls: current index in list listing_nodes
        - listing_nodes: original nodes to transform
        - current_index_cf: current index in list current_formula
        - current_formula: formula constructed
        - path: path where to save the nodes
        - nb_formula: counter of the number of formulas saved so far
        - length of the target formula
    
    
    E.g. implies may have two left children and one right child.
    The generated formulas are two: 1) first left child + implies + the right child
                                    2) second left child + implies + the right child
    '''  

    
    if current_formula is None: 
        return
    
    #If one single left child exists: add it and iterate
    if  len(listing_nodes[current_index_ls].leftchild_index)==1:
        current_formula.append(tree.BinaryTreeNode(listing_nodes[listing_nodes[current_index_ls].leftchild_index[0]].data))
        current_formula[-1].level_grammar = listing_nodes[listing_nodes[current_index_ls].leftchild_index[0]].level_grammar
        current_formula[-1].parent = current_formula[current_index_cf]
        current_formula[current_index_cf].leftchild = current_formula[-1]
        
        #If the added node is a predicate, end the recursion
        # if 'predicate' in current_formula[-1].data : return database, current_formula
        if 'predicate' not in current_formula[-1].data: 
            [database, current_formula] = creation_node(listing_nodes[current_index_ls].leftchild_index[0],\
                              listing_nodes, len(current_formula)-1, current_formula, database, length)
    
    #If one single right child exists: add it and iterate
    if len(listing_nodes[current_index_ls].rightchild_index)==1:
        current_formula.append(tree.BinaryTreeNode(listing_nodes[listing_nodes[current_index_ls].rightchild_index[0]].data))
        current_formula[-1].level_grammar = listing_nodes[listing_nodes[current_index_ls].rightchild_index[0]].level_grammar
        current_formula[-1].parent = current_formula[current_index_cf]
        current_formula[current_index_cf].rightchild = current_formula[-1]
        
        #If the added node is a predicate, end the recursion
        # if 'predicate' in current_formula[-1].data: return database, current_formula
        if 'predicate' not in current_formula[-1].data: 
            [database,current_formula] =creation_node(listing_nodes[current_index_ls].rightchild_index[0],\
                              listing_nodes,len(current_formula)-1, current_formula, database, length)
        
    #If it is an unsuccessful node (Because it ends but is its not a predicate)
    if  (len(listing_nodes[current_index_ls].rightchild_index)==0 \
            and len(listing_nodes[current_index_ls].leftchild_index)==0 \
            and 'predicate' not in listing_nodes[current_index_ls].data): return database, None
     
        
    #If there is a fork with several left children
    if len(listing_nodes[current_index_ls].leftchild_index)>1:
       for left_children_index in listing_nodes[current_index_ls].leftchild_index:
            new_formula = copy.deepcopy(current_formula)
            
            new_formula.append(tree.BinaryTreeNode(listing_nodes[left_children_index].data))
            new_formula[-1].level_grammar = listing_nodes[left_children_index].level_grammar
            new_formula[-1].parent = new_formula[current_index_cf]
            new_formula[current_index_cf].leftchild = new_formula[-1]
            
            [database,new_formula] = creation_node(left_children_index, listing_nodes, \
                          len(new_formula)-1,new_formula, database, length)  
            
            ##Climb back the tree to explore the right children of the nodes    
            index_to_use_ls = current_index_ls #index in the global list of nodes
            index_to_use_cf = current_index_cf #index in the local list of nodes (current formula)
            backup_new_formula = copy.deepcopy(new_formula)
            
            if new_formula is not None and len(new_formula) < length:
                
                # while True:#Climb back to the root
                #     #Loop to find a node with some right children
                    while True:
                        if index_to_use_cf == 0 or len(listing_nodes[index_to_use_ls].rightchild_index) > 0: break  # right children
                        else :#climb to a parent with some right children
                            index_to_use_cf = current_formula.index(current_formula[index_to_use_cf].parent)
                            index_to_use_ls = listing_nodes.index(listing_nodes[index_to_use_ls].parent)
                    
                    #Loop over the right children    
                    for right_children_index in listing_nodes[index_to_use_ls].rightchild_index:
                        new_formula = copy.deepcopy(backup_new_formula)
                        
                        new_formula.append(tree.BinaryTreeNode(listing_nodes[right_children_index].data))
                        new_formula[-1].level_grammar = listing_nodes[right_children_index].level_grammar
                        new_formula[-1].parent = new_formula[index_to_use_cf]
                        new_formula[index_to_use_cf].rightchild = new_formula[-1]
                        
                        [database, new_formula] = creation_node(right_children_index, listing_nodes, \
                               len(new_formula)-1,new_formula, database, length)
                    
                    # if index_to_use_cf == 0: break
           
                
    #If there is a fork with several right children
    if len(listing_nodes[current_index_ls].rightchild_index)>1:
        
        for right_children_index in listing_nodes[current_index_ls].rightchild_index:
            new_formula = copy.deepcopy(current_formula)
            
            new_formula.append(tree.BinaryTreeNode(listing_nodes[right_children_index].data))
            new_formula[-1].level_grammar = listing_nodes[right_children_index].level_grammar
            new_formula[-1].parent = new_formula[current_index_cf]
            new_formula[current_index_cf].rightchild = new_formula[-1]
            
            [database, new_formula] = creation_node(right_children_index, listing_nodes, \
                           len(new_formula)-1,new_formula, database, length)
           
    if current_formula is not None and len(current_formula) == length and sanity_check(current_formula):
        database.append(current_formula)
        # tree.print_tree(current_formula) 
    return database, current_formula




def sanity_check(formula):
    
    '''The function checks whether the formula is well constructed: i.e., 
     - all binary operators have two children,
     - all unary operator have one child,
     - no double negation is present'''
    if formula is not None: 
        for index, node in enumerate(formula):
            #Binary
            if       'or' in node.data \
                or 'implies'  in node.data \
                or 'equiv' in node.data \
                or 'and ' in node.data \
                or  'until' in node.data\
                or  'since' in node.data\
                or 'weak_u' in node.data:
                    
                 if node.leftchild is None or node.rightchild is None: return False
            #Unary
            elif 'not' in node.data\
                or 's_next' in node.data\
                or 'next' in node.data\
                or 's_prev' in node.data\
                or 'prev' in node.data\
                or  'eventually' in node.data\
                or  'always' in node.data\
                or  'once'in node.data\
                or  'historically'in node.data:
                    
                 if node.leftchild is None: return False
                 if 'not' in node.data and 'not' in node.leftchild.data: return False
        return True
    return False
 