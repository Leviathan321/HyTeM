#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# import itertools
import random
import copy

class BinaryTreeNode:
    def __init__(self, data):
        self.data = data
        self.leftchild = None
        self.rightchild = None
        self.parent = None 
        
        # The next attributes are needed to construct formulas when the grammar is given by the user
        self.leftchild_index = []
        self.rightchild_index = []
        
        self.level_grammar = None
        self.polarity = None

def sample_predicate(grammar_predicate, nodes, bool_different_variables, bool_same_number, second_variable):
   
    '''The function is used to sample predicates in case the user specifies if some 
    predicates need to be equal or they need to be different by imposing for examples:
        
        ['(predicate1 Until predicate2) or (always predicate1)'] for the weak until
        
    INPUTS:
        - grammar_predicate: list of predicates levels = [ predicate0, predicate1, predicate2, ... ]
                             where predicatei = [[rel_i_0, var_i_0, pair_traces_i_0], [rel_i_1, var_i_1, pair_traces_i_1],...], where
                                  rel : relational signs admitted
                                  var : admitted variables
                                  pair_traces: pair of traces admitted
                                                         
            
        - nodes: nodes (tree-structured) representing the formula;
            
        - bool_different_variables: True if a predicate can involve different variables (e.g. pi(x) = pi(y));
                
        -bool_same_number: True if whenever the same predicate name is used (e.g. 'predicate2'), all
                            nodes with the same name has to be equal;
                            BUT the presence of two different names (e.g. 'predicate1' and 'predicate2')
                            does not necessarily implies they are different.
        - second_variable: 'true' if the second variable should be True,
                           'actions' if the second variable should be a number representing a certain action (see dining philophers)
                           'variable' for any other variable
        
    '''
    
    number_items = []
   
    #number_items stores the number of concrete predicates for each level in grammar_predicate
    for level in range(len(grammar_predicate)):
        aux = []
        for sub_level in range(len(grammar_predicate[level])):
            aux.append(len(grammar_predicate[level][sub_level][0])*len(grammar_predicate[level][sub_level][1])*len(grammar_predicate[level][sub_level][2]))
        number_items.append(aux)
    
    list_bool = [False]*len(nodes) # if it is already apperead
    list_values = ['1']*len(nodes) 
    
    for index_node, node in enumerate(nodes): # Loop over the different nodes to modify all the strings 'predicate..'
        for i in range(len(grammar_predicate)):#Loop over grammar levels to find the number next to 'predicate' (e.g. 'predicate1')
            
          if f'predicate{i}' in node.data: #If the node is a predicate
                  
               #If assigning same predicate to same name and the current predicate name has already been assigned to a concrete predicate
               if  bool_same_number and list_bool[index_node] == True: node.data = list_values[index_node] #just copy a previous predicate
               
               else: #(i.e., bool_same_number == False or (bool_same_number == True and list_bool[index_node] == False))
                       
                     level_index = i 
                     
                     #Sample the sub-level in grammar_predicate[level_index] according to the number of predicates in each sub-level             
                     if len(grammar_predicate[level_index]) > 1:
                            random_num = random.randrange(sum(number_items[level_index]))#sample a number in the total number of predicates (in that level)
                            
                            for j in range(1, len(grammar_predicate[level_index])+1):
                                if random_num < sum(number_items[level_index][:j]): 
                                    sub_index = j - 1 #grammar level
                                    break
                    
                     else: sub_index = 0 #grammar sub-level
                    
                     # Sample the sign symbol in the grammar level '[level_index]' and sub-level 'sub_index'
                     sign = random.sample(grammar_predicate[level_index][sub_index][0], k=1)[0]
                     
                     # Sample the first variable
                     variable1 = random.sample(grammar_predicate[level_index][sub_index][1], k=1)[0]
                     ##If different variables are admitted in the predicates: 
                     if bool_different_variables == True: ##Different variables admitted in the predicates
                         variable2 = random.sample(grammar_predicate[level_index][sub_index][1], k=1)[0]
                     else: variable2 = variable1
                     
                     #Sample the trace variable
                     trace_variables = random.sample(grammar_predicate[level_index][sub_index][2], k=1)[0]
                     
                     
                     ##Predicate in the form of formula string
                     
                     if  second_variable == 'true': node.data = f'pi[{trace_variables[0]}][{variable1}] {sign} True '
                     elif second_variable == 'variable': 
                         if len(trace_variables) < 2: node.data = f'pi[{trace_variables[0]}][{variable1}] {sign} pi[{trace_variables[0]}][{variable2}] '
                         else: node.data = f'pi[{trace_variables[0]}][{variable1}] {sign} pi[{trace_variables[1]}][{variable2}] '
                     else:
                           action = random.sample(second_variable, k=1)[0]
                           node.data = f'pi[{trace_variables[0]}][{variable1}] {sign} {action} '
                    
                     
                     list_values[index_node] = node.data
                     list_bool[index_node] = True
            
    return nodes

def sample_predicate_STL(grammar_predicate, nodes, bool_different_variables):
    
    '''The function samples among predicates defined by the user. This is used for STL formulas:
        the user encodes a finite set of predicate templates in the grammar
        and the method can choose only among those.
        
        
        INPUTS:
            - grammar_predicate: in this case, it is a list of possible predicates 
            - nodes: nodes of the formula, where the 'predicate' need to be replaced by a sampled concrete predicate
            - bool_different_variables: NOT USED 
        
        
    '''
    
    #Count the number of predicates in the formula
    numb_predicates = 0
    for node in nodes:
        if 'predicate' in node.data:  numb_predicates +=1
        
    if numb_predicates <= len(grammar_predicate):
        #Sample numb_predicates different predicates (WITHOUT REPLACEMENT)
        sampled_predicates = random.sample(grammar_predicate, k = numb_predicates)
    else:
         #Sample numb_predicates WITH REPLACEMENT
         sampled_predicates = random.choices(grammar_predicate, k=numb_predicates)
    
    numb_predicates = 0
    #Replace the predicates with the sample concrete predicates
    for node in nodes:
         if 'predicate' in node.data:
             
             #Predicate in the form of formula string
             node.data = sampled_predicates[numb_predicates]
             numb_predicates +=1
    
    return nodes
    
    

def sample_quantifiers(grammar_quantifiers, min_nb_quantifiers, max_nb_quantifiers):
    
    '''Standard: grammar_quantifiers = [ ['E psi0', 'A psi0', 'phi0' ]]'''
    
    # For quantifiers: 
    # symbol 0 corresponds to EXISTS
    # symbol 1 corresponds to FOR ALL
    
    if min_nb_quantifiers is None: min_nb_quantifiers = 2
    if max_nb_quantifiers is None: max_nb_quantifiers = 10
    
        
    list_possible_strings = []
    current_formula = [ 0 ]
    level_grammar = 0
        
    list_possible_formulas = build_quantifiers(grammar_quantifiers, level_grammar, current_formula, list_possible_strings, min_nb_quantifiers, max_nb_quantifiers)
    
    if len(list_possible_formulas)==0: print('\nIt was not possible to find a string of quantifiers with admissible length!\n')
    
    word_quantifiers = random.sample(list_possible_formulas, k=1)[0]
    
    quantifiers = []
    #The loop starts from 1, because at index 0 there is only the number of quantifiers of the formula
    for i in range(1, len(word_quantifiers)): 
        for item in  word_quantifiers[i][0]: #contains the quantifiers
            if item == 'E': quantifiers.append(0)
            elif item == 'A': quantifiers.append(1)
        
    spec_level = int(word_quantifiers[-1][0][-1]) # the level in the spec grammar where the selected quantifiers point to
    
    # word_quantifiers[1:] is outputted because it will be used when changing the quantifiers; it indicates the level in the
    # quantifiers grammar for each quantifier
    
    #spec_level indicates the level in the grammar_spec that has to be called when using this quantifier 
    return quantifiers, word_quantifiers, spec_level

def build_quantifiers(grammar_quantifiers, level_grammar, current_formula, list_possible_formulas, min_nb_quantifiers, max_nb_quantifiers):
    
    '''Given a grammar, the function constructs all the strings of quantifiers that derive from the grammar 
    and whose length is within the given limits. 
    
    INPUTS:
        - grammar_quantifiers : the grammar that defines the quantifiers given as a list. 
            The upper level is defined by the the first element of the list. 
            If the all Hyper-LTL is considered: grammar_quantifiers = [ ['E psi0', 'A psi0', 'phi0' ]]
            
        - level_grammar : the level of the grammar, indicated as index of the list grammar_quantifiers
        
        - current_formula: the current quantifiers in a list of pairs the form: [[quantifier, quantifier_level], [], []]
        
        - list_possible_formulas: list of admissible formulas find so far having the form of list of current_formulas
        
        - min_nb_quantifiers: minimum number of quantifiers 
        
        - max_nb_quantifiers: maximum number of quantifiers   
    '''
   
    #For all possible options in the current grammar level
    for index, item in enumerate( grammar_quantifiers[level_grammar]): 
        
        current_quantifiers = current_formula.copy()
        
        if 'phi' in item: #phi indicates the END of the recursion over the quantifiers
        
            #If the number of quantifiers is acceptable
            if (min_nb_quantifiers <= item.count('A')+ item.count('E') + current_quantifiers[0] <= max_nb_quantifiers): #ok to consider this string as acceptable
                # spec_level = int(item[(item.index('phi')+3)])#level of specification in grammar_structure that is required by these quantifiers
                current_quantifiers[0] += item.count('A')+ item.count('E')
                current_quantifiers.append([item, level_grammar])
                list_possible_formulas.append(current_quantifiers)
                # print(current_quantifiers)
                
        else:      
            
            if item.count('A')+ item.count('E') + current_quantifiers[0] <= max_nb_quantifiers: # try another iteration on the quantifier
                
                current_quantifiers[0] += item.count('A')+ item.count('E')
                spec_level = int(item[(item.index('psi')+3)])
                current_quantifiers.append([item[:-4], level_grammar]) # -4 because we keep only 'A' and 'E' from item (= eliminate psi* ) that corresponds to 4 characters
                list_possible_formulas = build_quantifiers(grammar_quantifiers,spec_level , 
                   current_quantifiers,
                  list_possible_formulas,min_nb_quantifiers, 
                   max_nb_quantifiers )
        
    return list_possible_formulas
                



def sample_new_quantifiers(quantifiers, quantifiers_levels, grammar_quantifiers):
    
    '''The function samples new quantifiers from the given ones, such that:
        
        - the final number of new_quantifiers is the same as the original ones;
        - the final new_quantifiers point at the same grammar level as the original quantifiers '''
    
    list_possible_formulas = []
    current_formula = [0 ]
    
    selected_index = random.randrange(1, len(quantifiers_levels)) #selet an index from where to apply changes 
    level_grammar = quantifiers_levels[selected_index][1] # the second component stores the quantifier level of the current quantifier
    
    # Take into account cases in qwhich one elements contain multilple quantifiers [[ 'A',0] , ['EE',0]]
    length_new_quant = 0 #length of the new string of quantifiers
    
    for i in range(selected_index, len(quantifiers_levels)):
        length_new_quant +=  (quantifiers_levels[i][0].count('A') + quantifiers_levels[i][0].count('E'))
    
    list_possible_formulas = build_quantifiers(grammar_quantifiers, level_grammar, 
            current_formula, list_possible_formulas, length_new_quant, length_new_quant)
    
    
    while True: 
        word_quantifiers = random.sample(list_possible_formulas, k=1)[0]
        #If the new quantifiers point to the same spec level 
        if int(word_quantifiers[-1][0][-1]) == int(quantifiers_levels[-1][0][-1]): 
            break
            
    
    #Generation of the new quantifiers and new quantifier levels 
      
    new_quantifiers = quantifiers[:selected_index - 1]
    new_quantifiers_levels = quantifiers_levels[:selected_index] 
    
    
    for i in range(1, len(word_quantifiers)): 
        new_quantifiers_levels.append(word_quantifiers[i])
        for item in  word_quantifiers[i][0]: #contains the quantifiers
            if item == 'E': new_quantifiers.append(0)
            elif item == 'A': new_quantifiers.append(1)
         
    return new_quantifiers, new_quantifiers_levels


# level_grammar = 0  
# grammar_quantifiers = [ ['E psi0', 'A psi0', 'phi0' ,'E E phi1'],  ['A A phi1']]
# min_nb_quantifiers = 2
# max_nb_quantifiers = 3
# quantifiers, quantifiers_levels, spec_level = sample_quantifiers(grammar_quantifiers, min_nb_quantifiers, max_nb_quantifiers)
    
# new_quantifiers, new_quantifiers_level = sample_new_quantifiers(quantifiers, quantifiers_levels, grammar_quantifiers)  

def check_double_negation(nodes):
    
    '''The function checks if in nodes, one 'not' node has as parent a not; 
    in that case, the output is True, otherwise it is False.'''

    
    for i in range(len(nodes)):
        if 'not' in nodes[i].data and nodes[i].parent is not None and 'not' in nodes[i].parent.data: 
            return True
       
    return False

def print_tree(tree_formula):
    for i, node in enumerate(tree_formula):
        print(f'\n\n\nindex={i}, node={node.data}')
        if node.rightchild is not None: print(f'rightchild={tree_formula.index(node.rightchild)}')
        if node.leftchild is not None: print(f'leftchild={tree_formula.index(node.leftchild)}')
        if node.parent is not None: 
            print(f'parent={tree_formula.index(node.parent)}')
        # print(f'rightchild_index={node.rightchild_index}')
        # print(f'leftchild_index={(node.leftchild_index)}')
       
            
def tree_to_formula(nodes):
    
    ## HYPOTHESIS:
    ## nodes[0] = root
    
    start = '' 
    formula = f'{translate_node(start, nodes,  0 )}'
    
    return formula
    
def translate_node(formula, nodes, i):
    
    '''Implication and biimplication here are expressed as combination of or/not so that they can be evaluated by the eval function-- 
      because not having the time, RTAMT is not applicable here'''
    
    if 'or' in nodes[i].data or 'and' in nodes[i].data:
        return f'{formula}( {translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) {nodes[i].data} ( {translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} )'
    
    elif 'not' in nodes[i].data or \
         'next'in nodes[i].data or\
         's_next'in nodes[i].data or\
         'prev'in nodes[i].data or\
         's_prev' in nodes[i].data:
                 return  f'{formula}({nodes[i].data}( {translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} ))'
   
    elif  'historically' in nodes[i].data \
          or 'once' in nodes[i].data \
          or 'eventually' in nodes[i].data \
          or  'always'  in nodes[i].data:
        return  f'{formula}({nodes[i].data}({translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} ))'
    
    elif 'implies' in nodes[i].data:
          return f'{formula}( ( {translate_node(formula,nodes, nodes.index(nodes[i].rightchild))}) or (not({translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} )) )'
    
    elif  'equiv' in nodes[i].data:
          return f'{formula}( ( ({translate_node(formula,nodes, nodes.index(nodes[i].rightchild))}) and({translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) ) or ( (not({translate_node(formula,nodes, nodes.index(nodes[i].leftchild))})) and (not({translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} ) )) )'
      # return f'{formula}( ( {translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} or not({translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) ) and ( {translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} or not({translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} ) ) )'
     
    elif  'until' in nodes[i].data \
            or 'weak_u' in nodes[i].data \
            or 'since' in nodes[i].data :
        return f'{formula}( ( {translate_node(formula,nodes, nodes.index(nodes[i].leftchild))}){nodes[i].data}({translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} ) )'
 
    else: #PREDICATE
        return nodes[i].data
    
    
    
def NATURAL_tree_to_formula(nodes):
    
    ## HYPOTHESIS:
    ## nodes[0] = root
    
    start = '' 
    formula = f'{NATURAL_translate_node(start, nodes,  0 )}'
    
    return formula
    
def NATURAL_translate_node(formula, nodes, i):
    
    '''As function translate_node but with implies and biimplication explicit (not expressed as combination of or/not)'''
    

    if 'or' in nodes[i].data or 'and' in nodes[i].data:
        return f'{formula}( {NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) {nodes[i].data} ( {NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} )'
    
    elif 'not' in nodes[i].data or\
         'next' in nodes[i].data or\
         's_next' in nodes[i].data or\
         'prev' in nodes[i].data or\
         's_prev' in nodes[i].data:
             return  f'{formula}{nodes[i].data}( {NATURAL_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} )'
    
    elif  'historically' in nodes[i].data\
            or 'once' in nodes[i].data\
            or 'eventually' in nodes[i].data\
            or 'always' in nodes[i].data:
               
        return  f'{formula}{nodes[i].data}( {NATURAL_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} )'
    
    elif 'implies' in nodes[i].data:
          return f'{formula}( {NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) --> ( {NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} )'
    
    elif 'equiv' in nodes[i].data:
          return f'{formula}( {NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) <--> ( {NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} )'
  
    elif  'until' in nodes[i].data\
          or 'weak_u' in nodes[i].data \
            or 'since' in nodes[i].data:
          return f'{formula}( {NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) {nodes[i].data} ( {NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} )'
  
    else: #PREDICATE
        return nodes[i].data
    


def select_node(nodes):
    
    '''The function uniformly randomly selects a node in a formula.
    
        INPUT: 
            - nodes: formula expressed as a tree (list of nodes)
            
        OUTPUTS:
            - indices_old_node: indices of the node (and all its descendants) selected that will be replaced
            - length_node: length of the chain of descendants of the selected node
    '''
    ## !? 
    # ROOT cannot be selected--> start from 1
    # ROOT can be selected--> start from 0
    if len(nodes ) == 1: selected_index = 0
    else: selected_index = random.randrange(1,len(nodes))
    #old_node stores the indices of the nodes that will be replace 
    indices_old_node = [selected_index]
    #length of the chain of descendants of nodes[selected_index]
    length_node = 1
     
    to_be_examined = [selected_index]
    
    while True:
        if len(to_be_examined)==0: break
        
        if nodes[to_be_examined[0]].leftchild is not None:
            indices_old_node.append(nodes.index(nodes[to_be_examined[0]].leftchild))
            to_be_examined.append(nodes.index(nodes[to_be_examined[0]].leftchild))
            length_node +=1
        
        if nodes[to_be_examined[0]].rightchild is not None:
            indices_old_node.append(nodes.index(nodes[to_be_examined[0]].rightchild))
            to_be_examined.append(nodes.index(nodes[to_be_examined[0]].rightchild))
            length_node +=1
            
        to_be_examined.pop(0)
        
    return indices_old_node, length_node



def change_node(nodes, indices_old_node, new_node):
    
    '''The function replaces the nodes corresponding to indices_old_nodes in
    original_nodes with new_node'''
    
    right_child_bool = False
    
    is_root = True
    
    new_nodes = copy.deepcopy(nodes)
    
    ## REMOVE the bound from the selected node and its parent
    #if the node selected to be removed (together with all its descendants), it is not the root: 
    if new_nodes[indices_old_node[0]].parent is not None:
        # if the node randomly selected was a right child
        if new_nodes[indices_old_node[0]].parent.rightchild == new_nodes[indices_old_node[0]]:
            #remove the bound: the selected node is no more the right child of its parent
            new_nodes[indices_old_node[0]].parent.rightchild = None
            right_child_bool = True
            
        # if the node randomly selected was a left child
        elif new_nodes[indices_old_node[0]].parent.leftchild == new_nodes[indices_old_node[0]]:
            #remove the bound: the selected node is no more the left child of its parent
            new_nodes[indices_old_node[0]].parent.leftchild = None
            
        new_parent = new_nodes.index(new_nodes[indices_old_node[0]].parent)
        is_root = False
            
    ## REMOVE all elements having index in indices_old_node from original_nodes
    indices_old_node.sort(reverse=True)
    for item in indices_old_node: new_nodes.pop(item)
    
    ## ADD elements in new_node to original_nodes
    for i, item in enumerate(new_node): 
        new_nodes.append(item)
        if i == 0 and (is_root==False): 
            if right_child_bool: new_nodes[new_parent].rightchild = new_nodes[-1]
            else: new_nodes[new_parent].leftchild = new_nodes[-1]
            new_nodes[-1].parent = new_nodes[new_parent]
    
    return new_nodes


# ##OLD FORMAT:

# def OLD_tree_to_formula(nodes):
    
#     ## HYPOTHESIS:
#     ## nodes[0] = root
    
#     start = '' 
#     formula = f'{OLD_translate_node(start, nodes,  0 )}'
    
#     return formula
    
# def OLD_translate_node(formula, nodes, i):
    
#     '''Implication and biimplication here are expressed as combination of or/not so that they can be evaluated by the eval function-- 
#      because not having the time, RTAMT is not applicable here'''
    
#     if nodes[i].data == 'or':
#         return f'{formula}( {OLD_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) or ( {OLD_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} )'
    
#     elif nodes[i].data == 'and':
#         return f'{formula}( {OLD_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) and ( {OLD_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} )'
      
#     elif nodes[i].data == 'not':
#         return  f'{formula}(not( {OLD_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} ))'
#     #next
#     elif nodes[i].data == 'next':
#         return  f'{formula}(next( {OLD_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} ))'
#     #strong next
#     elif nodes[i].data == 's_next':
#         return  f'{formula}(s_next( {OLD_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} ))'
    
#     elif nodes[i].data == 'historically':
#         return  f'{formula}(historically( {OLD_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} ))'
    
#     elif nodes[i].data == 'once':
#         return  f'{formula}(once( {OLD_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} ))'
   
#     elif nodes[i].data == 'prev':
#         return  f'{formula}(prev( {OLD_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} ))'
#     elif nodes[i].data == 's_prev':
#         return  f'{formula}(s_prev( {OLD_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} ))'
    
    
#     elif nodes[i].data == '->':
#          return f'{formula}( ( {OLD_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))}) or (not({OLD_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} )) )'
    
#     elif nodes[i].data == '<->':
#          return f'{formula}( ( ({OLD_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))}) and({OLD_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) ) or ( (not({OLD_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))})) and (not({OLD_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} ) )) )'
#       # return f'{formula}( ( {translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} or not({translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) ) and ( {translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} or not({translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} ) ) )'
      
#     #Time operators 
#     elif nodes[i].data == 'eventually':
#         return  f'{formula}(eventually( {OLD_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} ))'
    
#     elif nodes[i].data == 'always':
#         return  f'{formula}(always( {OLD_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} ))'
    
#     elif nodes[i].data == 'until':
#         return f'{formula}( ( {OLD_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))}) until ({OLD_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} ) )'
    
#     elif nodes[i].data == 'weak_u': #weak until
#         return f'{formula}(( ( {OLD_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))}) until ({OLD_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} ) ) or (always({OLD_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))}) ) )'
    
#     elif nodes[i].data == 'since':
#         return f'{formula}( ( {OLD_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))}) since ({OLD_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} ) )'
    
    
#     else: #PREDICATE
#         return nodes[i].data
    
    
    
# def OLD_NATURAL_tree_to_formula(nodes):
    
#     ## HYPOTHESIS:
#     ## nodes[0] = root
    
#     start = '' 
#     formula = f'{OLD_NATURAL_translate_node(start, nodes,  0 )}'
    
#     return formula
    
# def OLD_NATURAL_translate_node(formula, nodes, i):
    
#     '''As function translate_node but with implies and biimplication explicit (not expressed as combination of or/not)'''
    

#     if nodes[i].data == 'or':
#         return f'{formula}( {OLD_NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) or ( {OLD_NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} )'
    
#     elif nodes[i].data == 'and':
#         return f'{formula}( {OLD_NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) and ( {OLD_NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} )'
      
#     elif nodes[i].data == 'not':
#         return  f'{formula}not( {OLD_NATURAL_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} )'
#     #weak next
#     elif nodes[i].data == 'next':
#         return  f'{formula}next( {OLD_NATURAL_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} )'
#     #strong next
#     elif nodes[i].data == 's_next':
#         return  f'{formula}nextS( {OLD_NATURAL_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} )'
    
#     elif nodes[i].data == 'historically':
#         return  f'{formula}historically( {OLD_NATURAL_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} )'
    
#     elif nodes[i].data == 'once':
#         return  f'{formula}once( {OLD_NATURAL_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} )'
    
#     elif nodes[i].data == 'prev':
#         return  f'{formula}prev( {OLD_NATURAL_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} )'
#     elif nodes[i].data == 's_prev':
#         return  f'{formula}s_prev( {OLD_NATURAL_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} )'
    
#     ##Binary
    
#     elif nodes[i].data == '->':
#           return f'{formula}( {OLD_NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) --> ( {OLD_NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} )'
    
#     elif nodes[i].data == '<->':
#           return f'{formula}( {OLD_NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) <--> ( {OLD_NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} )'
  
#     elif nodes[i].data == 'eventually':
#         return  f'{formula}Eventually( {OLD_NATURAL_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} )'
    
#     elif nodes[i].data == 'always':
#         return  f'{formula}Always( {OLD_NATURAL_translate_node(formula, nodes, nodes.index(nodes[i].leftchild))} )'
    
#     elif nodes[i].data == 'until':
#           return f'{formula}( {OLD_NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) Until ( {OLD_NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} )'
  
#     elif nodes[i].data == 'weak_u':
#           return f'{formula}( {OLD_NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) WeakUntil ( {OLD_NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} )'
  
#     elif nodes[i].data == 'since':
#           return f'{formula}( {OLD_NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].leftchild))} ) Since ( {OLD_NATURAL_translate_node(formula,nodes, nodes.index(nodes[i].rightchild))} )'
  
  
#     else: #PREDICATE
#         return nodes[i].data
    



# # LENGTH = 3
# def create_tree_for_RNI():
#     '''Relational Non Interference'''
    
#     nodes = []
    
#     nodes.append(BinaryTreeNode('implies'))
    
#     nodes.append(BinaryTreeNode('predicate0'))
    
#     nodes[0].leftchild = nodes[1]
#     nodes[1].parent = nodes[0]
    
#     nodes.append(BinaryTreeNode('predicate0'))
    
#     nodes[0].rightchild = nodes[2]
#     nodes[2].parent = nodes[0]
    
#     return nodes


# # LENGTH =3
# def create_tree_for_GNI():
#     '''Generalized Non Interference'''
    
#     nodes = []
    
#     nodes.append(BinaryTreeNode('and'))
    
#     nodes.append(BinaryTreeNode('predicate0'))
    
#     nodes[0].leftchild = nodes[1]
#     nodes[1].parent = nodes[0]
    
#     nodes.append(BinaryTreeNode('predicate0'))
    
#     nodes[0].rightchild = nodes[2]
#     nodes[2].parent = nodes[0]
    
#     return nodes