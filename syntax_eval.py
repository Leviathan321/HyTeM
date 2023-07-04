#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
# sys.path.insert(0, 'rtamt/')

import itertools
import numpy as np

import rtamt
from rtamt.spec.stl.discrete_time.specification import Semantics

# import SMT
import math
import random
# import time
# import data_collection as dc

# import pickle
# import tree



def evaluate_rob_boolean(pi, formula, bool_interface_aware, input_variables, output_variables):
    
    '''
    For LTL in which we only want satisfaction or violation
    
    Computes the boolean satisfaction/violation of a (quantifier-free) formula 
    with respect to a specific tuple of traces. 
    It uses the RTAMT tool with boolean satisfaction: outputs are only +inf or -inf
    
    -bool_interface_aware: boolean, whether to use or not the prior knowledge of input-output variables
    -input_variables: None if bool_interface_aware is False, else the list of input variables
    -output_variables: None if bool_interface_aware is False, else the list of output variables
    
    '''
    
        
    # Transform pi[i][j] into pi_i_j because the tool does not support the brackets
    # new_formula = formula.replace('[','_')
    # new_formula = new_formula.replace(']','')
    
    new_formula = formula.replace('][','_')
    new_formula = new_formula.replace('pi[','pi_')
    # new_formula = new_formula.replace('] ',' ')
    
    # Indices of where the 'pi' start
    variable_indices = [i for i in range(len(new_formula)) if new_formula.startswith('pi', i)]
    
    
    spec = rtamt.STLSpecification(language=rtamt.Language.PYTHON,semantics=rtamt.Semantics.OUTPUT_ROBUSTNESS)
    # spec.semantics = Semantics.STANDARD
    
    spec.name = 'STL Discrete-time Offline monitor'
    
    dataSet = { 'time': list(np.arange(len(pi[0][0])))}
    
    already_declared = []
    
    for item in  variable_indices:
        indx1 = new_formula[item + 3]# +3 because it is the length of the string 'pi_'
        end = 4 # the index has at least length 1 ( = 4-3 )
        while True : # find the end of first index 
            if new_formula[item + end] == '_': break #the first index ends when '_' appears
            end += 1
        indx1 = int(new_formula[item + 3: item + end]  )
        
        start = item + end + 1 #from the position of '_' + 1 (start of second index)
        end = start + 1 
        while True: # find the end of second index 
            if new_formula[end] == ']' or new_formula[end] == ')': break #the second index ends when ' ' or ')' appear
            end += 1
        
        indx2 = int(new_formula[start : end])
        new_formula = new_formula[:end]+ ' ' + new_formula[end +1:]
        
        #Avoid already declared pairs of indices
        if f'{indx1}{indx2}' not in already_declared:
        
            spec.declare_var(f'pi_{indx1}_{indx2}', 'float')
            
            if bool_interface_aware == False: spec.set_var_io_type(f'pi_{indx1}_{indx2}', 'input')
            
            else:
               
                if indx2 in input_variables :    spec.set_var_io_type(f'pi_{indx1}_{indx2}', 'input')
                elif indx2 in output_variables:  spec.set_var_io_type(f'pi_{indx1}_{indx2}', 'output')
                else:
                    print(f'\nVariable {indx2} not specified either as input or output variable!\n')
                    sys.exit()
                    
            dataSet.update({f'pi_{indx1}_{indx2}' : pi[int(indx1)][int(indx2)]})
            
            already_declared.append(f'{indx1}{indx2}')
    
    if 'True' in new_formula:
        spec.declare_const('True', 'int', 1)
        spec.set_var_io_type('True', 'input')
            
        
    spec.spec = new_formula
    
    try:
        spec.parse()
        
    except rtamt.STLParseException as err:
        print(new_formula)
        print('STL Parse Exception: {}'.format(err))
        sys.exit()
        
    rho = spec.evaluate(dataSet)
    
   
    # JUST CHECK SATISFACTION VS VIOLATION
    if bool_interface_aware == False:
        
        if rho[0][1] == 0: print('\n\n\n\n\n\n\n\n\n\n\n\n\nROB = 0 !!!!\n\n\\n\n\n\n\n\n\n\n\n\n\n')
    
        if rho[0][1] >= 0: return 'sat'
        
        else: return 'unsat'
        
    #CHECK NON - VACOUSLY SATISFACTION    
    elif bool_interface_aware == True:
        if rho[0][1] == np.inf: 
            # print('skipped')
            return 'skipped'
    
        elif rho[0][1] >= 0:
            # if rho[0][1] == 0: print('che fare?')
            # print('sat')
            return 'sat'
        
        elif rho[0][1]< 0: return 'unsat'
            

def evaluate_rob_quantitative(pi, formula):
    
    '''Computes the QUANTITATIVE satisfaction/violation of a (quantifier-free) formula 
    with respect to a specific tuple of traces. 
    It uses the RTAMT tool with quantitative robustness 
    
    
     ---> DECIDE HOW TO HANDLE ROB = 0 <---
     
     
     '''
    # Transform pi[1][2] into pi_1_2 because the tool does not support the brackets
    # new_formula = formula.replace('[','_')
    # new_formula = new_formula.replace(']','')
    
    new_formula = formula.replace('][','_')
    new_formula = new_formula.replace('pi[','pi_')
    # new_formula = new_formula.replace('] ',' ')
    
    # Indices of where the 'pi' start
    variable_indices = [i for i in range(len(new_formula)) if new_formula.startswith('pi', i)]
    
    spec = rtamt.STLSpecification(language=rtamt.Language.PYTHON)
    spec.name = 'STL Discrete-time Offline monitor'
    
    dataSet = { 'time': list(np.arange(len(pi[0][0])))}
    
    already_declared = []
    
    for item in  variable_indices:
        indx1 = new_formula[item + 3] # +3 because it is the length of the string 'pi_'
        end = 4 # the index has at least length 1 ( = 4-3 )
        while True : # find the end of first index 
            if new_formula[item + end] == '_': break #the first index ends when '_' appears
            end += 1
        indx1 = int(new_formula[item + 3: item + end]  )
        
        start = item + end + 1 #from the position of '_' + 1 (start of second index)
        end = start + 1 
        while True: # find the end of second index 
            if new_formula[end] == ']' or new_formula[end] == ')': break #the second index ends when ' ' or ')' appear
            end += 1
        
        indx2 = int(new_formula[start : end])
        new_formula = new_formula[:end]+ ' ' + new_formula[end +1:]
        
        #Avoid already declared pairs of indices
        if f'{indx1}{indx2}' not in already_declared:
        
            spec.declare_var(f'pi_{indx1}_{indx2}', 'float')
        
            dataSet.update({f'pi_{indx1}_{indx2}' : pi[int(indx1)][int(indx2)]})
            
            already_declared.append(f'{indx1}{indx2}')
            
    spec.spec = new_formula
    spec.semantics = Semantics.STANDARD
    
    try:
        spec.parse()
        
    except rtamt.STLParseException as err:
        print(new_formula)
        print('STL Parse Exception: {}'.format(err))
        sys.exit()
        
    rho = spec.evaluate(dataSet)
    
    return rho[0][1]
    

def compute_cost_positive(traces, formula, quantifiers, bool_temporal_operators,  bool_interface_aware, input_variables, output_variables):
    
    '''The function comuputes the POSITIVE cost of a given formula.
    
    INPUTS:
        - traces: set of traces w.r.t which the score has to be evaluated.
        - formula: formula of the formula of which the score has to be computed
        - quantifiers: quantofoers to be applied to the formulas's trace variables
        - bool_temporal operators: 
            If bool_temporal_operators == True: there are temporal operators in the grammar --> use RTAMT
            If bool_temporal_operators == False: there are NOT temporal operators in the grammar --> use python function evaluate
       
    '''
    
    m = len(quantifiers)
    numb_traces = len(traces)
    
    results = []
    
    ## REMIND TO CHANGE ALSO INSIDE THE LAST LOOP OF THE FUNCTION!!!
    ##Compare pi also with itself
    # indices_traces_combo = [item for item in itertools.product(range(numb_traces), repeat=m)]
    
    ##Do not compare pi  with itself
    indices_traces_combo = [item for item in itertools.permutations(range(numb_traces), m)]
    
    ##STUDY SATISFACTION/VIOLATION OF EACH COMBINATION OF TRACES
    for ind in indices_traces_combo:
        
        pi = [traces[item] for item in ind] #pi is used inside eval(formula)
        
        #Positive --> Counts the number of satisfaction 
        
        # if SMT.evaluate_satisfaction(z3formula, pi) == 'sat': results.append(1) #satisfaction
        
        if bool_temporal_operators == False: 
            if eval(formula): results.append(1) #satisfaction
            else: results.append(0) #violation
        #Use RTAMT
        elif bool_temporal_operators == True:
            if evaluate_rob_boolean(pi, formula, bool_interface_aware, input_variables, output_variables) == 'sat': results.append(1) #satisfaction
            else: results.append(0) #violation
    
    sat = results.copy()
      
    ## EVALUATE FORMULA WITH QUANTIFIERS
    for i in range(0,m): #loop over the quantifiers
        aux = []
        
        #POSITIVE --> counts the number of satisfaction 
        # Answer the question: 
        # Which is the minimum number of combinations (of m traces) outcomes that I need to change to get a violation?
           
        # If the evaluation with the formula itself is not admitted
        for j in range(0, int(math.factorial(numb_traces)/math.factorial(numb_traces - (m-i-1)))):
            #EXISTS
            if quantifiers[-i-1]==0:  aux.append(sum( sat[j*(numb_traces-(m - i - 1)): (j+1)*(numb_traces-(m - i - 1))]) ) #(#traces - #left quantifiers)
            #FORALL
            elif quantifiers[-i-1]== 1: aux.append(min( sat[j*(numb_traces-(m - i - 1)): (j+1)*(numb_traces-(m - i - 1))]) )
        
        # If the evaluation with the formula itself is admitted
        # for j in range(0, numb_traces**(m-i-1)):
            # if quantifiers[-i-1]==0:  aux.append(sum( sat[j*numb_traces: (j+1)*numb_traces]) )
            # elif quantifiers[-i-1]== 1: aux.append(min( sat[j*numb_traces: (j+1)*numb_traces]) )
          
        sat = aux.copy()
    
    return sat[0]
    

def compute_cost_negative(traces, formula, quantifiers, bool_temporal_operators,  bool_interface_aware, input_variables, output_variables):
    
    '''The function comuputes the NEGATIVE cost of a given formula.
    
    INPUTS:
        - traces: set of traces w.r.t which the score has to be evaluated.
        - formula:  formula of which the score has to be computed
        - quantifiers: quantifiers to be applied to the formulas's trace variables
        - bool_temporal operators: 
            If bool_temporal_operators == True: there are temporal operators in the grammar --> use RTAMT
            If bool_temporal_operators == False: there are NOT temporal operators in the grammar --> use python function evaluate
       
    '''
    
    m = len(quantifiers)
    numb_traces = len(traces)
    
    results = []
    
    ## REMIND TO CHANGE ALSO INSIDE THE LAST LOOP OF THE FUNCTION!!!
    #Compare pi also with itself
    # indices_traces_combo = [item for item in itertools.product(range(numb_traces), repeat=m)]
    
    ##Do not compare pi  with itself
    indices_traces_combo = [item for item in itertools.permutations(range(numb_traces), m)]
    
    ##STUDY SATISFACTION/VIOLATION OF EACH COMBINATION OF TRACES
    for ind in indices_traces_combo:
        
        pi = [traces[item] for item in ind] #pi is used inside eval(formula)
        
        if bool_temporal_operators == False:
            if eval(formula): results.append(0) #satisfaction
            else: results.append(-1) #violation
            
        elif bool_temporal_operators == True:
            if evaluate_rob_boolean(pi, formula, bool_interface_aware, input_variables, output_variables) == 'sat': results.append(0) #satisfaction
            else: results.append(-1) #violation
    
    sat = results.copy()
      
    ## EVALUATE FORMULA WITH QUANTIFIERS
    for i in range(0,m): #loop over the quantifiers
        
        aux = []
        
        #NEGATIVE --> counts the number of violations 
        # Answer the question: 
        # Which is the minimum number of combinations (of m traces) outcomes that I need to change to get satisfaction?
         
        # If the evaluation with the trace itself is not admitted
        for j in range(0, int(math.factorial(numb_traces)/math.factorial(numb_traces - (m-i-1)))):
            #EXISTS
            if quantifiers[-i-1]==0:  aux.append(max( sat[j*(numb_traces-(m - i - 1)): (j+1)*(numb_traces-(m - i - 1))]) ) # change: every time we study (numb_traces -1)
            #FORALL
            elif quantifiers[-i-1]== 1: aux.append(sum( sat[j*(numb_traces-(m - i - 1)): (j+1)*(numb_traces-(m - i - 1))]) )
        
        #If the evaluation with the trace itself is admitted   
        # for j in range(0, numb_traces**(m-i-1)):
            # if quantifiers[-i-1]==0:  aux.append(max( sat[j*numb_traces: (j+1)*numb_traces]) ) # change: every time we study (numb_traces -1)
            # elif quantifiers[-i-1]== 1: aux.append(sum( sat[j*numb_traces: (j+1)*numb_traces]) )
            
            
        sat = aux.copy()
    
    return sat[0]
    

def compute_cost(traces, formula, quantifiers, bool_temporal_operators,  bool_interface_aware, input_variables, output_variables):
    
    '''The function computes the positive cost if the most internal quantifier is EXISTS,
       and the negative cost if it is FOR ALL.
       In case the previous result is 0, the computation is repeated with the other measure of cost.
       
       If bool_temporal_operators == True: there are temporal operators in the grammar --> use RTAMT
       If bool_temporal_operators == False: there are NOT temporal operators in the grammar --> use python function evaluate
       '''
    # z3formula = SMT.tree_to_Z3(nodes)
    
    if quantifiers[-1] == 0: # if the most internal quantifier is EXISTS
    
        cost = compute_cost_positive(traces, formula, quantifiers, bool_temporal_operators,  bool_interface_aware, input_variables, output_variables)    
        if cost == 0: cost = compute_cost_negative(traces, formula, quantifiers, bool_temporal_operators,  bool_interface_aware, input_variables, output_variables)
    
    elif quantifiers[-1] == 1: # if the most internal quantifier is FOR ALL
    
        cost = compute_cost_negative(traces, formula, quantifiers, bool_temporal_operators,  bool_interface_aware, input_variables, output_variables)    
        if cost == 0: cost = compute_cost_positive(traces, formula, quantifiers, bool_temporal_operators,  bool_interface_aware, input_variables, output_variables)
    
    return cost


def efficient_monitoring(traces, formula, quantifiers, bool_temporal_operators, n_subset,  bool_interface_aware, input_variables, output_variables):
    
    '''The function computes whether the set of traces satisfies the formula with the given quantifiers.
    In this case, the cost is not computed: but only the boolean satisfaction of the hyperformula.
    In this case, the forall is stopped as soon as one tuple is violated, 
    and the exists is stopped as soon as one tuple is satisfied.
    
    - n_subset : total number of tuples to examine  or None. 
    If n_subset is False, we use it to return the counter of LTL/STL monitors by setting bool_return_counter to True
    
    OUTPUT: +1 (satisfied) or -1 (violated). With n_subset = False, OUTPUT: +-1 , counter (number of calls to ltl monitor)'''
    
    
    bool_return_counter = False
    
    #!!! Indices are tought for not evaluating the formula with the same trace repeated more than once
    m = len(quantifiers)
    numb_traces = len(traces)
    
    ##Do not compare pi  with itself
    indices_traces_combo = [item for item in itertools.permutations(range(numb_traces), m)]
    
    #Compare pi also with itself
    # indices_traces_combo = [item for item in itertools.product(range(numb_traces), repeat=m)]
    
    #Consider ONLY a subset of tuples to examine
    if n_subset is None: sampled_tuples = indices_traces_combo
    elif n_subset is False: 
        bool_return_counter = True
        sampled_tuples = indices_traces_combo
    else: sampled_tuples = random.sample(indices_traces_combo,numb_traces)
    
    # random.seed(int(time.time()))
    
    ind = [item for item in range(m)]
    
    index_to_skip = None #index to skip in case a satisfied/violated tuple has been found (E/F, respectively)
    value_to_skip  = None
    sat_skip = None #either 0 or 1
    
    counter = 0
    results = [[]] * m #Inizialization of the vector of the results. Each component refers to a quantifiers
    
    ## STUDY SATISFACTION/VIOLATION OF EACH COMBINATION OF TRACES
    #For each combinaion of traces
    for ind in indices_traces_combo:
        
        #If we can skip (due to already satisfied exists or already violated forall)
        if index_to_skip is not None: 
            #If we are still on the index to skip
            if ind[index_to_skip] == value_to_skip: 
                results[-1] = results[-1] + [sat_skip]
            #If we have skipped the index to skip    
            elif ind[index_to_skip] != value_to_skip:
                 index_to_skip = None
                 
        #If there are NOT already satisfied-exists or already violated-forall  
        if index_to_skip is None:
            #If the current tuple has been sampled to be evaluated
            if ind in sampled_tuples: 
            
                #Evaluate the satisfaction/violation of the current tuple of traces
                pi = [traces[item] for item in ind] #pi is used inside eval(formula)
                
                if bool_temporal_operators == False: 
                    if eval(formula): results[-1] = results[-1] + [1] #satisfaction
                    else: results[-1] = results[-1] + [0] #violation
                #Use RTAMT
                elif bool_temporal_operators == True:
                    counter += 1
                    
                    if counter > 200000: 
                        print('Runout')
                        if bool_return_counter : return 'runout', counter #-1#
                        else: return 'runout'
                    
                    rho = evaluate_rob_boolean(pi, formula,  bool_interface_aware, input_variables, output_variables)
                    if  rho == 'sat': results[-1] = results[-1] + [1] #satisfaction
                    elif rho == 'unsat': results[-1] = results[-1] + [0] #violation
                    else:  results[-1] = results[-1] + ['skipped']
            
            #If the current tuple has NOT been sampled to be evaluated -> skip the monitor  
            else:  results[-1] = results[-1] + ['skipped']
                
        #Loop over quantifiers starting from the most internal one
        for i in range(m-1 , 0 , -1): 
            # If the list is completed
            if len(results[i]) == numb_traces - i: ###! Change if same trace can be repeated (== numb_traces)
            
                if quantifiers[i] == 0: #EXISTS
                    if 1 in results[i]: results[i-1] = results[i-1] + [1] #satisfied
                    else: results[i-1] = results[i-1] + [0] #violated
                
                elif quantifiers[i] == 1: #FORALL 
                    if 0 in results[i]: results[i-1] = results[i-1] + [0] #violated
                    else: results[i-1] = results[i-1] + [1] #satisfied
                
                results[i] = []
                
            # If the list can be completed automatically with success (EXISTS) 
            elif (index_to_skip is None or index_to_skip > i-1) and quantifiers[i] == 0 and 1 in results[i]:  
                sat_skip = 1
                index_to_skip = i - 1
                value_to_skip = ind[i - 1]
            
            # If the list can be completed automatically with insuccess (FORALL) 
            elif (index_to_skip is None or index_to_skip > i-1) and quantifiers[i] == 1 and 0 in results[i]:  
                sat_skip = 0
                index_to_skip = i - 1
                value_to_skip = ind[i - 1]
            
            else: #Continue normally
                break
        #Check on the most external quantifier (early stop for EXISTS)
        if quantifiers[0] == 0 and 1 in results[0]: 
            if bool_return_counter: return 1, counter
            else: return 1 #satisfaction
        #Check on the most external quantifier (early stop for FORALL)
        elif quantifiers[0] == 1 and 0 in results[0]: 
            if bool_return_counter: return -1, counter
            else: return -1 #violation
            
    if quantifiers[0] == 0 :#EXISTS
        
        if 1 in results[0]: 
            if bool_return_counter: return 1, counter
            else: return 1 #satisfaction
        else: 
             if bool_return_counter: return -1, counter
             else: return -1 #violation
    
    elif quantifiers[0] == 1 : #FORALL
        if 0 in results[0]: 
             if bool_return_counter: return -1, counter
             else: return -1 #violation
        else: 
             if bool_return_counter: return 1, counter
             else: return 1 #satisfaction


   
def compute_robustness_STL(traces, formula, quantifiers):
    
    '''The function computes the robustness of the pair (quantifiers, formula) with respect to the current traces'''
    
    m = len(quantifiers)
    numb_traces = len(traces)
    
    results = []
    
    
    ## REMIND TO CHANGE ALSO INSIDE THE LAST LOOP OF THE FUNCTION!!!
    ## Compare pi also with itself
    # indices_traces_combo = [item for item in itertools.product(range(numb_traces), repeat=m)]
    
    ## Do not compare pi  with itself
    indices_traces_combo = [item for item in itertools.permutations(range(numb_traces), m)]
    
    ##STUDY SATISFACTION/VIOLATION OF EACH COMBINATION OF TRACES
    for ind in indices_traces_combo:
        pi = [traces[item] for item in ind]
        #Use RTAMT
        results.append(evaluate_rob_quantitative(pi, formula)) #append the computed robustness
        
    sat = results.copy()
      
    ## EVALUATE FORMULA WITH QUANTIFIERS
    for i in range(0,m): #loop over the quantifiers
        
        aux = []
        
        # If the evaluation with the formula itself is not admitted
        for j in range(0, int(math.factorial(numb_traces)/math.factorial(numb_traces - (m-i-1)))):
            #EXISTS
            if quantifiers[-i-1]==0:  aux.append(max( sat[j*(numb_traces-(m - i - 1)): (j+1)*(numb_traces-(m - i - 1))]) ) #(#traces - #left quantifiers)
            #FORALL
            elif quantifiers[-i-1]== 1: aux.append(min( sat[j*(numb_traces-(m - i - 1)): (j+1)*(numb_traces-(m - i - 1))]) )
        
        # If the evaluation with the formula itself is admitted
        # for j in range(0, numb_traces**(m-i-1)):
            #EXISTS
            # if quantifiers[-i-1]==0:  aux.append(max( sat[j*numb_traces: (j+1)*numb_traces]) )
            #FORALL
            # elif quantifiers[-i-1]== 1: aux.append(min( sat[j*numb_traces: (j+1)*numb_traces]) )
          
        sat = aux.copy()
    
    return sat[0]

def lookup_table_uniform_distance(traces):
    
    '''The function computes the lookup table that associate each pair of
        traces with their uniform distancee.
        Thee matrix is symmetric with diagonal of zeros'''
        
    numb_traces = len(traces)
    table = np.zeros((numb_traces, numb_traces))
    
    for i in range(numb_traces):
        matrix_trace_i = np.matrix(traces[i])
        
        for j in range(i+1, numb_traces):
           table[i,j] = np.linalg.norm(matrix_trace_i-np.matrix(traces[j]), ord = np.inf ) #uniform norm . multidimensional?
           table[j,i] = table[i,j]
    
    return table


def monitor_hyperSTL_correctness(traces, formula, quantifiers):
       
    '''The function computes the robustness of the pair (quantifiers, formula) with respect to the current traces
    using the correctness theorem for the quantitative semantics of STL:
        
    If w1 satsifies phi with robustness rho and |w1-w2|<rho then w2 satisfies phi as well'''
    
    
    ### !!!! SET TO TRUE TO RUN stl_monitoring.py for salability results
    bool_return_counter = False
    
    m = len(quantifiers)
    numb_traces = len(traces) 
    # print(formula)
    
    #Table of distances between pairs of traces
    lookuptable = lookup_table_uniform_distance(traces)
    
    #Table of the robustness
    rob_ind_table = [ ] # indices of the tuples of traces
    rob_rho_table = [ ] # indices of values of robustness
    
    ##Do not compare pi  with itself
    indices_traces_combo = [item for item in itertools.permutations(range(numb_traces), m)]
    
    #Compare pi also with itself
    # indices_traces_combo = [item for item in itertools.product(range(numb_traces), repeat=m)]
    
    ind = [item for item in range(m)]
    
    index_to_skip = None #index to skip in case a satisfied/violated tuple has been found (E/F, respectively)
    value_to_skip  = None
    sat_skip = None #either 0 or 1
    
    counter_monitor = 0
    counter_skip_correctness =0
    results = [[]] * m #Inizialization of the vector of the results. Each component refers to a quantifiers
    
    ## STUDY SATISFACTION/VIOLATION OF EACH COMBINATION OF TRACES
    #For each combinaion of traces
    for ind in indices_traces_combo:
        
        #If we can skip (due to already satisfied exists or already violated forall)
        if index_to_skip is not None: 
            #If we are still on the index to skip
            if ind[index_to_skip] == value_to_skip: 
                results[-1] = results[-1] + [sat_skip]
            #If we have skipped the index to skip    
            elif ind[index_to_skip] != value_to_skip:
                 index_to_skip = None
                 
        #If there are NOT already satisfied-exists or already violated-forall  
        if index_to_skip is None:
            
         ##If the current tuple has been sampled to be evaluated
         #if ind in indices_traces_combo: 
            
            skip_rob = False #if true: the evaluation of the robustness can be skipped,
            # because the correctness theorem can be applied
             
            #Compare the current indices with the indices of previously computed robustness
            for index_corr, item_corr in enumerate(rob_ind_table):
                
                #Indices of components with different values
                different_components = [i for i, (a, b) in enumerate(zip(item_corr, ind)) if a != b]
                
                #If there is a tuple of traces with all but one equal components
                if len(different_components) == 1:
                    
                   # Distance between the two traces corresponding to the index that differs in the tuple
                   distance = lookuptable[ item_corr[different_components[0]] ,ind[different_components[0]]]
                    
                   #If the distance is smaller than the abs value of the robustness
                   if distance < abs(rob_rho_table[index_corr]):
                       
                       if rob_rho_table[index_corr] >=0: results[-1] = results[-1] + [1] #satisfaction
                       elif rob_rho_table[index_corr]<0: results[-1] = results[-1] + [0] #violation
                       #skip the computation of the robustness for the current tuple of traces
                       skip_rob = True
                       counter_skip_correctness += 1
                       break
                
            #(Correctness theorem cannot be applied) The robustness needs to be computed
            if skip_rob == False:
                
                #Evaluate the satisfaction/violation of the current tuple of traces
                pi = [traces[item] for item in ind] #pi is used inside eval(formula)
                counter_monitor += 1
                if counter_monitor > 200000: 
                    if bool_return_counter : return 'runout', counter_monitor, counter_skip_correctness
                    else: return 'runout' #-1#
                
                rho = evaluate_rob_quantitative(pi, formula)
                
                #Update table of the robustness
                rob_ind_table.append(ind)
                rob_rho_table.append(rho)
                
                if  rho >= 0: results[-1] = results[-1] + [1] #satisfaction
                elif rho < 0: results[-1] = results[-1] + [0] #violation
            
            #If the current tuple has NOT been sampled to be evaluated -> skip the monitor  
            # else:  results[-1] = results[-1] + ['skipped']
                
        #Loop over quantifiers starting from the most internal one
        for i in range(m-1 , 0 , -1): 
            # If the list is completed
            if len(results[i]) == numb_traces - i: ###! Change if same trace can be repeated (== numb_traces)
            
                if quantifiers[i] == 0: #EXISTS
                    if 1 in results[i]: results[i-1] = results[i-1] + [1] #satisfied
                    else: results[i-1] = results[i-1] + [0] #violated
                
                elif quantifiers[i] == 1: #FORALL 
                    if 0 in results[i]: results[i-1] = results[i-1] + [0] #violated
                    else: results[i-1] = results[i-1] + [1] #satisfied
                
                results[i] = []
                
            # If the list can be completed automatically with success (EXISTS) 
            elif (index_to_skip is None or index_to_skip > i-1) and quantifiers[i] == 0 and 1 in results[i]:  
                sat_skip = 1
                index_to_skip = i - 1
                value_to_skip = ind[i - 1]
            
            # If the list can be completed automatically with insuccess (FORALL) 
            elif (index_to_skip is None or index_to_skip > i-1) and quantifiers[i] == 1 and 0 in results[i]:  
                sat_skip = 0
                index_to_skip = i - 1
                value_to_skip = ind[i - 1]
            
            else: #Continue normally
                break
        #Check on the most external quantifier (early stop for EXISTS)
        if quantifiers[0] == 0 and 1 in results[0]: 
            if bool_return_counter :return  1 , counter_monitor, counter_skip_correctness
            else: return 1 #satisfaction
        #Check on the most external quantifier (early stop for FORALL)
        elif quantifiers[0] == 1 and 0 in results[0]: 
            if bool_return_counter : return -1, counter_monitor, counter_skip_correctness
            else: return -1 #violation
            
    if quantifiers[0] == 0 :#EXISTS
        
        if 1 in results[0]: 
            if bool_return_counter : return 1, counter_monitor, counter_skip_correctness
            else:return 1 #satisfaction
        else: 
            if bool_return_counter : return -1, counter_monitor, counter_skip_correctness
            else: return -1 #violation
    
    elif quantifiers[0] == 1 :#FORALL
        if 0 in results[0]: 
             if bool_return_counter : return  -1, counter_monitor, counter_skip_correctness
             else: return -1 #violation
        else: 
             if bool_return_counter :  return 1, counter_monitor, counter_skip_correctness
             else: return 1 #satisfaction