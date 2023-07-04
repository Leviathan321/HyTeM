#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
import pandas as pd
import math
import numpy as np

import sys
import rtamt
# from rtamt.spec.stl.discrete_time.specification import Semantics

import os
import log_dining_philosophers


def data_sum_start():
    
    dataset = [  [2,3,5],  [3,2, 5], 
                 [1,9,10], [9,1,10], 
                 [2,4,6],  [4,2,6], 
                 [2,6,8],  [6,2,8],
                 [4,4,8],  
                 [1,0,1],  [0,1,1],
                 [0,3,3],  
                 [5,3,8], 
                 [8,8,16], [9,7,16]]
    
    return dataset

def data_sum(n, seed):
    
    '''The function samples n/2 pairs of integer values in the range (0, 1000)'''
    n = int(n/2)
    # dataset = data_sum_start()
    dataset = []
    par_range = [0,1000]
    random.seed(seed)
    for i in range(n):
        addend1 = random.randrange(par_range[0], par_range[1]+1)
        addend2 = random.randrange(par_range[0], par_range[1]+1)
        result = addend1 + addend2
        
        dataset.append([addend1, addend2, result])
        dataset.append([addend2, addend1, result])

    return dataset

def data_random_format_sum(n, seed):
    
    dataset = []
    par_range = [0,1000]
    random.seed(seed)
    for i in range(n):
        addend1 = random.randrange(par_range[0], par_range[1]+1)
        addend2 = random.randrange(par_range[0], par_range[1]+1)
        result =  random.randrange(par_range[0], par_range[1]+1)
        
        dataset.append([addend1, addend2, result])

    return dataset


def data_trig(n, seed):
   
    dataset = []
    random.seed(seed)
    par_range = [-1,1]
     
    for i in range(0,n):
        value = random.uniform(par_range[0], par_range[1])
        
        dataset.append([value, math.acos(value), math.asin(value), math.cos(value), math.sin(value)])

    return dataset

def data_random_format_trig(n, seed):
   
    dataset = []
    random.seed(seed)
    par_range = [-1,1]
     
    for i in range(0,n):
        value = random.uniform(par_range[0], par_range[1])
        dataset.append([value, random.uniform(par_range[0], par_range[1]), random.uniform(par_range[0], par_range[1]), random.uniform(par_range[0], par_range[1]),random.uniform(par_range[0], par_range[1])])

    return dataset

            
def serial_adder(input1, input2):
    
    '''The two inputs are assumed to have the same length '''
    
    ##Reverse the two inputs?
    
    output = []
    
    carry_on = False # 'riporto' set to False
    
    for i in range(len(input1)):
        
        if   input1[i] == 0 and input2[i] == 0 and carry_on == False: output.append(0)
            
            
        elif input1[i] == 0 and input2[i] == 1 and carry_on == False: output.append(1)
            
        elif input1[i] == 1 and input2[i] == 0 and carry_on == False: output.append(1)
        
        elif input1[i] == 1 and input2[i] == 1 and carry_on == False:
            output.append(0)
            carry_on = True
            
        elif input1[i] == 0 and input2[i] == 0 and carry_on == True:
            output.append(1)
            carry_on = False
            
        elif input1[i] == 0 and input2[i] == 1 and carry_on == True:  output.append(0)
            
        elif input1[i] == 1 and input2[i] == 0 and carry_on == True:  output.append(0)
            
        elif input1[i] == 1 and input2[i] == 1 and carry_on == True:  output.append(1)
            
        else: print('Error in the definition of the sum')
        
        
    #When carry-on is true in the end --> add another 1
    # if carry_on == True: output.append(1) ##With this: ERROR: The input pi_2_2 does not have the same number of samples as time
    
    
     ##Reverse the output?
    return output



def data_serial_adder(length_time, numb_traces, seed):

    dataset = []
    
    random.seed(seed)
    
    #For one trace
    for i_trace in range(1):
        input1 = []
        input2 = []
        output = []
           
        #For each time frame
        for i_time in range(length_time):
            
            #Generation of the two inputs
            input1.append(random.randint(0, 1))
            input2.append(random.randint(0, 1))
            
        output = serial_adder(input1, input2)
        
        dataset.append([input1 , input2 , output])   ## EXTRA 1
    
    #One with same inputs --> same outputs
    dataset.append([input1 , input2 , output])       ## EXTRA 1
    
    #One with inputs swapped --> same outputs
    dataset.append([input2 , input1 , output])      ## EXTRA 1
    
    #Same first input
    input2 = []   
    for i_time in range(length_time):
        input2.append(random.randint(0, 1))
    
    output = serial_adder(input1, input2)
    dataset.append([input1 , input2 , output])       ## EXTRA 1
    
    #Same second input
    input1 = []
    for i_time in range(length_time):
            input1.append(random.randint(0, 1))
            
    output = serial_adder(input1, input2)
    dataset.append([input1 , input2 , output])       ## EXTRA 1
    
    #generate the other 'normal' traces
    for i_trace in range(0, numb_traces - 6):
        input1 = []
        input2 = []
        output = []
           
        #For each time frame
        for i_time in range(length_time):
            
            #Generation of the two inputs
            input1.append(random.randint(0, 1))
            input2.append(random.randint(0, 1))
            
        output = serial_adder(input1, input2)
        dataset.append([input1 , input2 , output]) 
        
    #One same inputs --> same outputs    
    dataset.append([input1 , input2 , output])      ## EXTRA 1
       
    return dataset

def data_serial_adder_manual():
    
    dataset = [[  [0,0,0,0,1,0,0,1,0,1], [0,1,0,1,1,0,1,1,0,1] , [0,1,0,1,0,1,1,0,1,0]],
               [  [0,0,0,0,1,0,0,1,0,1], [0,1,0,1,1,0,1,1,0,1] , [0,1,0,1,0,1,1,0,1,0]],
               [  [0,0,0,0,1,1,1,1,0,1], [0,1,0,1,1,0,1,1,0,1] , [0,1,0,1,0,0,1,1,1,0]],
               [  [0,0,0,0,1,0,0,1,0,1], [1,1,0,1,1,0,1,1,0,1] , [1,1,0,1,0,1,1,0,1,0]],
               [  [0,0,0,0,1,0,0,1,0,1], [1,1,0,1,1,0,1,1,0,1] , [1,1,0,1,0,1,1,0,1,0]],
               [  [1,1,0,0,1,0,0,1,0,1], [0,0,0,1,1,0,1,1,0,1],  [1,1,0,1,0,1,1,0,1,0]],
               [  [0,0,0,0,1,0,0,1,0,1], [0,1,0,1,1,0,1,1,0,1] , [0,1,0,1,0,1,1,0,1,0]],
               [  [0,0,0,0,1,0,0,1,0,0], [0,1,0,1,1,0,1,1,0,0] , [0,1,0,1,0,1,1,0,1,0]],
               [  [0,0,0,0,1,0,0,1,0,1], [0,1,0,1,1,0,1,1,0,0] , [0,1,0,1,0,1,1,0,1,1]],
               ]
    
    return dataset


def data_next(length_time, numb_traces, seed):

    dataset = []
    
    random.seed(seed)
    
    #For each execution/trace
    for i_trace in range(numb_traces):
        data_input = []
           
        #For each time frame
        for i_time in range(length_time):
            #Generation of the input
            data_input.append(random.randint(0, 1))
            
        ##PADDING WITH 1 ZERO TO HAVE THE TWO COMPONENTS WITH THE SAME LENGTH
        ## THE NEXT FORMULA IS STILL VALID
        dataset.append([data_input, data_input[1:]+[0]])
        
    return dataset


def data_next_next(length_time, numb_traces, seed):

    dataset = []
    random.seed(seed)
    
    #For each execution/trace
    for i_trace in range(numb_traces):
        data_input = []
           
        #For each time frame
        for i_time in range(length_time):
            #Generation of the input
            data_input.append(random.randint(0, 1))
        
        ##PADDING WITH 2 ZEROS TO HAVE THE TWO COMPONENTS WITH THE SAME LENGTH
        ## THE NEXT NEXT FORMULA IS STILL VALID
        dataset.append([data_input, data_input[2:]+[0,0]])
        
         
    # random.seed(int(time.time()))
    return dataset




def data_always(length_time, numb_traces, seed):
    
    dataset = []
    random.seed(seed)
    
    spec = rtamt.STLSpecification(language=rtamt.Language.PYTHON,semantics=rtamt.Semantics.STANDARD)
    spec.name = 'STL Discrete-time Offline monitor'
    
    spec.spec = 'always(x)'
    spec.declare_var('x', 'float')
    
    #For each execution/trace
    for i_trace in range(numb_traces):
        data_input = []
           
        #For each time frame
        for i_time in range(length_time):
            #Generation of the input
            data_input.append(random.randint(0, 1))
        dataSet = { 'time': list(np.arange(len(data_input)))}
        dataSet.update({'x' : data_input})
        
        try:
            spec.parse()
        
        except rtamt.STLParseException as err:
            
            print('STL Parse Exception: {}'.format(err))
            sys.exit()
        
        rho = spec.evaluate(dataSet)
        
        dataset.append([data_input, [rho[i][1] for i in range(len(rho))] ])
    
    return dataset   

    
def data_eventually(length_time, numb_traces, seed):
    
    dataset = []
    random.seed(seed)
    
    spec = rtamt.STLSpecification(language=rtamt.Language.PYTHON,semantics=rtamt.Semantics.STANDARD)
    spec.name = 'STL Discrete-time Offline monitor'
    
    spec.spec = 'eventually(x)'
    spec.declare_var('x', 'float')
    
    #For each execution/trace
    for i_trace in range(numb_traces):
        data_input = []
           
        #For each time frame
        for i_time in range(length_time):
            #Generation of the input
            data_input.append(random.randint(0, 1))
            
        dataSet = { 'time': list(np.arange(len(data_input)))}
        dataSet.update({'x' : data_input})
        
        try:
            spec.parse()
        
        except rtamt.STLParseException as err:
            
            print('STL Parse Exception: {}'.format(err))
            sys.exit()
        
        rho = spec.evaluate(dataSet)
        
        dataset.append([data_input, [rho[i][1] for i in range(len(rho))] ])
    
    return dataset   

    


def data_until(length_time, numb_traces, seed):

    dataset = []
    
    random.seed(seed)
    
    spec = rtamt.STLSpecification(language=rtamt.Language.PYTHON,semantics=rtamt.Semantics.STANDARD)
    spec.name = 'STL Discrete-time Offline monitor'
    
    spec.spec = '( ( x ) until (y) )'
    spec.declare_var('x', 'float')
    spec.declare_var('y', 'float')
    
    #For one trace
    for i_trace in range(numb_traces):
        input1 = []
        input2 = []
           
        #For each time frame
        for i_time in range(length_time):
            
            #Generation of the two inputs
            input1.append(random.randint(0, 1))
            input2.append(random.randint(0, 1))
           
        dataSet = { 'time': list(np.arange(len(input1)))}
        dataSet.update({'x' : input1})
        dataSet.update({'y' : input2})
        
        try:
            spec.parse()
        
        except rtamt.STLParseException as err:
            
            print('STL Parse Exception: {}'.format(err))
            sys.exit()
        
        rho = spec.evaluate(dataSet)
        
        dataset.append([input1, input2, [rho[i][1] for i in range(len(rho))] ])
    
    return dataset



def data_random_binary(length_time, numb_traces, seed):
    dataset = []
    random.seed(seed)

    # For each execution/trace
    for i_trace in range(numb_traces):
        data_input = []
        data_output = []

        # For each time frame
        for i_time in range(length_time):
            # Generation of the input
            data_input.append(random.randint(0, 1))
            data_output.append(random.randint(0, 1))

        dataset.append([data_input, data_output])

    return dataset
    

def data_random_three(length_time, numb_traces, seed):
    dataset = []
    random.seed(seed)

    # For each execution/trace
    for i_trace in range(numb_traces):
        data_input1 = []
        data_input2 = []
        data_output = []

        # For each time frame
        for i_time in range(length_time):
            # Generation of the input
            data_input1.append(random.randint(0, 1))
            data_input2.append(random.randint(0, 1))
            data_output.append(random.randint(0, 1))

        dataset.append([data_input1, data_input2, data_output])

    return dataset




def CollectCarlaData():
    
    # path = 'data/Carla/'
    # dataset_name = 'data'
    
    path = 'data/DataCarla/'
    dataset_name = 'disttrigger_'
    
    traces = []
    list_length = []
    
    n = 800
    for p in range(0, n): 
    
        data = pd.read_csv(f"{path}{dataset_name}{p}.csv", header = None)
        
        ## scenario = Day + Clear + Adult WITHOUT CRASH
        if data[1][0].strip('][').split(', ')[0] == 'day'\
           and data[1][1].strip('][').split(', ')[0] == 'clear'\
           and data[1][2].strip('][').split(', ')[0] == 'adult':
               
            crash = (data[1][5].strip('][').split(', '))[0]
            
            if crash == 'False': crash = False
            elif crash == 'True': crash = True
            
            if crash == False:
                # time = [ float(item) for item in data[1][6].strip('][').split(', ')]
                ego_y = [ float(item) for item in data[1][8].strip('][').split(', ')]
                ego_v_y = [ float(item) for item in data[1][10].strip('][').split(', ') ]
                # ped_y = [ float(item) for item in data[1][14].strip('][').split(', ') ]
                
                # traces.append([time ,crash,  ego_y, ego_v_y, ped_y  ])
                traces.append([ ego_v_y , ego_y ])
                if len(ego_v_y) <= len(ego_y): list_length.append(len(ego_v_y))
                else: list_length.append(len(ego_y))
     
    #Adjust traces such that they all have the same length
    
    min_length = min(list_length)
    
    print(min_length)
    for item in traces:
        #cut the traces to the min_length
        item[0] = item[0][0:min_length]
        item[1] = item[1][0:min_length]
        
    
    return traces

def collect_parking_data(numb_traces, seed):

    np.random.seed(seed)
    
    path = '../../../data/DataCarla/'
    # path = 'data/DataCarla/'
    dataset_name = 'disttrigger_'
    
    traces = []
    list_length = []
    
    n = 800
    for p in range(0, n): 
    
        data = pd.read_csv(f"{path}{dataset_name}{p}.csv", header = None)
        
        ## scenario = Day + Clear + Adult WITHOUT CRASH
        if data[1][0].strip('][').split(', ')[0] == 'day'\
           and data[1][1].strip('][').split(', ')[0] == 'clear'\
           and data[1][2].strip('][').split(', ')[0] == 'adult':
               
            crash = (data[1][5].strip('][').split(', '))[0]
            
            if crash == 'False': crash = False
            elif crash == 'True': crash = True
            
            if crash == False:
                # time = [ float(item) for item in data[1][6].strip('][').split(', ')]
                ego_y = [ float(item) for item in data[1][8].strip('][').split(', ')]
                ego_v_y = [ float(item) for item in data[1][10].strip('][').split(', ') ]
                
                traces.append([ ego_v_y , ego_y ])
                if len(ego_v_y) <= len(ego_y): list_length.append(len(ego_v_y))
                else: list_length.append(len(ego_y))
     
    #Adjust traces such that they all have the same length
    min_length = min(list_length)
    
    for item in traces:
        #cut the traces to the min_length
        item[0] = item[0][0:min_length]
        item[1] = item[1][0:min_length]
        
    #Shuffle and take two dijoint sets of numb_traces each
    np.random.shuffle(traces)
    
    dataset1 = traces[:numb_traces]
    dataset2 = traces[numb_traces:]
    
    return dataset1, dataset2, min_length # 171 : min_length


def data_random_binary_real(length_time, numb_traces, seed):
    dataset = []
    random.seed(seed)

    # For each execution/trace
    for i_trace in range(numb_traces):
        data_input = []
        data_output = []

        # For each time frame
        for i_time in range(length_time):
            # Generation of the input
            data_input.append(random.uniform(0,50))
            data_output.append(random.uniform(0 ,50 ))

        dataset.append([data_input, data_output])

    return dataset


def data_SoC():
    
    traces = []
    path = 'data/mirco-controller/'
    
    list_files = os.listdir(f'{path}')
    list_files.sort()
    maximum = 0 
    
    for trace_id in range(100):
        # print(trace_id)
        aux = []
       
        for i_file, file in enumerate(list_files):
            print(i_file)
            current_dataset = pd.read_csv(f"{path}/{file}", header = None, sep = '\t')
            
            aux.append(current_dataset[2][np.where(current_dataset[0]==trace_id)[0][:15884]].values) #15884 is the minimum length
            # maximum_el= np.max(current_dataset[2][np.where(current_dataset[0]==trace_id)[0][:15884]].values)
            # maximum= max(maximum_el, maximum)

        traces.append(aux)
    
    
    return traces 

def data_microcontroller(numb_traces, seed):
    
    np.random.seed(seed)
    
    traces = []
    path = '../../../data/mirco-controller/'
    
    list_files = os.listdir(f'{path}')
    list_files.sort()
    
    for trace_id in range(100):
        aux = []
       
        for i_file, file in enumerate(list_files):
            # print(i_file)
            current_dataset = pd.read_csv(f"{path}/{file}", header = None, sep = '\t')
            
            aux.append(current_dataset[2][np.where(current_dataset[0]==trace_id)[0][:15884]].values)#15884 is the minimum length
        
        traces.append(aux)
    
    #Shuffle and take two dijoint sets of numb_traces each
    np.random.shuffle(traces)
    
    dataset1 = traces[:numb_traces]
    dataset2 = traces[numb_traces:]
    
    return dataset1, dataset2, 15884#15884 is the minimum length


def data_random_microcontroller(length_time, numb_traces, seed):
    dataset = []
    random.seed(seed)
    
    numb_variables = 35
    
    
    # For each execution/trace
    for i_trace in range(numb_traces):
        data = [[]] * numb_variables
        for var in range(numb_variables):
            # For each time frame
            aux = []
            for i_time in range(length_time):
                # Generation of the input
                aux.append(random.randint(0,300)) #the maximum value in Soc is 255
            data[var] = aux.copy()
        dataset.append(data)

    return dataset


def data_philosophers():
    
    
    traces = log_dining_philosophers.dining_philosophers()
    n = len(traces[0]) # number of variables  
    data1 = []
    data2 = []
    data3 = []
    for var in range(n):
        data1.append(traces[0][var][:363]) #1088 is the totoal length. 1088/3 = 363
        data2.append(traces[0][var][363:726]) #1088 is the totoal length. 1088/3 = 363
        data3.append(traces[0][var][725:]) #1088 is the totoal length. 1088/3 = 363
    #To have the same format as log_dining philosophers
    data1 = [data1] 
    data2 = [data2] 
    data3 = [data3] 
   
    return data1, data2, data3
    
    
def data_philosophers_random(length_time, numb_traces, seed):   #363 
    
    dataset = []
    random.seed(seed)
    
    numb_variables = 5
    
    # For each execution/trace
    for i_trace in range(numb_traces):
        data = [[]] * numb_variables
        for var in range(numb_variables):
            # For each time frame
            aux = []
            for i_time in range(length_time):
                # Generation of the input
                aux.append(random.randint(0,2)) #the values of the actions are 0, 1,2
            data[var] = aux.copy()
        dataset.append(data)

    return dataset
    
    
    
    
    
    
    
    
