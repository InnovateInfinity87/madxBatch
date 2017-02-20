# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 11:59:12 2015
@author: mfraser adapted by wvandepo
MAD-X universal Twiss output reader/plotter
needs python 2.7 or up, dictionary comprehension
"""

# import packages required
import numpy as np
import os  
import Constants as c


# function to read twiss file into a dictionary
def TwissRead(filename):
    # open file
    f = open(filename,'r')

    # assign parameters and data variables to dictionaries: look for title header using '@' and '*'
    parameters = {}
    for line in f:
        line = line.strip()
        columns = line.split()
        if columns[0] is '@':
            parameters[columns[1]] = columns[3]
        if columns[0] is '*':
            variables = columns
            variables.remove('*')
            variable_index = {variables[x]:x for x in range(len(variables))}
            variable_data = {variables[x]:[] for x in range(len(variables))}
            f.next() # skip line under header titles
            break

    for line in f:
        line = line.strip()
        columns = line.split()
        for index in variable_index:
            if is_number(columns[variable_index[index]]):
                variable_data[index].append(float(columns[variable_index[index]]))
            else:
                variable_data[index].append(columns[variable_index[index]])
    
    # convert data into numpy arrays
    for index in variable_index:
        if all([isinstance(x, float) for x in variable_data[index]]):
            variable_data[index] = np.asarray(variable_data[index])
                
    f.close()
    
    return variable_data, parameters
    
def TrackRead():
    Tracks = {}
    current=0
    for fn in os.listdir(c.tracksdir):
        #print fn
        batch =  int(fn.split(".")[1][5:])
        if(current!=batch):
            print 'Processing batch '+str(batch)
            current=batch
        place =  int(fn.split(".")[2][3:])
        particle =  int(fn.split(".")[3][1:])
        Tracks[(batch,place,particle)]=SingleTrack(fn)
    return Tracks

def SingleTrack(filename):
    # open file
    f = open(filename,'r')
    
    for line in f:
        line = line.strip()
        columns = line.split()
        if columns[0] is '*':
            variables = columns
            variables.remove('*')
            variable_index = {variables[x]:x for x in range(len(variables))}
            variable_data = {variables[x]:[] for x in range(len(variables))}
            f.next() # skip line under header titles
            break 
        
    for line in f:
        line = line.strip()
        columns = line.split()
        for index in variable_index:
            if is_number(columns[variable_index[index]]):
                variable_data[index].append(float(columns[variable_index[index]]))
            else:
                variable_data[index].append(columns[variable_index[index]])
    
    # convert data into numpy arrays
    for index in variable_index:
        if all([isinstance(x, float) for x in variable_data[index]]):
            variable_data[index] = np.asarray(variable_data[index])
                
    f.close()
    
    return variable_data
    
def SingleLossFile(filename):
    # open file
    f = open(filename,'r')
    index=False
    
    for line in f:
        line = line.strip()
        columns = line.split()
        if columns[0] is '*':
            index=True
            variables = columns
            variables.remove('*')
            variable_index = {variables[x]:x for x in range(len(variables))}
            variable_data = {variables[x]:[] for x in range(len(variables))}
            f.next() # skip line under header titles
            break
    if(not index):
        return -1

    for line in f:
        #print "testdata"
        line = line.strip()
        columns = line.split()
        for index in variable_index:
            if is_number(columns[variable_index[index]]):
                variable_data[index].append(float(columns[variable_index[index]]))
            else:
                variable_data[index].append(columns[variable_index[index]])
    
    # convert data into numpy arrays
    for index in variable_index:
        if all([isinstance(x, float) for x in variable_data[index]]):
            variable_data[index] = np.asarray(variable_data[index])
                
    f.close()
    
    return variable_data
    
def is_number(num):
    try:
        float(num)
        return True
    except ValueError:
        return False
        