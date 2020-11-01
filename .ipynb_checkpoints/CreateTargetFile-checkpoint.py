# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 18:10:05 2020

@author: Boris
"""

import glob

def readConditions(sdrf_file):
    handle = open(sdrf_file, 'r')
    next(handle)
    out = {}
    for line in handle:
        out[line.split('\t')[0].split()[0]] = line.split('\t')[1].split()[1:]
    handle2 = open("output.txt", "w")
    handle2.write("FileName" + '\t' + "DiseaseStatus" + '\n')
    for element in out:
        line = str(str(element) + '\t' + str(out[element]) + '\n')
        handle2.write(line)
    handle2.close()
    return out

readConditions("E-GEOD-45666.sdrf.txt")

final_out = {}
for file in glob.glob("GSM*"):
    #print(file)
    fullname = str(file)
    name = fullname.split('_')[0]
    template = readConditions("E-GEOD-45666.sdrf.txt")
    final_out[fullname] = template[name]
    
handle3 = open("targets.txt", "w")
handle3.write("FileName" + '\t' + "DiseaseStatus" + '\n')
for element in final_out:
    lijst = final_out[element]
    string = " "
    string = string.join(lijst)
    handle3.write(element + '\t' + string + '\n')
handle3.close()
    
    
    
        
        
