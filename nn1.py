#!/bin/python3
import os, sys 
# from ... import LHEImport
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
import vector
import mplhep as hep
import xml.etree.ElementTree as ET
from LHEImport.LHEImport2 import read_lhe, tohdf5

current_dir = os.path.dirname(os.path.realpath(__file__))
lhe_dir = os.path.join(current_dir, 'mil_lhe_files/')

#
# # for file in os.listdir(lhe_dir):
# #     data = read_lhe(lhe_dir +'/' + file)
# #     tohdf5(data, current_dir + '/hdf5/240121_2', key=file.split('.')[0])
#
# data = read_lhe(lhe_dir+'run_1.lhe')
data = read_lhe('./event_files/lhe/eft/run_09.lhe')
df = pd.DataFrame(data)
print(df[0][0].__dict__)



# eventdict={}
# for event, element in ET.iterparse("./mil_lhe_files/run_1.lhe", events=["end"]):
#     if element.tag == "initrwgt":
#         for initrwgtel in element:
#             if initrwgtel.tag == "weightgroup":
#                 eventdict["weightinfo"]={}
#                 for weightgroupel in initrwgtel:
#                     if weightgroupel.tag=="weight":
#                         id = str(weightgroupel.attrib["id"])
#                         eventdict["weightinfo"][id] = str(weightgroupel.text).split(' #')[0].split(' ')[-2:]
# #
# print(eventdict)

        
