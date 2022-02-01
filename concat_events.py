from LHEImport.LHEImport2 import read_lhe, tohdf5
import os
import pandas as pd
DIR = "./mil_lhe_files/"
FILE = DIR+"run_1.lhe"
# for file in os.listdir(DIR):
#     data = read_lhe(DIR+file)
#     tohdf5(data, "./mil_lhe_files/mil", key=f"{file.split('.')[0]}")


data = read_lhe(FILE)
# tohdf5(data, "./mil_lhe_files/mil", key="hello")
# eventinfo= [e.eventinfo for e in events]
# particles = [e.particles for e in events]

events = [d for d in data]
eventinfo= [e.eventinfo for e in events]
particles = [e.particles for e in events]
weights = [e.weights for e in events]
weightinfo = [e.weightinfo for e in events]

print(weightinfo)
