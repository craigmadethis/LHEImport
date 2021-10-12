import math
import pandas as pd
import importlib.resources as pkg_resources

def convert_px_py_pz(px,py,pz):
    # transverse momentum
    pt = math.sqrt(math.fsum([math.pow(px, 2), math.pow(py, 2)]))
    # total momentum
    p = math.sqrt(math.fsum([math.pow(pt, 2), math.pow(pz, 2)]))
    if pt != 0:
        eta = math.asinh(pz / pt)
        phi = math.asin(py / pt)
    else:
        eta = math.copysign(float('inf'), pz)
        phi = 0
    return pt, eta, phi

def pdgid_to_string(pdgid):
    stream = pkg_resources.open_text(__package__, 'pdgid_string.csv')
    pdgid_data = pd.read_csv(stream)
    pdgid_data=pdgid_data.set_index('ID')
    return pdgid_data.loc[pdgid]['Name'], pdgid_data.loc[pdgid]['Latex']
