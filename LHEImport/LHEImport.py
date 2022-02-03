import xml.etree.ElementTree as ET
# import math
import vector
import pandas as pd
import numpy as np
# import importlib.resources as pkg_resources

class LHEEvent:
    def __init__(self,  eventinfo, particles, weights=None, attributes=None, weightinfo=None):
        self.weightinfo = weightinfo
        self.eventinfo = eventinfo
        self.particles = particles
        self.weights = weights
        self.attributes = attributes


class LHEEventInfo:
    fieldnames=["nparticles","pid", "weight", "scale", "aqed", "aqcd"]
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self,k,v)
    @classmethod
    def fromstring(cls, string):
        return cls(**dict(zip(cls.fieldnames, map(float,string.split()))))
class LHEParticle(object):
    def __init__(self, **kwargs):
        self.px = 0.0
        self.py = 0.0
        self.pz = 0.0
        self.e = 0.0
        self.m = 0.0
        self.id= 0
        self.__dict__.update(**kwargs)
        self.pdgid=int(self.id)
        self.fourvec = vector.obj(px=self.px,
                                  py=self.py,
                                  pz=self.pz,
                                  E=self.e)
        self.pt=self.fourvec.pt
        self.phi=self.fourvec.phi
        self.eta=self.fourvec.eta
        # self.pdgid_string, self.pdgid_latex =pdgid_to_string(self.pdgid)

    fieldnames = [ "id",
        "status",
        "mother1",
        "mother2",
        "color1",
        "color2",
        "px",
        "py",
        "pz",
        "e",
        "m",
        "lifetime",
        "spin",
    ]

    @classmethod
    def fromstring(cls,string):
        return cls(**dict(zip(cls.fieldnames, map(float, string.split()))))

    def __str__(self):
        return "Particle, PDGID{0}".format( self.pdgid)

def read_lhe(filepath):
    weightdict={}
    for event, element in ET.iterparse(filepath, events=["end"]):
        if element.tag == "initrwgt":
            for initrwgtel in element:
                if initrwgtel.tag == "weightgroup":
                    for weightgroupel in initrwgtel:
                        if weightgroupel.tag=="weight":
                            id = str(weightgroupel.attrib["id"])
                            weightdict[id] = str(weightgroupel.text).split(' #')[0].split(' ')[-2:]
        if element.tag == "event":
            ## here we're not extracting the info block
            eventdict={}
            data = element.text.split("\n")[1:-1]
            eventdata, particles = data[0], data[1:]
            # extracting event info 
            eventdict["eventinfo"] = LHEEventInfo.fromstring(eventdata)
            # extracting weights and attributes
            eventdict["weights"] = {}
            # eventdict["attributes"] = element.attrib
            eventdict["particles"] = []
            for p in particles:
                if not p.strip().startswith("#"):
                    eventdict["particles"] += [LHEParticle.fromstring(p)]
            for sub in element:
                if sub.tag =="rwgt":
                    for r in sub:
                        if r.tag=="wgt":
                            eventdict["weights"][r.attrib["id"]]=float(r.text.strip())
            yield LHEEvent(
                    eventinfo = eventdict["eventinfo"],
                    particles=eventdict["particles"],
                    weights = eventdict["weights"],
                    weightinfo=weightdict,
                    # attributes=eventdict["attributes"]
                    )

def tohdf5(data, filename, key, limit_events=False):
    events = [d for d in data]
    eventinfo= [e.eventinfo for e in events]
    particles = [e.particles for e in events]
    weights = [e.weights for e in events]
    weightinfo = [e.weightinfo for e in events]
    if limit_events:
        df = pd.DataFrame({'event_info':eventinfo[:int(len(events)*0.1)],
                           'particles':particles[:int(len(events)*0.1)],
                           'weights':weights[:int(len(events)*0.1)],
                           })
        df.to_hdf(f"{filename}.h5", key=f"{key}")
    else:
        df = pd.DataFrame({'event_info':eventinfo, 'particles':particles,
                           'weights':weights})
        df.to_hdf(f"{filename}.h5", key=f"{key}")

def extractparams(data, filename, key): 
    events = [d for d in data]
    # eventinfo = [e.eventinfo for e in events]
    particles= [e.particles for e in events]
    weights = [e.weights for e in events]
    # weightinfo = [e.weightinfo for e in events]
    df = pd.DataFrame({'particles':particles, 'weights':weights})
    df['pt_z'] = df.apply(lambda r: ptot(r['particles'],23), axis=1)
    ## extracting eta(Z)
    df['eta_z'] = df.apply(lambda r: eta(r['particles'],23), axis=1)
    ## calc delta phi from two leptons from the Z, in this case the mu+ and mu-
    df['deltaphi_ll_Z'] = df.apply(lambda r: deltaphi(r['particles'], 13, -13), axis=1)

    # t's
    ## extracting pt(t)
    df['pt_t'] = df.apply(lambda r: ptot(r['particles'],6 ), axis=1)
    ## extracting eta(t)
    df['eta_t'] = df.apply(lambda r: eta(r['particles'],6), axis=1)

    # t~'s
    ## extracting pt(t~)
    df['pt_tbar'] = df.apply(lambda r: ptot(r['particles'],-6 ), axis=1)
    ## extracting eta(t~)
    df['eta_tbar'] = df.apply(lambda r: eta(r['particles'],-6), axis=1)

    # deltaR
    df['dR_t_z'] = df.apply(lambda r: deltaR(r['particles'], 6, 23), axis=1)

    # cos theta star z
    df['cosstar'] = df.apply(lambda r: cosstarzlep(r['particles']), axis=1)
    df2 = df.drop('particles', axis=1)
    # df.to_hdf(f"{filename}.h5", key=f"{key}")
    df2.to_hdf(f"{filename}.h5", key=f"{key}")

def ptot(particles, particle_pdgid):
    '''
    determine total transverse momentum of given particle, Z by default
    use with apply function for dataframes
    kwargs: - particles expects a pd dataframe column
            - particle_pdgid is the pdgid of the individual particle
    '''
    for p in particles: 
        if abs(p.pdgid) == particle_pdgid:
            return p.fourvec.pt

def event_weight(events):
    '''
    return weights from event objects when given column of events
'''
    for event in events: 
        return event.weight

def eta(particles, particle_pdgid):
    '''
    determine eta of given particle
    use with apply function for dataframes
    kwargs: - particles expects a pd dataframe column
            - particle_pdgid is the pdgid of the individual particle
    '''
    for p in particles: 
        if abs(p.pdgid) == particle_pdgid:
            return p.fourvec.eta

def deltaphi(particles, pdgid1, pdgid2):
    '''
    determine the difference in phi between two given particles, identified by their pdgids
    kwargs: - particles: a pd dataframe column
            - pdgid1: particle 1
            - pdgid2: particle 2
    '''
    particle_list=[]
    for p in particles:
        if p.id==pdgid1 or p.id==pdgid2:
            particle_list.append(p)
    return particle_list[0].fourvec.deltaphi(particle_list[1].fourvec)

def deltaeta(particles, pdgid1, pdgid2):
    '''
    determine the difference in eta between two given particles, identified by their pdgids
    kwargs: - particles: a pd dataframe column
            - pdgid1: particle 1
            - pdgid2: particle 2
    '''
    particle_list=[]
    for p in particles:
        if p.id==pdgid1 or p.id==pdgid2:
            particle_list.append(p)
    return particle_list[0].fourvec.deltaeta(particle_list[1].fourvec)

def deltaR(particles, pdgid1, pdgid2):
    '''
    determine the difference in eta between two given particles, identified by their pdgids
    kwargs: - particles: a pd dataframe column
            - pdgid1: particle 1
            - pdgid2: particle 2
    '''
    particle_list=[]
    for p in particles:
        if p.id==pdgid1 or p.id==pdgid2:
            particle_list.append(p)
    return particle_list[0].fourvec.deltaR(particle_list[1].fourvec)

def listparticles(particles): 
    '''
    takes the first row of a dataframe and outputs an array of pdgids for _all_ involved particles
    has it's flaws but often useful for sanity checks
    '''
    all_pdgids = []
    for particle in particles[0]:
        if particle.pdgid not in all_pdgids:
            all_pdgids.append(particle.pdgid)
    return all_pdgids

def particlebypdgid(particles, pdgid):
    '''
    given a list of particles and a single pdgid, the vector object of the particle will be returned

    '''
    for p in particles: 
        if p.pdgid == pdgid:
            return p.fourvec


def cosstarzlep(particles):
    '''
    the cosine of the angle
    between the direction of the Z boson in the detector reference 
    frame, and the direction of the negatively-charged lepton from
    the Z boson decay in the rest frame of the Z boson

    to do this we need the Z fourvec
    identify -ve lepton (+ pdgid bc leptons are -ve)
    apply boost_p4(four_vector): change coordinate system 
    using another 4D vector as the difference
    typically apply the negative 4 vec?
    '''

    for p in particles: 
        if p.pdgid == 23:
            z = p.fourvec
        elif p.pdgid == 13: 
            mu_p = p.fourvec

    mu_p_boost = mu_p.boost_p4(z)
    return np.cos(z.deltaangle(mu_p_boost))


# def pdgid_to_string(pdgid):
#     stream = pkg_resources.open_text(__package__, 'pdgid_string.csv')
#     pdgid_data = pd.read_csv(stream)
#     pdgid_data=pdgid_data.set_index('ID')
#     return pdgid_data.loc[pdgid]['Name'], pdgid_data.loc[pdgid]['Latex']


