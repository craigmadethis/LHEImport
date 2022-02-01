import xml.etree.ElementTree as ET
# import math
import vector
import pandas as pd
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
    eventdict={}
    eventdict["weightinfo"]={}
    for event, element in ET.iterparse(filepath, events=["end"]):
        if element.tag == "initrwgt":
            for initrwgtel in element:
                if initrwgtel.tag == "weightgroup":
                    for weightgroupel in initrwgtel:
                        if weightgroupel.tag=="weight":
                            id = str(weightgroupel.attrib["id"])
                            eventdict["weightinfo"][id] = str(weightgroupel.text).split(' #')[0].split(' ')[-2:]

        elif element.tag == "event":
            # eventdict={}
            ## here we're not extracting the info block
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
    yield LHEEvent(eventinfo = eventdict["eventinfo"],
            particles=eventdict["particles"],
            weights = eventdict["weights"],
            weightinfo=eventdict["weightinfo"],
            # attributes=eventdict["attributes"]
            )

def tohdf5(data, filename, key, limit_events=False):
    events = [d for d in data]
    eventinfo= [e.eventinfo for e in events]
    particles = [e.particles for e in events]
    weights = [e.weights for e in events]
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

# def pdgid_to_string(pdgid):
#     stream = pkg_resources.open_text(__package__, 'pdgid_string.csv')
#     pdgid_data = pd.read_csv(stream)
#     pdgid_data=pdgid_data.set_index('ID')
#     return pdgid_data.loc[pdgid]['Name'], pdgid_data.loc[pdgid]['Latex']
