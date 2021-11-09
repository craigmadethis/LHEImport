import xml.etree.ElementTree as ET
import math

class LHEEvent:
    def __init__(self, eventinfo, particles, weights=None, attributes=None):
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


def read_lhe(filepath):
    for event, element in ET.iterparse(filepath, events=["end"]):
        if element.tag == "event":
            eventdict={}
            ## here we're not extracting the info block
            data = element.text.split("\n")[1:-1]
            eventdata, particles = data[0], data[1:]
            eventdict["eventinfo"] = LHEEventInfo.fromstring(eventdata)
            eventdict["particles"] = []
            # eventdict["weights"] = {}
            # eventdict["attributes"] = element.attrib
            for p in particles:
                if not p.strip().startswith("#"):
                    eventdict["particles"] += [Particle.fromstring(p)]
            yield LHEEvent(eventinfo = eventdict["eventinfo"],
                           particles=eventdict["particles"],
                           # weights = eventdict["weights"],
                           # attributes=eventdict["attributes"]
                           )


class Particle(object):
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
    def __init__(self, pdgid=0, status=0,
                 # initial_state=False,
                 # final_state=False,
                 **kwargs):
        self.pdgid=int(pdgid)
        # build this function (guessing it takes pdgids and checks dict vals to give particle name)
        # self.name=pdgid_to_string(self.pdgid)
        # self.name=self.pdgid
        self.status=int(status)
        # self.final_state = final_state
        # self.initial_state = initial_state
        for k in ['pt', 'eta', 'phi', 'px', 'py', 'pz', 'energy', 'mass']:
        # for k in ['px', 'py', 'pz', 'energy', 'mass']:
            # __dict__ a method of adding params to a class
            self.__dict__[k] = 0.0
        self.__dict__.update(**kwargs)
        if all([k in kwargs for k in ['px', 'py','pz']]):
            pt, eta, phi = self.convert_px_py_pz(float(self.px), float(self.py), float(self.pz))
            self.__dict__['pt'] = pt
            self.__dict__['eta'] = eta
            self.__dict__['phi'] = phi
    ## when you print the object this is what is returned
    # def __str__(self):
    #     return "Particle {0}, PDGID{1}".format(self.barcode, self.pdgid)

    @classmethod
    def fromstring(cls,string):
        return cls(**dict(zip(cls.fieldnames, map(float, string.split()))))

    @classmethod
    def convert_px_py_pz(cls,px,py,pz):
        # transverse momentum
        pt = math.sqrt(math.fsum([math.pow(px, 2), math.pow(py, 2)]))
        # total momentum
        ptot = math.sqrt(math.fsum([math.pow(pt, 2), math.pow(pz, 2)]))
        if pt != 0:
            eta = math.asinh(pz / pt)
            phi = math.asin(py / pt)
        else:
            eta = math.copysign(float('inf'), pz)
            phi = 0
        dphi = 
        return pt, eta, phi

