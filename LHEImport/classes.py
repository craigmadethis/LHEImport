from LHEImport.defs import convert_px_py_pz
class Particle(object):
    
    def __init__(self, barcode, pdgid=0, status=0,
                 # initial_state=False,
                 # final_state=False,
                 **kwargs):
    
        self.barcode=int(barcode)
        self.pdgid=int(pdgid)
        # build this function (guessing it takes pdgids and checks dict vals to give particle name)
        # self.name=pdgid_to_string(self.pdgid)
        # self.name=self.pdgid
        self.status=int(status)
        # self.final_state = final_state
        # self.initial_state = initial_state
    
        for k in ['pt', 'eta', 'phi', 'px', 'py', 'pz', 'energy', 'mass']:
            # __dict__ a method of adding params to a class
            self.__dict__[k] = 0.0
        self.__dict__.update(**kwargs)
        if all([k in kwargs for k in ['px', 'py','pz']]):
            pt, eta, phi = convert_px_py_pz(float(self.px), float(self.py), float(self.pz))
            self.__dict__['pt'] = pt
            self.__dict__['eta'] = eta
            self.__dict__['phi'] = phi
    
    ## when you print the object this is what is returned
    def __str__(self):
        return "Particle {0}, PDGID{1}".format(self.barcode, self.pdgid)
