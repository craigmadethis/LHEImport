from lxml import etree as ET
from LHEImport.classes import Particle
from LHEImport.defs import convert_px_py_pz
class LHEImport(object):
    '''
    main class to parse a LHE file

    user passes event number to return that event (first event = 1)
    if unassigned, return first event in file
    '''

    def __init__(self,filename):
        '''
        params
        ------
        filename: str
            Input filename.
        event_num: int, optional
            index of event to parse
        '''
        self.filename = filename
        # self.barcode = barcode
        self.events = []
        self.root,self.init_ind, self.event_ind, self.num_events = self.readfile()
        self.event_num = self.event_ind
        # self.Particle = self.Particle()

    def __str__(self):
        return "LHEParser: %s".format(self.filename)

    def readfile(self):
        tree = ET.parse(self.filename)
        root = tree.getroot()
        tags = [r.tag for r in root]
        init_ind = tags.index('init')
        event_ind = tags.index('event')
        num_events = tags.count('event')
        return (root,
                init_ind,
                event_ind,
                num_events)

    # def parse(self):
    #     '''
    #     parse contents of input file and extract the particles
    #     Returns
    #     -------
    #     Event
    #         Event object containning info about the event
    #     list[NodeParticle]
    #         Collection of Node particles to be assigned to a graph
    #     '''

    #     # parsing the xml
    #     tree = ET.parse(self.filename)
    #     root = tree.getroot()
    #     tags = [r.tag for r in root]

    #     # determining index of init block
    #     init_ind = tags.index('init')

    #     ## obtain user event
    #     event_ind = init_ind
    #     event_ind = tags.index('event',self.event_num)

    #     init = self.parse_init_line(root[init_ind].text)
    #     event, node_particles = self.parse_event_text(root[event_ind].text)
    #     return event, node_particles, init
    def importevents(self, limit_events=False):
        events = []
        init = self.parse_init_line(self.root[self.init_ind].text)
        if limit_events:
            for i in range(self.event_ind, int((self.num_events+1)*0.01)):
                events.append(self.parse_event_text(self.root[i].text))
        else:
            for i in range(self.event_ind, self.num_events+1):
                events.append( self.parse_event_text(self.root[i].text))
        data = {'stats': init, 'eventdata': events}
        return data


    def map_columns_to_dict(self, fields, line, delim=None):
        '''
        Splits a line into fields, stores them in a dict
        '''
        parts = line.strip().split(delim)[0:len(fields)+1]
        return {k: v.strip() for k,v in zip(fields, parts)}

    def parse_event_text(self, text):
        event = None
        final_particles = []
        counter = 1
        for line in text.split('\n'):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if not event:
                event = self.parse_event_line(line)
            else:
                node_particle = self.parse_particle_line(line, barcode=counter)
                final_particles.append(node_particle)
                counter += 1
        return {'eventstats':event, 'final_particles':final_particles}

    def parse_event_line(self, line):
        fields = ["num_particles", "proc_id", "weight", "scale", "aQED", "aQCD"]
        contents = self.map_columns_to_dict(fields, line)
        return contents

    def parse_init_line(self, line):
        '''
        determined the init fields from package
        https://github.com/scikit-hep/pylhe/blob/master/src/pylhe/__init__.py
        '''
        fields = ["beamA", "beamB", "E BeamA", "E BeamB",
                  "PDFgroupA","PDFgroupB",
                  "PDFsetA","PDFsetB","weighting_strat", "numProcesses",
                  "XSec", "error", "unit weight", "procId"]
        contents = self.map_columns_to_dict(fields,line)
        return contents

    def parse_particle_line(self,line,barcode):
        fields = ["pdgid", "status", "parent1", "parent2", "col1", "col2",
                      "px", "py", "pz", "energy", "mass", "lifetime", "spin"]
        contents_dict = self.map_columns_to_dict(fields,line)
        p = Particle(barcode=barcode,
                     pdgid=int(contents_dict['pdgid']),
                     status=int(contents_dict['status']),
                     px=float(contents_dict['px']),
                     py=float(contents_dict['py']),
                     pz=float(contents_dict['pz']),
                     energy=float(contents_dict['energy']),
                     mass=float(contents_dict['mass']))

        return p







