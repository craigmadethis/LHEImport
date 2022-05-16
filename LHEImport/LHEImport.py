import xml.etree.ElementTree as ET
# import math
import vector
import pandas as pd
import numpy as np
import os 
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.model_selection import train_test_split
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

    fieldnames = [ 
        "id",
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

def gen_weight_dict(filepath):
    '''
    provide an lhe file and generate a {coeff: value} dict for the rwgt scheme
    returns: dictionary of weighting coefficents, {coeff:val}
    '''
    weightdict={}
    for event, element in ET.iterparse(filepath, events=["end"]):
        if element.tag == "initrwgt":
            for initrwgtel in element:
                if initrwgtel.tag == "weightgroup":
                    for weightgroupel in initrwgtel:
                        if weightgroupel.tag=="weight":
                            id = str(weightgroupel.attrib["id"])
                            if weightgroupel.text != None:
                                split = str(weightgroupel.text).split(' #')[0].split(' ')[-2:]
                                # weightdict[id] = str(weightgroupel.text).split(' #')[0].split(' ')[-2:]
                                coeff = split[0]
                                val = split[1]
                                weightdict[id] = {coeff: val}
                            else:
                                weightdict[id] = {None: None}
            else:
                break # break here as we don't want to loop over the entire file
    return weightdict

def read_events(filepath):
    '''
    input: filepath of lhe file
    returns: pd.Dataframe with cols: [eventinfo, particles and weights], each row is an individual event
    '''
    events =[]
    for event, element in ET.iterparse(filepath, events=["end"]):
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
            events.append([
                    eventdict["eventinfo"],
                    eventdict["particles"],
                    eventdict["weights"],
                    ])
    params = ['eventinfo', 'particles', 'weights']
    df = pd.DataFrame(events, columns=params)
    return df

def read_weights(filepath):
    '''
    input: filepath of lhe file
    returns: pd.Dataframe with cols: [eventinfo, particles and weights], each row is an individual event
    '''
    events =[]
    for event, element in ET.iterparse(filepath, events=["end"]):
        if element.tag == "event":
            eventdict={}
            ## here we're not extracting the info block
            # extracting event info 
            # extracting weights and attributes
            eventdict["weights"] = {}
            # eventdict["attributes"] = element.attrib
            for sub in element:
                if sub.tag =="rwgt":
                    for r in sub:
                        if r.tag=="wgt":
                            eventdict["weights"][r.attrib["id"]]=float(r.text.strip())
            events.append([
                    eventdict["weights"],
                    ])
    params = ['weights']
    df = pd.DataFrame(events, columns=params)
    return df



def extract_params(df): 
    '''
    accept a dataframe with eventinfo, particles and weights
    return dataframe with predefined parameters specific to current project
    '''
    df['pt_z'] = df.apply(lambda r: ptot(r['particles'],23), axis=1)
    ## extracting eta(Z)
    df['eta_z'] = df.apply(lambda r: eta(r['particles'],23), axis=1)
    df['phi_z'] = df.apply(lambda r: phi(r['particles'],23), axis=1)
    ## calc delta phi from two leptons from the Z, in this case the mu+ and mu-
    df['deltaphi_ll_Z'] = df.apply(lambda r: deltaphi(r['particles'], 13, -13), axis=1)

    # t's
    ## extracting pt(t)
    df['pt_t'] = df.apply(lambda r: ptot(r['particles'],6 ), axis=1)
    ## extracting eta(t)
    df['eta_t'] = df.apply(lambda r: eta(r['particles'],6), axis=1)
    df['phi_t'] = df.apply(lambda r: phi(r['particles'],6), axis=1)

    # t~'s
    ## extracting pt(t~)
    df['pt_tbar'] = df.apply(lambda r: ptot(r['particles'],-6 ), axis=1)
    ## extracting eta(t~)
    df['eta_tbar'] = df.apply(lambda r: eta(r['particles'],-6), axis=1)
    df['phi_tbar'] = df.apply(lambda r: phi(r['particles'],-6), axis=1)

    # mu's
    ## extracting pt(mu)
    df['pt_mu'] = df.apply(lambda r: ptot(r['particles'],13 ), axis=1)
    ## extracting eta(mu)
    df['eta_mu'] = df.apply(lambda r: eta(r['particles'],13), axis=1)
    df['phi_mu'] = df.apply(lambda r: phi(r['particles'],13), axis=1)

    # mu~'s
    ## extracting pt(mu~)
    df['pt_mubar'] = df.apply(lambda r: ptot(r['particles'],-13 ), axis=1)
    ## extracting eta(mu~)
    df['eta_mubar'] = df.apply(lambda r: eta(r['particles'],-13), axis=1)
    df['phi_mubar'] = df.apply(lambda r: phi(r['particles'],-13), axis=1)

    # w's
    ## extracting pt(w)
    df['pt_w'] = df.apply(lambda r: ptot(r['particles'],24 ), axis=1)
    df['mt_w'] = df.apply(lambda r: mT(r['particles'],24 ), axis=1)
    ## extracting eta(w)
    df['eta_w'] = df.apply(lambda r: eta(r['particles'],24), axis=1)
    df['phi_w'] = df.apply(lambda r: phi(r['particles'],24), axis=1)

    # w~'s
    ## extracting pt(w~)
    df['pt_wbar'] = df.apply(lambda r: ptot(r['particles'],-24 ), axis=1)
    ## extracting eta(w~)
    df['eta_wbar'] = df.apply(lambda r: eta(r['particles'],-24), axis=1)
    df['phi_wbar'] = df.apply(lambda r: phi(r['particles'],-24), axis=1)

    # b's
    ## extracting pt(b)
    df['pt_b'] = df.apply(lambda r: ptot(r['particles'],5 ), axis=1)
    ## extracting eta(b)
    df['eta_b'] = df.apply(lambda r: eta(r['particles'],5), axis=1)
    df['phi_b'] = df.apply(lambda r: phi(r['particles'],5), axis=1)

    # b~'s
    ## extracting pt(b~)
    df['pt_bbar'] = df.apply(lambda r: ptot(r['particles'],-5 ), axis=1)
    ## extracting eta(b~)
    df['eta_bbar'] = df.apply(lambda r: eta(r['particles'],-5), axis=1)
    df['phi_bbar'] = df.apply(lambda r: phi(r['particles'],-5), axis=1)

    # deltaR
    df['dR_t_z'] = df.apply(lambda r: deltaR(r['particles'], 6, 23), axis=1)

    # cos theta star z
    df['cosstar'] = df.apply(lambda r: cosstarzlep(r['particles']), axis=1)
    df2 = df.drop(['eventinfo'], axis=1)


    ## could maybe not drop particles, that way I could add features if I wanted?
    # df.to_hdf(f"{filename}.h5", key=f"{key}")
    return df2


def extract_params2(df): 
    '''
    accept a dataframe with eventinfo, particles and weights
    return dataframe with predefined parameters specific to current project
    '''
    df['pt_mumu'] = df.apply(lambda r: ptot(r['particles'], [13,-13]), axis=1)
    ## extracting eta(Z)
    df['eta_mumu'] = df.apply(lambda r: eta(r['particles'],[13,-13]), axis=1)
    df['phi_mumu'] = df.apply(lambda r: phi(r['particles'],[13,-13]), axis=1)
    ## calc delta phi from two leptons from the Z, in this case the mu+ and mu-
    df['deltaphi_mumu'] = df.apply(lambda r: deltaphi(r['particles'], 13, -13), axis=1)
    df['deltaeta_mumu'] = df.apply(lambda r: deltaeta(r['particles'], 13, -13), axis=1)
    df['s_mumu'] = df.apply(lambda r: invmass(r['particles'], [13, -13]), axis=1)

    # t's
    ## extracting pt(t)
    df['pt_Wb'] = df.apply(lambda r: ptot(r['particles'],[24, 5] ), axis=1)
    df['eta_Wb'] = df.apply(lambda r: eta(r['particles'],[24, 5] ), axis=1)
    df['phi_Wb'] = df.apply(lambda r: phi(r['particles'],[24, 5] ), axis=1)
    df['deltaphi_wb'] = df.apply(lambda r: deltaphi(r['particles'], 24, 5), axis=1)
    df['deltaeta_wb'] = df.apply(lambda r: deltaeta(r['particles'], 24, 5), axis=1)
    df['mass_wb'] = df.apply(lambda r: mass(r['particles'], [24, 5]), axis=1)
    df['s_wb'] = df.apply(lambda r: invmass(r['particles'], [24, 5]), axis=1)
    ## extracting eta(t)

    # t~'s
    ## extracting pt(t~)
    df['pt_Wb-'] = df.apply(lambda r: ptot(r['particles'],[-24, -5] ), axis=1)
    df['eta_Wb-'] = df.apply(lambda r: eta(r['particles'],[-24, -5] ), axis=1)
    df['phi_Wb-'] = df.apply(lambda r: phi(r['particles'],[-24, -5] ), axis=1)
    df['deltaphi_wb-'] = df.apply(lambda r: deltaphi(r['particles'], -24, -5), axis=1)
    df['deltaeta_wb-'] = df.apply(lambda r: deltaeta(r['particles'], -24, -5), axis=1)
    df['s_wb-'] = df.apply(lambda r: invmass(r['particles'], [-24, -5]), axis=1)

    # mu's
    ## extracting pt(mu)
    df['pt_mu'] = df.apply(lambda r: ptot(r['particles'],[13] ), axis=1)
    ## extracting eta(mu)
    df['eta_mu'] = df.apply(lambda r: eta(r['particles'],[13]), axis=1)
    df['phi_mu'] = df.apply(lambda r: phi(r['particles'],[13]), axis=1)

    # mu~'s
    ## extracting pt(mu~)
    df['pt_mubar'] = df.apply(lambda r: ptot(r['particles'],[-13] ), axis=1)
    ## extracting eta(mu~)
    df['eta_mubar'] = df.apply(lambda r: eta(r['particles'],[-13]), axis=1)
    df['phi_mubar'] = df.apply(lambda r: phi(r['particles'],[-13]), axis=1)

    # w's
    ## extracting pt(w)
    df['pt_w'] = df.apply(lambda r: ptot(r['particles'],[24] ), axis=1)
    df['mt_w'] = df.apply(lambda r: mT(r['particles'],[24]), axis=1)
    ## extracting eta(w)
    df['eta_w'] = df.apply(lambda r: eta(r['particles'],[24]), axis=1)
    df['phi_w'] = df.apply(lambda r: phi(r['particles'],[24]), axis=1)

    # w~'s
    ## extracting pt(w~)
    df['pt_wbar'] = df.apply(lambda r: ptot(r['particles'],[-24]), axis=1)
    ## extracting eta(w~)
    df['eta_wbar'] = df.apply(lambda r: eta(r['particles'],[-24]), axis=1)
    df['phi_wbar'] = df.apply(lambda r: phi(r['particles'],[-24]), axis=1)

    # b's
    ## extracting pt(b)
    df['pt_b'] = df.apply(lambda r: ptot(r['particles'],[5]), axis=1)
    ## extracting eta(b)
    df['eta_b'] = df.apply(lambda r: eta(r['particles'],[5]), axis=1)
    df['phi_b'] = df.apply(lambda r: phi(r['particles'],[5]), axis=1)

    # b~'s
    ## extracting pt(b~)
    df['pt_bbar'] = df.apply(lambda r: ptot(r['particles'],[-5] ), axis=1)
    ## extracting eta(b~)
    df['eta_bbar'] = df.apply(lambda r: eta(r['particles'],[-5]), axis=1)
    df['phi_bbar'] = df.apply(lambda r: phi(r['particles'],[-5]), axis=1)

    # deltaR
    df['dR_t_z'] = df.apply(lambda r: deltaR(r['particles'], 6, 23), axis=1)
    df['dR_mumu'] = df.apply(lambda r: deltaR(r['particles'], 13, -13), axis=1)
    df['dR_wb'] = df.apply(lambda r: deltaR(r['particles'], 24, 5), axis=1)
    df['dR_wb-'] = df.apply(lambda r: deltaR(r['particles'], -24, -5), axis=1)

    # cos theta star z
    df['cosstar'] = df.apply(lambda r: cosstarzlep(r['particles']), axis=1)
    df2 = df.drop(['eventinfo', 'weights', 'particles'], axis=1)


    ## could maybe not drop particles, that way I could add features if I wanted?
    # df.to_hdf(f"{filename}.h5", key=f"{key}")
    return df2




def ptot(particles: list, particle_pdgid: list):
    '''
    determine total transverse momentum of given particle, Z by default
    use with apply function for dataframes
    kwargs: - particles expects a pd dataframe column
            - particle_pdgid is an array with pdgids of particles
    '''
    p_vec = vector.obj(rho=0,eta = 0,phi=0,e = 0)
    for p in particles: 
        if p.pdgid in particle_pdgid:
            p_vec += p.fourvec
    return p_vec.pt

def invmass(particles: list, particle_pdgid: list):
    p_vec = vector.obj(rho=0, eta=0, phi=0, e=0)
    for p in particles: 
        if p.pdgid in particle_pdgid:
            p_vec += p.fourvec
    return np.sqrt(p_vec.E2 + p_vec.mag2)

def mass(particles: list, particle_pdgid: list):
    '''
    determine total transverse momentum of given particle, Z by default
    use with apply function for dataframes
    kwargs: - particles expects a pd dataframe column
            - particle_pdgid is an array with pdgids of particles
    '''
    p_vec = vector.obj(rho=0,eta = 0,phi=0,e = 0)
    for p in particles: 
        if p.pdgid in particle_pdgid:
            p_vec += p.fourvec
    return p_vec.m

def event_weight(events: list):
    '''
    return weights from event objects when given column of events
'''
    for event in events: 
        return event.weight

def eta(particles: list, particle_pdgid: list):
    '''
    determine eta of given particle
    use with apply function for dataframes
    kwargs: - particles expects a pd dataframe column
            - particle_pdgid is the pdgid of the individual particle
    '''
    p_vec = vector.obj(rho=0,eta = 0,phi=0,e = 0)
    for p in particles: 
        if p.pdgid in particle_pdgid:
            p_vec += p.fourvec
    return p_vec.eta

def mT(particles: list, particle_pdgid: list):
    '''
    determine eta of given particle
    use with apply function for dataframes
    kwargs: - particles expects a pd dataframe column
            - particle_pdgid is the pdgid of the individual particle
    '''
    p_vec = vector.obj(rho=0,eta = 0,phi=0,e = 0)
    for p in particles: 
        if p.pdgid in particle_pdgid:
            p_vec += p.fourvec
    return p_vec.mt

def phi(particles: list, particle_pdgid: list):
    '''
    determine eta of given particle
    use with apply function for dataframes
    kwargs: - particles expects a pd dataframe column
            - particle_pdgid is the pdgid of the individual particle
    '''
    p_vec = vector.obj(rho=0,eta = 0,phi=0,e = 0)
    for p in particles: 
        if p.pdgid in particle_pdgid:
            p_vec += p.fourvec
    return p_vec.phi

def deltaphi(particles: list, pdgid1: int, pdgid2: int):
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

def deltaeta(particles: list, pdgid1: int, pdgid2: int):
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

def deltaR(particles: list, pdgid1: int, pdgid2: int):
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

def listparticles(particles: list): 
    '''
    takes the first row of a dataframe and outputs an array of pdgids for _all_ involved particles
    has it's flaws but often useful for sanity checks
    '''
    all_pdgids = []
    for event in particles:
        for particle in event:
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


def cosstarzlep(particles: list):
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

def getSample(df, samplesize: int, weights: list, label: str):
    '''
    from provided dataframe, fetch a sample of size samplesize given a set of weights
    df: pd.dataframe(data)
    samplesize: int, size of sample required
    weights: array of size df
    label: label to give each of the values
    '''
    pseudo_idx = np.random.choice(df.shape[0], size=samplesize, replace=True, p = weights/weights.sum() )
    dff = df.iloc[pseudo_idx].copy(deep=True)
    dff.reset_index()
    dff['label'] = np.repeat(label, len(dff))
    return dff

def create_weight_df(keys, weights):
    '''
    Given keys relating to rwgt id, a dataframe of event weights are generated for each eft coefficient. 
    keys: keys from a weight dictionary, e.g. weight_dict.keys()
    weights: weight column of a df
    '''
    df = pd.DataFrame()
    for key in keys:
        df[f'{key}'] = np.array([d.get(key) for d in weights])
    return df

def filter_dataset(sample,IQR_percent=2):
        '''
        here just cleaning up the pt outliers, also cleaning the eta outliers as there were some extreme outliers. 
        '''
        variables = sample.columns
        for var in np.asarray(variables):
            if var!= 'label':
                Q1 = sample[str(var)].quantile(0.25)
                Q3 = sample[(var)].quantile(0.75)
                IQR = Q3 - Q1
                upper = Q3 + IQR_percent*IQR
                lower = Q1 - IQR_percent*IQR
                sample = sample.loc[(sample[var] >= lower) & (sample[var] <= upper) & (sample[var] != np.inf) & (sample[var] != np.nan)]
        return sample

def minmax(data, new_min=0, new_max=1):
    '''
    For a given dataseries, will rescale setting the lowest value in the series to 0 and greatest to 1
    '''
    return (data-data.min())/(data.max()- data.min()) * (new_max - new_min) + new_min


def rescale_df(data, columns):
    X_scaled = pd.DataFrame()
    for column in np.array(columns):
        df = pd.DataFrame({column: minmax(data[column])})
        X_scaled = pd.concat([df, X_scaled], axis=1)
    X_scaled = X_scaled[data.columns]
    return X_scaled 

def train_model(df, weights, eft_label, sample_size, num_layers =3, num_nodes=256, num_epochs=50, batchsize=512, earlystopping=True, save_model=False, model_dir='', model_id=1):
    '''
    df: DataFrame containing all parameters
    eft_label: rwgt_22 etc
    ''' 
    dataset = get_sample_set(df, sample_size, weights, eft_label)
    # dataset['eft'] = dataset.apply(lambda r: 0 if (r['label'] == 'rwgt_1') or (r['label'] == 0) else 1, axis=1)
    
    print('Samples generated')
    
    X = dataset.copy()
    y = X.pop('eft')
    labels = X.pop('label')
    
    #rescaling
    X_scaled = rescale_df(X, X.columns)
    
    #train test split
    train_X, val_X, train_y, val_y = train_test_split(X_scaled, y, test_size = 0.2)
    input_shape = [train_X.shape[1]]
    
    # model creation
    model = keras.Sequential()
    model.add(layers.BatchNormalization(input_shape=input_shape))
    for _ in range(num_layers):
        model.add(layers.Dense(num_nodes, activation='relu'))
        model.add(layers.Dropout(0.3))
    model.add(layers.Dense(1, activation='sigmoid'))
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['binary_accuracy'])

    early_stopping = keras.callbacks.EarlyStopping(
            monitor='val_loss',
            patience = 0.5*num_epochs,
            min_delta=0.001, 
            restore_best_weights=True,
            )

# history = model.fit(train_X, train_y, validation_data=(val_X, val_y), batch_size=512, epochs=30, callbacks=[early_stopping])
    print(f"Training model {model_id} with {num_layers} layers and {num_nodes} nodes for {num_epochs} epochs")
    history = model.fit(
            train_X, 
            train_y, 
            validation_data=(val_X, val_y), 
            batch_size=batchsize, 
            epochs=num_epochs, 
            callbacks=[early_stopping], 
            verbose=0,
            )
    history_df = pd.DataFrame(history.history)
    if save_model == True:
        if not os.path.exists(model_dir+"/models"):
            os.makedirs(model_dir+"/models/")
        if not os.path.exists(model_dir+"/history"):
            os.makedirs(model_dir+"/history")
        model.save(model_dir+f"models/{num_layers}_{num_nodes}_{model_id}.h5")
        model.save(model_dir+f"models/{num_layers}_{num_nodes}_{model_id}.h5")
        history_df.to_hdf(model_dir+f"history/{num_layers}_{num_nodes}_{model_id}.h5", key="all")
        print('Model and history saved!')
    return model, history_df


def eval_models(df, weights, model_dir, samples, weight_label, roc=False):
    aucs_all = []
    rocs = []
    modelnames =[]
    modelfiles = []
    for file in os.listdir(model_dir):
        if ".h5" in file:
            modelfiles.append(file)
    for modelname in modelfiles:
        # sm_sample = getSample(df, samples, weights['rwgt_1'], 'rwgt_1')
        # eft_sample = getSample(df, samples, weights[weight_label], weight_label)

        # sm_filtered = filter_dataset(sm_sample)
        # eft_filtered = filter_dataset(eft_sample)
        # sm_filtered['eft'] = np.repeat(0, len(sm_filtered))
        # eft_filtered['eft'] = np.repeat(1, len(eft_filtered))
        # dataset = pd.concat([sm_filtered, eft_filtered])
        # print(dataset['label'].describe())
        auc_model = []
        roc_model = []
        modelnames.append(model_dir+modelname)
        for weight in weight_label:
            eval_weights = ['rwgt_1', weight]
            if weight != 'rwgt_1':
                dataset = get_sample_set(df, samples, weights, eval_weights)
                X = dataset.copy()
                y_testing = X.pop('eft')
                labels = X.pop('label')
                                     
                X_testing = rescale_df(X, X.columns)
                print('Samples generated')
            
                print('loading ', model_dir+modelname)
                model = keras.models.load_model(model_dir+modelname)
                print(model_dir+modelname+' loaded')
                # print(model.summary())
                print('Making predictions...')
                y_preds = model.predict(X_testing)
                y_true = y_testing
                roc_stats = roc_curve(y_true, y_preds)
                auc_val = auc(roc_stats[0], roc_stats[1])
                auc_model.append(auc_val)
                roc_model.append(roc_stats)
                print('AUC Calculated')
        aucs_all.append(auc_model)
        rocs.append(roc_model)
    aucs_df = pd.DataFrame({'Model': modelnames, 'AUC':aucs_all, 'roc_stats': rocs})
    return aucs_df,rocs

def get_sample_set(df, sample_size, weights, weight_list):
    dataset =pd.DataFrame()
    eft_samples = pd.DataFrame()
    sm_samples = pd.DataFrame()
    for weight in weight_list:
        sample = getSample(df, sample_size, weights[weight], weight)
        sample = filter_dataset(sample)
        if weight == 'rwgt_1':
            # sm_sample = getSample(df, sample_size, weights[weight], weight)
            sample['eft'] = np.repeat(0, len(sample))
            sm_samples = pd.concat([sm_samples, sample])
        else:
            sample['eft'] = np.repeat(1, len(sample))
            eft_sample = sample
            eft_samples = pd.concat([eft_samples, eft_sample])
    if len(weight_list) > 2:
        eft_samples = eft_samples.sample(len(sm_samples))
    dataset = pd.concat([sm_samples,eft_samples])
    # print(dataset)
    return dataset

