import os
from pyqtgraph import configfile
from neuroanalysis.data.loaders.loaders import ExperimentLoader


class Acq4ExperimentLoader(ExperimentLoader):
    
    @classmethod
    def get_site_info(site_path):
        index = os.path.join(site_path, '.index')
        if not os.path.isfile(index):
           return 
        return configfile.readConfigFile(index)['.']

    @classmethod
    def get_slice_info(site_path):
        index = os.path.join(site_path, '..', '.index')
        if not os.path.isfile(index):
           return 
        return configfile.readConfigFile(index)['.']

    @classmethod
    def get_expt_info(site_path):
        index = os.path.join(site_path,'../..' '.index')
        if not os.path.isfile(index):
           return 
        return configfile.readConfigFile(index)['.']