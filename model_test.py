import os
import h5py
import pyqtgraph as pg

## dataset testing
from neuroanalysis.data.loaders.mies_dataset_loader import MiesNwbLoader
#from neuroanalysis.data.loaders.acq4_dataset_loader import Acq4DatasetLoader
from neuroanalysis.data.dataset import Dataset
from optoanalysis.analyzers import OptoBaselineAnalyzer
from aisynphys.analyzers import MPBaselineAnalyzer
from neuroanalysis.miesnwb import MiesNwb

## experiment testing
from neuroanalysis.data.experiment import AI_Experiment
from neuroanalysis.data.loaders.opto_experiment_loader import OptoExperimentLoader


pg.dbg()


f = "/Users/meganbkratz/Code/ai/example_data/data/2019-06-13_000/slice_000/site_000/2019_06_13_exp1_TH-compressed.nwb"
f2 = "/Users/meganbkratz/Code/ai/example_data/2019_06_24_131623-compressed.nwb"
f3 = "/Users/meganbkratz/Documents/ManisLab/L4Mapping/ExcitationProfileData/2012.11.09_000/slice_000/cell_004"

#hdf = h5py.File(f, 'r')

mies_nwb = Dataset(loader=MiesNwbLoader(f2, baseline_analyzer_class=MPBaselineAnalyzer))
opto_nwb = Dataset(loader=MiesNwbLoader(f, baseline_analyzer_class=OptoBaselineAnalyzer))
#acq4_dataset = Dataset(loader=Acq4DatasetLoader(f3))

expt = AI_Experiment(loader=OptoExperimentLoader(site_path=os.path.split(f)[0]))


