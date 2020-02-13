

class DatasetLoader():
    """An abstract base class for Dataset loaders."""
    
    def get_dataset_name(self):
        """Return a string with the name of this dataset."""
        raise NotImplementedError("Must be implemented in subclass.")

    def get_sync_recordings(self, dataset):
        """Return a tuple (list of SyncRecordings, list of RecordingSequences)."""
        raise NotImplementedError("Must be implemented in subclass.")

    def get_recordings(self, sync_rec):
        """Return a dict of {device: Recording}"""
        raise NotImplementedError("Must be implemented in subclass.")

    def get_tseries_data(self, tseries):
        """Return a numpy array of the data in the tseries."""
        raise NotImplementedError("Must be implemented in subclass.")

    def load_stimulus(self, recording):
        """Return an instance of stimuli.Stimulus"""
        raise NotImplementedError("Must be implemented in subclass.")

    def load_stimulus_items(self, recording):
        """Return a list of Stimulus instances. 
        Used with LazyLoadStimulus to parse stimuli when they are needed."""
        raise NotImplementedError("Must be implemented in subclass.")

    def load_test_pulse(self, recording):
        """Return a PatchClampTestPulse."""
        raise NotImplementedError("Must be implemented in subclass.")

    def find_nearest_test_pulse(self, recording):
        raise NotImplementedError("Must be implemented in subclass.")

    def get_baseline_regions(self, recording):
        raise NotImplementedError("Must be implemented in subclass.")


class ExperimentLoader():
    """An abstract base class for Experiment loaders."""

    def __init__(self):
        self._expt = None ### need loader to have access to experiment because the experiment must be handed to Cell objects when they are created.

    def load(self, meta_info=None):
        """Return a tuple of (electrodes, cells, pairs), where each element is an 
        OrderedDict of {electrode_id:Electrode}, {cell_id:Cell} or {pair_id:Pair}.
        meta_info - an optional dict of meta_info to be attached to electrodes, cells or pairs."""
        raise NotImplementedError('Must be implemented in subclass')

    @property
    def expt(self):
        """Return the experiment object this loader loads data for."""
        if self._expt is None:
            raise Exception('No experiment set. During Experiment __init__, please call Loader.set_expt(expt)')
        return self._expt

    def set_expt(self, expt):
        self._expt = expt

    def get_site_info(self):
        """Return a dictionary of meta_info about the site."""
        raise NotImplementedError('Must be implemented in subclass')

    def get_slice_info(self):
        """Return a dictionary of meta_info about the slice."""
        raise NotImplementedError('Must be implemented in subclass')

    def get_expt_info(self):
        """Return a dictionary of meta_info about the experiment."""
        raise NotImplementedError('Must be implemented in subclass')
        






