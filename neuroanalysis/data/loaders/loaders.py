

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

    def load(self, expt, meta_info=None):
        """Return a tuple of (electrodes, cells, pairs), where each element is an 
        OrderedDict of {electrode_id:Electrode}, {cell_id:Cell} or {pair_id:Pair}.
        expt - the experiment instance -- we need this because we need to pass it to cells when we instantiate them.
        meta_info - an optional dict of meta_info to be attached to electrodes, cells or pairs."""
        raise NotImplementedError('Must be implemented in subclass')

    def verify(self):
        """Perform some basic checks to ensure this experiment was acquired correctly.
        Raise Exception if there is a problem."""
        raise NotImplementedError('Must be implemented in subclass')

    # def find_files(self):
    #     """Return a dict of {'key': path} for important files in this experiment."""
    #     raise NotImplementedError('Must be implemented in subclass')

    # def get_timestamp(self):
    #     """Return a timestamp corresponding to the approximate start of data acquisition for this experiment."""
    #     raise NotImplementedError('Must be implemented in subclass')

    def get_info(self, meta_info=None):
        """Return a dict of meta_info about this experiment.
        meta-info - an optional dict of meta_info that was passed to the Experiment upon initialization"""
        raise NotImplementedError('Must be implemented in subclass')

    def get_ext_id(self):
        """Return a unique identification string for this experiment."""
        raise NotImplementedError('Must be implemented in subclass')

    def get_ephys_data(self):
        """Return a Dataset object containing the ephys data collected with this experiment."""
        raise NotImplementedError('Must be implemented in subclass')

    def get_surface_depth(self):
        """Return the depth of the surface of the slice, or None. Needed for Cells to calculate cell.depth"""
        raise NotImplementedError('Must be implemented in subclass')

    def get_slice(self):
        """Return a Slice object representing the slice used in this experiment."""
        raise NotImplementedError('Must be implemented in subclass')

    def get_datetime(self):
        """Return a datetime object corresponding to the beginning of data acquisition for this experiment."""
        raise NotImplementedError('Must be implemented in subclass')

    def get_target_region(self):
        """Return a string containing the name of the brain region targeted in this experiment."""
        raise NotImplementedError('Must be implemented in subclass')

    def get_target_temperature(self):
        """Return the intended temperature of the experiment in C (as a float), or None."""
        raise NotImplementedError('Must be implemented in subclass')

