import os, glob, re, math, datetime
import json, yaml
from collections import OrderedDict
from .ai_experiment_loader import AI_ExperimentLoader
from .mies_dataset_loader import MiesNwbLoader
from neuroanalysis.data.cell import Cell
from neuroanalysis.data.electrode import Electrode
from neuroanalysis.data.dataset import Dataset
from neuroanalysis.data.pair import Pair
from neuroanalysis.data.slice import AI_Slice
from pyqtgraph import configfile
from optoanalysis.analyzers import OptoBaselineAnalyzer
from optoanalysis import data_model


class OptoExperimentLoader(AI_ExperimentLoader):
    
    def __init__(self, load_file=None, site_path=None, meta_info=None):
        AI_ExperimentLoader.__init__(self, load_file=load_file, site_path=site_path)

        self._cnx_file = load_file
        self._meta_info = meta_info

    @property
    def cnx_file(self):
        if self._cnx_file is None:
            self._cnx_file = self.find_connections_file()
        return self._cnx_file

    def load(self, expt):
        """Return a tuple of (electrodes, cells, pairs), where each element is an 
        OrderedDict of {electrode_id:Electrode}, {cell_id:Cell} or {pair_id:Pair}.
        Parameters
        ----------
        meta_info - an optional dict of meta_info to be attached to electrodes, cells or pairs.
        """

        if self.cnx_file == "not found":
            electrodes, cells = self.load_without_connection_file(expt)
        else:
            version = self.get_cnx_file_version(self.cnx_file)
            if version == 0:
                electrodes, cells = self.load_markPoints_connection_file(expt)
            else:
                electrodes, cells = self.load_mosaiceditor_connection_file(expt)

        self.process_meta_info(electrodes, cells, self._meta_info)

        pairs = self.create_pairs(cells)

        return (electrodes, cells, pairs)

    def get_ephys_data(self):
        nwb = self.get_ephys_file()
        if nwb is not None:
            #return OptoNwb(nwb)
            return Dataset(loader=MiesNwbLoader(nwb, baseline_analyzer_class=OptoBaselineAnalyzer))

    # def find_files(self):
    #     files = AI_ExperimentLoader.find_files(self)
    #     files['connections'] = self.cnx_file if self.cnx_file is not None else self.find_connections_file()
    #     return files

    def get_ext_id(self):
        if self.cnx_file != 'not found':
            return os.path.split(self.cnx_file)[1].partition('_connections')[0] ### connections file name minus '_connections.json'
        elif self._meta_info is not None:
            return self._meta_info['experiment'].partition('_connections')[0]
        else:
            raise Exception("Can't make name for experiment: (\n\tload_file:%s, \n\tsite_path:%s, \n\tmeta_info:%s)"%(self._cnx_file, self.site_path, self._meta_info))

    def get_surface_depth(self):
        if self.site_path is None:
            return None

        mp_files = sorted(glob.glob(os.path.join(self.site_path, 'MultiPatch_*.log')))

        lines = []
        for mplog in mp_files:
            lines.extend([l for l in open(mplog, 'r').readlines() if 'surface_depth_changed' in l])
        if len(lines) == 0:
            return None
        line = lines[-1].rstrip(',\r\n')
        return json.loads(line)['surface_depth']

    def get_last_modification_time(self):
        if self.site_path is not None:
            files = [
                self.site_path,
                self.cnx_file,
                self.get_mosaic_file(),
                os.path.join(self.site_path, '.index'),
                os.path.join(self.site_path, '..', '.index'),
                os.path.join(self.site_path, '..', '..', '.index'),
                ]
        else:
            files = [
                self.cnx_file
                ]

        mtime = 0
        for f in files:
            if f is None or not os.path.exists(f):
                continue
            mtime = max(mtime, os.stat(f).st_mtime)
        
        return datetime.datetime.fromtimestamp(mtime)

    def get_rig_name(self):
        rig = AI_ExperimentLoader.get_rig_name(self)
        if rig is None:
            data = self.get_ephys_data()
            if data is not None:
                rig = data_model.get_rig_from_nwb(data)

        return rig

    def get_target_region(self):
        #### Should we hardcode in a region, like V1?
        if self._expt.slice is not None and self._expt.slice.lims_record is not None:
            return AI_ExperimentLoader.get_target_region(self)
        else:
            return None

    def get_rig_operator(self):
        op = AI_ExperimentLoader.get_rig_operator(self)
        if op is None:
            return self.get_ext_id().split('_')[-1]

    def get_info(self):
        info = AI_ExperimentLoader.get_info(self)
        info['additional_info'] = self._meta_info
        return info
                    


#### private functions:

    def get_electrodes_from_dataset(self):
        """Return a list of electrode ids found in .nwb file or None if no .nwb file is available."""

        dataset = self.get_ephys_data()
        if dataset is None:
            return

        devs = list(dataset.devices)
        return [x for x in devs if type(x) == type(1)]



    def process_meta_info(self, electrodes, cells, meta_info):
        """Process optional meta_info dict that is passed in at the initialization of expt.
        """
        ## Need to load: presynapticCre, presynapticEffector, [class, reporter, layer for each headstage], internal
        if meta_info is None:
            print('Warning: Loading %s without meta info.' % self._expt)

        preEffector = meta_info.get('presynapticEffector', '').lower()
        for e_id, elec in electrodes.items():
            n = e_id[-1]
            cell = elec.cell
            if cell._morphology.get('initial_call') is None:
                cell._morphology['initial_call'] = meta_info.get('HS%s_class'%n)
            if cell._target_layer is None:
                cell._target_layer = meta_info.get('HS%s_layer'%n).strip().strip('L').strip('l')
            if meta_info.get('HS%s_reporter'%n, '').lower() == 'positive':
                cell._cre_type = meta_info.get('presynapticCre','').lower()
            self.label_cell(cell, preEffector, meta_info.get('HS%s_reporter'%n, '').lower() == 'positive')
            elec._internal_solution = meta_info.get('internal', '').lower()
            if len(meta_info.get('distances', [])) > 0: ### we may not have distance measurements for all cells
                dist = [e for e in meta_info.get('distances') if e.get('headstage')==n]
                if len(dist) == 1:
                    if cell._distance_to_pia is None:
                        try:
                            cell._distance_to_pia = float(dist[0]['toPia'])*1e-6
                        except ValueError:
                            if dist[0]['toPia'] != '':
                                raise
                    if cell._distance_to_wm is None:
                        try:
                            cell._distance_to_wm = float(dist[0]['toWM'])*1e-6
                        except ValueError:
                            if dist[0]['toWM'] != '':
                                raise
                    if dist[0]['photostim_as'] != '': ## this electrode targeted a photostimulated point, need to merge cells
                        name = 'Point '+str(int(dist[0]['photostim_as']))
                        ps_cell = cells.pop(name)
                        old_name = cell.cell_id
                        cell.cell_id = name+'/'+old_name
                        cell.info.update(stim_point_ext_id=name, stim_point_info={
                            'position':ps_cell.position, 
                            'target_layer':ps_cell._target_layer, 
                            'percent_depth':ps_cell._percent_depth,
                            'distance_to_pia':ps_cell._distance_to_pia,
                            'distance_to_wm':ps_cell._distance_to_wm})
                        cells.pop(old_name)
                        cells[cell.cell_id]=cell
                elif len(dist) > 1:
                    raise Exception('Something is wrong.')
            if cell._percent_depth is None and cell._distance_to_pia is not None and cell._distance_to_wm is not None:
                cell._percent_depth = cell._distance_to_pia/(cell._distance_to_pia + cell._distance_to_wm)

        for i, cell in cells.items():
            if not cell.has_readout and cell.has_stimulation:
                cell._cre_type = meta_info.get('presynapticCre', '').lower()
                self.label_cell(cell, preEffector, positive=True) ## assume all non-patched stimulated cells are positive for preEffector


    def label_cell(self, cell, preEffector, positive=True):
        """Populate appropriate labels for a cell positive for the preEffector."""

        ## Wanted to separate this function out so it is easier to change/amend the logic later

        ## if the cell is positive for the preEffector, populate genetic labels
        if positive:
            cell.labels[preEffector] = True
            if preEffector == 'ai167':
                cell.labels['tdTomato'] = True
            elif preEffector == 'ai136':
                cell.labels['EYFP'] = True

        ## if the cell is patched (has an electrode), populate dyes 
        ##    -- this assumes which dye is used for the whole experiment based on the color of the preEffector
        if cell.electrode is not None:
            if preEffector == 'ai167':
                cell.labels['AF488'] = True
            elif preEffector == 'ai136':
                cell.labels['AF594'] = True


    def load_markPoints_connection_file(self, expt):
        with open(self.cnx_file, 'r') as f:
            exp_json = json.load(f)

        ## load stim point positions
        tseries_keys=[key for key in exp_json.keys() if 'TSeries' in key]
        points = OrderedDict()
        for tseries in tseries_keys:
            for point, data in exp_json[tseries]['MarkPoints'].items():
                pos = (data['x_pos']*1e-6, data['y_pos']*1e-6, data['z_pos']*1e-6)
                points[point] = {'pos':pos}


        ##load electrode positions and connections
        nwb_electrodes = self.get_electrodes_from_dataset()
        if nwb_electrodes is not None:
            names = ['electrode_%i'%dev_id for dev_id in nwb_electrodes]
        else:
            names = [key for key in exp_json['Headstages'].keys() if key.startswith('electrode')]

        electrodes = OrderedDict()
        cells = OrderedDict()
        HS_keys=[key for key in exp_json['Headstages'].keys() if key.startswith('electrode')]
        for headstage in names:
            data = exp_json['Headstages'].get(headstage)
            elec = Electrode(headstage, start_time=None, stop_time=None, device_id=headstage[-1])
            electrodes[headstage] = elec
            cell = Cell(expt, headstage, elec)
            elec.cell = cell
            cell.has_readout = True
            cell.has_stimulation = True ## it's possible to interogate patched cell pairs, even if that's not employed often
            cells[cell.cell_id]=cell
            if data is not None:
                cell.position = (data['x_pos']*1e-6, data['y_pos']*1e-6, data['z_pos']*1e-6) #covert from um to m
                cell.angle = data['angle']
                for p in points.keys():
                    points[p][headstage] = data['Connections'][p]



        ## check points for overlappingness (same as MP_trim in old analysis)
        distance = 10e-6
        skip=[]
        for i, p1 in enumerate(points.keys()):
            for p2 in list(points.keys())[i+1:len(points)]:
                x_dif = points[p2]['pos'][0] - points[p1]['pos'][0]
                y_dif = points[p2]['pos'][1] - points[p1]['pos'][1]
                z_dif = points[p2]['pos'][2] - points[p1]['pos'][2]
                xy_dif=math.sqrt(x_dif**2+y_dif**2)
                xyz_dif=math.sqrt(x_dif**2+y_dif**2+z_dif**2)
                if xyz_dif < distance: 
                    same = True
                    for hs in list(HS_keys):
                        if points[p1][hs] != points[p2][hs]:
                            same=False
                    if same:
                        skip.append(p1)

        ## create cells for points that were not overlapping
        for point, data in points.items():
            if point not in skip:
                cell = Cell(expt, point, None)
                cell.position = data['pos']
                cell.has_readout = False
                cell.has_stimulation = True
                cell.info.update(stim_point_ext_id=point)
                cells[cell.cell_id]=cell

        return (electrodes, cells)


    def load_mosaiceditor_connection_file(self, expt):
        with open(self.cnx_file, 'r') as f:
            exp_json = json.load(f)

        ## create Cells for stimulation points
        cells = OrderedDict()
        for name, data in exp_json['StimulationPoints'].items():
            if data['onCell']:
                cell = Cell(expt, name, None)
                cell.position = tuple(data['position'])
                cell.has_readout = False
                cell.has_stimulation = True
                cell.info.update(stim_point_ext_id=name)
                target_layer = data.get('target_layer')
                if target_layer is not None:
                    target_layer = target_layer.strip().strip('L').strip('l')
                cell._target_layer = target_layer
                cell._percent_depth = data.get('percent_depth')
                cell._distance_to_pia = data.get('distance_to_pia')
                cell._distance_to_wm = data.get('distance_to_wm')
                cells[cell.cell_id] = cell

        ## create Cells for recorded cells
        nwb_electrodes = self.get_electrodes_from_dataset()
        if nwb_electrodes is not None:
            names = ['electrode_%i'%dev_id for dev_id in nwb_electrodes]
        else:
            names = exp_json['Headstages'].keys()

        electrodes = OrderedDict()
        for name in names:
            data = exp_json['Headstages'].get(name, None)
            elec = Electrode(name, start_time=None, stop_time=None, device_id=name[-1])
            electrodes[name] = elec
            cell = Cell(expt, name, elec)
            elec.cell = cell
            cell.has_readout = True
            cell.has_stimulation = True
            cells[cell.cell_id] = cell
            if data is not None:
                cell.position = (data['x_pos'], data['y_pos'], data['z_pos'])
                cell.angle = data['angle']
                target_layer = data.get('target_layer')
                if target_layer is not None:
                    target_layer = target_layer.strip().strip('L').strip('l')
                cell._target_layer = target_layer
                cell._percent_depth = data.get('percent_depth')
                cell._distance_to_pia = data.get('distance_to_pia')
                cell._distance_to_wm = data.get('distance_to_wm')
            
        return (electrodes, cells)

    def load_without_connection_file(self, expt):
        cells = OrderedDict()
        electrodes = OrderedDict()

        stim_log_points = self.get_points_from_photostim_log()
        for name, data in stim_log_points.items():
            if data['onCell']:
                cell = Cell(expt, name, None)
                cell.position = tuple(data['position'])
                cell.has_readout = False
                cell.has_stimulation = True
                cell.info.update(stim_point_ext_id=name)
                cells[cell.cell_id] = cell

        nwb_electrodes = self.get_electrodes_from_dataset()
        if nwb_electrodes is not None:
            names = ['electrode_%i'%dev_id for dev_id in nwb_electrodes]
        else:
            raise Exception("Couldn't find .nwb file in %s. .nwb is needed to load without a connections file.")
        for name in names:
            elec = Electrode(name, start_time=None, stop_time=None, device_id=name[-1])
            electrodes[name] = elec
            cell = Cell(expt, name, elec)
            elec.cell = cell
            cell.has_readout = True
            cell.has_stimulation = True
            cells[cell.cell_id] = cell

        return (electrodes, cells)


    def create_pairs(self, cells):

        pairs = OrderedDict()
        for preCell in cells.values():
            for postCell in cells.values():
                if preCell is postCell:
                    continue
                if preCell.has_stimulation and postCell.has_readout:
                    pair = Pair(preCell, postCell)
                    pairs[(preCell.cell_id, postCell.cell_id)] = pair

        if self.cnx_file != "not found":
            with open(self.cnx_file, 'r') as f:
                exp_json = json.load(f)
            for p in pairs.values():
                try:
                    p._connection_call = exp_json['Headstages'][p.postCell.electrode.electrode_id]['Connections'][p.preCell.info.get('stim_point_ext_id')]
                    p._probed = True
                except KeyError:
                    p._probed = False
        return pairs


    def find_connections_file(self):
        cnx_files = sorted(glob.glob(os.path.join(self.site_path, '*connections*.json')))
        #print('cnx_files:', cnx_files, "path:", site_path)
        if len(cnx_files) == 1:
            return cnx_files[0]
        elif len(cnx_files) == 0:
            #raise Exception("Could not find a connections file in %s." % self.site_path)
            name = os.path.join(self._meta_info['connections_dir'], self._meta_info['experiment'])
            if os.path.exists(name):
                return name
            else:
                return "not found"
        else:
            ### return the file with the highest version number. If there's still more than one raise an exception
            max_version = 0
            cnx_file = []
            for f in cnx_files:
                version = self.get_cnx_file_version(f)
                if version > max_version:
                    max_version = version
                    cnx_file = [f]
                elif version == max_version:
                    cnx_file.append(f)

            if len(cnx_file) == 1:
                return cnx_file[0]
            else:
                raise Exception("Found %i cnx files (version:%i) in %s. Not sure which one to use: %s" %(len(cnx_file), max_version, self.site_path, [os.path.split(f)[-1] for f in cnx_file]))

    def get_cnx_file_version(self, cnx_file):
        with open(cnx_file, 'r') as f:
            exp_json = json.load(f)
        return exp_json.get('version', 0)

    def load_stimulation_log(self):
        log_files = sorted(glob.glob(os.path.join(self.site_path, 'PhotoStimulationLog_*.log')))

        log = {}
        for log_file in log_files:
            with open(log_file, 'r') as f:
                for line in f.readlines():
                    log.update(json.loads(line.rstrip(',\r\n')))

        return log

    def get_points_from_photostim_log(self):
        log = self.load_stimulation_log()
        version = log.pop('version', -1)

        points = {}

        for stim in log.values():
            if version >= 2:
                pt_name = stim['stimulationPoint']['name']
                pos = stim['stimPointPos']
                onCell=stim['stimulationPoint']['onCell']
            else:
                pt_name = 'Point ' + str(stim['stimulationPoint'])
                pos = stim['pos']
                onCell=True ### just assume this is true, because we didn't start saving this info yet

            points[pt_name] = {'name':pt_name, 'position':pos, 'onCell':onCell}

        return points






                    
