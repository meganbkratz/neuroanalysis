import os, glob, re
import json
from collections import OrderedDict
from .loaders import ExperimentLoader
from .mies_dataset_loader import MiesNwbLoader
from neuroanalysis.data.cell import Cell
from neuroanalysis.data.electrode import Electrode
from neuroanalysis.data.dataset import Dataset
from neuroanalysis.data.pair import Pair
from pyqtgraph import configfile


class AI_ExperimentLoader(ExperimentLoader):

    def __init__(self, site_path=None, load_file=None):
        ExperimentLoader.__init__(self)
        if (load_file is not None) and (site_path is not None):
            raise Exception('Please specify either load_file OR site_path, not both.')
        elif (load_file is None) and (site_path is None):
            raise Exception('Please specify either a load_file or a site_path.')

        self.site_path = site_path

    def find_files(self):
        files = {
            'path': self.site_path,
            'ephys': self.get_ephys_file(),
            'multipatch_log': self.get_multipatch_log(),
        }
        return files

    def get_timestamp(self):
        info = self.get_site_info()
        if info is not None:
            return info.get('__timestamp__')
            
    def get_info(self, meta_info=None):
        info = {'additional_info':meta_info}
        if self.site_path is not None:
            info['day_info']=self.get_expt_info()
            info['slice_info']=self.get_slice_info()
            info['site_info']=self.get_site_info()
        return info

    def get_surface_depth(self):
        try:
            mplog = self.get_multipatch_log()
        except TypeError:
            return None
        if mplog is None:
            return None

        lines = [l for l in open(mplog, 'r').readlines() if 'surface_depth_changed' in l]
        if len(lines) == 0:
            return None
        line = lines[-1].rstrip(',\r\n')
        return json.loads(line)['surface_depth']

    ### private:

    def get_site_info(self):
        if self.site_path is None:
            return
        index = os.path.join(self.site_path, '.index')
        if not os.path.isfile(index):
           return 
        return configfile.readConfigFile(index)['.']

    def get_slice_info(self):
        if self.site_path is None:
            return
        index = os.path.join(self.site_path, '..', '.index')
        if not os.path.isfile(index):
           return 
        return configfile.readConfigFile(index)['.']

    def get_expt_info(self):
        if self.site_path is None:
            return
        index = os.path.join(self.site_path, '../..', '.index')
        if not os.path.isfile(index):
           return 
        return configfile.readConfigFile(index)['.']

    def get_multipatch_log(self):
        files = [p for p in os.listdir(self.site_path) if re.match(r'MultiPatch_\d+.log', p)]
        if len(files) == 0:
            raise TypeError("Could not find multipatch log file for %s" % self)
        if len(files) > 1:
            raise TypeError("Found multiple multipatch log files for %s" % self)
        return os.path.join(self.site_path, files[0])

    def get_ephys_file(self):
        """Return the name of the nwb file for this experiment."""
        p = self.site_path
        files = glob.glob(os.path.join(p, '*.nwb'))
        if len(files) == 0:
            files = glob.glob(os.path.join(p, '*.NWB'))

        if len(files) == 0:
            return None
        elif len(files) > 1:
            ephys_file = None
            # multiple NWB files here; try using the file manifest to resolve.
            manifest = os.path.join(p, 'file_manifest.yml')
            if os.path.isfile(manifest):
                manifest = yaml.load(open(manifest, 'rb'))
                for f in manifest:
                    if f['category'] == 'MIES physiology':
                        ephys_file = os.path.join(os.path.dirname(self.path), f['path'])
                        break
            if ephys_file is None:
                raise Exception("Multiple NWB files found for %s" % self.expt)
        else:
            ephys_file = files[0]
        return ephys_file


class OptoExperimentLoader(AI_ExperimentLoader):
    
    def __init__(self, cnx_file=None, site_path=None):
        AI_ExperimentLoader.__init__(self, site_path=site_path, load_file=cnx_file)

        self.cnx_file = cnx_file

    def load(self, expt, meta_info=None):
        """Return a tuple of (electrodes, cells, pairs), where each element is an 
        OrderedDict of {electrode_id:Electrode}, {cell_id:Cell} or {pair_id:Pair}.
        Parameters
        ----------
        meta_info - an optional dict of meta_info to be attached to electrodes, cells or pairs.
        """

        if self.cnx_file is None:
            self.cnx_file = self.find_connections_file()

        version = self.get_cnx_file_version(self.cnx_file)

        if version == 0:
            electrodes, cells, pairs = self.load_markPoints_connection_file(expt)
        else:
            electrodes, cells, pairs = self.load_mosaiceditor_connection_file(expt)

        self.process_meta_info(electrodes, cells, meta_info)

        return (electrodes, cells, pairs)

    def get_ephys_data(self, files):
        nwb = files.get('ephys')
        if nwb is not None:
            #return OptoNwb(nwb)
            return Dataset(loader=MiesNwbLoader(nwb))

    def find_files(self):
        files = AI_ExperimentLoader.find_files(self)
        files['connections'] = self.cnx_file if self.cnx_file is not None else self.find_connections_file()
        return files

    def get_uid(self):
        if self.cnx_file is None:
            self.cnx_file = self.find_connections_file()
        return os.path.split(self.cnx_file)[1].partition('_connections')[0] ### connections file name minus '_connections.json'

#### private functions:

    def process_meta_info(self, electrodes, cells, meta_info):
        """Process optional meta_info dict that is passed in at the initialization of expt.
        """
        ## Need to load: presynapticCre, presynapticEffector, [class, reporter, layer for each headstage], internal
        if meta_info is None:
            return

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
        electrodes = OrderedDict()
        cells = OrderedDict()
        HS_keys=[key for key in exp_json['Headstages'].keys() if key.startswith('electrode')]
        for headstage in HS_keys:
            data = exp_json['Headstages'][headstage]
            elec = Electrode(headstage, start_time=None, stop_time=None, device_id=headstage[-1])
            electrodes[headstage] = elec
            cell = Cell(expt, headstage, elec)
            elec.cell = cell
            cell.position = (data['x_pos']*1e-6, data['y_pos']*1e-6, data['z_pos']*1e-6) #covert from um to m
            cell.angle = data['angle']
            cell.has_readout = True
            cell.has_stimulation = True ## it's possible to interogate patched cell pairs, even if that's not employed often
            cells[cell.cell_id]=cell
            #for p, conn in data['Connections'].items():
            #    if conn is None: ## skip this one, it doesn't have a match in points and is a duplicate
            #        continue
            #    points[p][headstage] = conn
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
                    for hs in list(electrodes.keys()):
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
                cells[cell.cell_id]=cell

        pairs = self.create_pairs(cells, exp_json)

        return (electrodes, cells, pairs)


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
                target_layer = data.get('target_layer')
                if target_layer is not None:
                    target_layer = target_layer.strip().strip('L').strip('l')
                cell._target_layer = target_layer
                cell._percent_depth = data.get('percent_depth')
                cell._distance_to_pia = data.get('distance_to_pia')
                cell._distance_to_wm = data.get('distance_to_wm')
                cells[cell.cell_id] = cell

        ## create Cells for recorded cells
        electrodes = OrderedDict()
        for name, data in exp_json['Headstages'].items():
            elec = Electrode(name, start_time=None, stop_time=None, device_id=name[-1])
            electrodes[name] = elec
            cell = Cell(expt, name, elec)
            elec.cell = cell
            cell.position = (data['x_pos'], data['y_pos'], data['z_pos'])
            cell.angle = data['angle']
            cell.has_readout = True
            cell.has_stimulation = True
            target_layer = data.get('target_layer')
            if target_layer is not None:
                target_layer = target_layer.strip().strip('L').strip('l')
            cell._target_layer = target_layer
            cell._percent_depth = data.get('percent_depth')
            cell._distance_to_pia = data.get('distance_to_pia')
            cell._distance_to_wm = data.get('distance_to_wm')
            cells[cell.cell_id] = cell

        pairs = self.create_pairs(cells, exp_json)

        return (electrodes, cells, pairs)


    def create_pairs(self, cells, exp_json):

        pairs = OrderedDict()
        for preCell in cells.values():
            for postCell in cells.values():
                if preCell is postCell:
                    continue
                if preCell.has_stimulation and postCell.has_readout:
                    pair = Pair(preCell, postCell)
                    pairs[(preCell.cell_id, postCell.cell_id)] = pair

        for p in pairs.values():
            try:
                p._connection_call = exp_json['Headstages'][p.postCell.cell_id]['Connections'][p.preCell.cell_id]
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
            raise Exception("Could not find a connections file.")
        else:
            ### return the file with the highest version number. If there's still more than one return the file with the latest modification time
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
                mts = map(os.path.getmtime, cnx_file)
                i = mts.index(max(mts))
                return cnx_file[i]

    def get_cnx_file_version(self, cnx_file):
        with open(cnx_file, 'r') as f:
            exp_json = json.load(f)
        return exp_json.get('version', 0)
                    
