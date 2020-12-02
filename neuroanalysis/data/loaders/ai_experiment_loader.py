import os, glob, re, math, datetime
import json, yaml
from .loaders import ExperimentLoader
from neuroanalysis.data.slice import AI_Slice
from pyqtgraph import configfile
from aisynphys import lims

class AI_ExperimentLoader(ExperimentLoader):

    def __init__(self, load_file=None, site_path=None):
        ExperimentLoader.__init__(self)
        if (load_file is not None) and (site_path is not None):
            raise Exception('Please specify either load_file OR site_path, not both.')
        elif (load_file is None) and (site_path is None):
            raise Exception('Please specify either a load_file or a site_path.')

        self.site_path = site_path

    def get_datetime(self):
        """Return a datetime object corresponding to the beginning of data acquisition for this experiment."""
        ts = self.get_site_info().get('__timestamp__')
        if ts is not None:
            return datetime.datetime.fromtimestamp(ts)
        return None

    def get_target_temperature(self):
        """Return the intended temperature of the experiment in C (as a float), or None."""
        temp = self.get_expt_info().get('temperature')
        if temp is not None:
            temp = temp.lower().rstrip(' c').strip()
            if temp == 'rt':
                temp = 22.0
            elif temp == '':
                temp = None
            else:
                try:
                    temp = float(temp)
                except Exception:
                    raise ValueError('Invalid temperature: "%s"' % self.expt_info.get('temperature'))
        return temp
            
    def get_info(self):
        info={}
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

    def get_slice(self):
        if self.site_path is None:
            return None
        else:
            return AI_Slice.get(os.path.split(self.site_path)[0])

    def get_site_info(self):
        if self.site_path is None:
            return {}
        index = os.path.join(self.site_path, '.index')
        if not os.path.isfile(index):
           return  {}
        return configfile.readConfigFile(index)['.']

    def get_slice_info(self):
        if self.site_path is None:
            return {}
        index = os.path.join(self.site_path, '..', '.index')
        if not os.path.isfile(index):
           return {}
        return configfile.readConfigFile(index)['.']

    def get_expt_info(self):
        if self.site_path is None:
            return {}
        index = os.path.join(self.site_path, '..', '..', '.index')
        if not os.path.isfile(index):
           return {}
        return configfile.readConfigFile(index)['.']

    def get_multipatch_log(self):
        if self.site_path is None:
            return
        files = [p for p in os.listdir(self.site_path) if re.match(r'MultiPatch_\d+.log', p)]
        if len(files) == 0:
            raise TypeError("Could not find multipatch log file for %s" % self.site_path)
        if len(files) > 1:
            raise TypeError("Found multiple multipatch log files for %s" % self.site_path)
        return os.path.join(self.site_path, files[0])

    def get_ephys_file(self):
        """Return the name of the nwb file for this experiment."""
        if self.site_path is None:
            return
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

    def get_mosaic_file(self):
        """Return the .mosaic file for this experiment."""
        if self.site_path is None:
            return None
        sitefile = os.path.join(self.site_path, "site.mosaic")
        if not os.path.isfile(sitefile):
            sitefile = os.path.join(os.path.split(self.site_path)[0], "site.mosaic")
        if not os.path.isfile(sitefile):
            mosaicfiles = [f for f in os.listdir(self.site_path) if f.endswith('.mosaic')]
            if len(mosaicfiles) == 1:
                sitefile = os.path.join(self.site_path, mosaicfiles[0])
        if not os.path.isfile(sitefile):
            # print(os.listdir(self.path))
            # print(os.listdir(os.path.split(self.path)[0]))
            return None
        return sitefile

    def get_cluster_id(self):
        """Return the lims cell_cluster id"""
        cids = lims.expt_cluster_ids(self.get_specimen_name(), self.get_site_info()['__timestamp__'])
        if len(cids) == 0:
            return None
        if len(cids) > 1:
            raise Exception("Experiment %s has multiple LIMS clusters." % self.get_ext_id())
        return cids[0]

    def get_specimen_name(self):
        return self._expt.slice.lims_specimen_name

    def get_rig_operator(self):
        return self.get_expt_info().get('rig_operator', None)

    def get_rig_name(self):
        """Return the name of the rig used to acquire this experiment.
        """
        rig_name = self.get_expt_info().get('rig_name', None)
        # if rig_name is None:
        #     path = self.original_path.lower()
        #     m = re.search(r'\/(mp\d)', self.original_path)
        #     if m is not None:
        #        rig_name = m.groups()[0]
        return rig_name

    def get_target_region(self):
        if self._expt.slice.lims_record['organism'] == 'mouse':
            # mouse: look up in acq4 metadata
            rgn = self.get_expt_info().get('target_region', None)
            corrected = {'V1': 'VisP'}.get(rgn, rgn)
            return corrected
        else:
            # human: read from LIMS
            return self._expt.slice.lims_record['structure']

    def get_project_name(self):
        """The name of the project to which this experiment belongs.
        """
        return self.get_slice_info().get('project', None)

    @property
    def original_path(self):
        """The original path where this experiment was acquired. 
        """
        ss = os.path.join(self.site_path, 'sync_source')
        if os.path.isfile(ss):
            return open(ss, 'rb').read().decode('latin1')
        else:
            return self.site_path

    @property
    def relative_path(self):
        """The path of this experiment relative to the data repository it lives in.
        """
        repo_path = os.path.abspath(os.path.join(self.site_path, '..', '..', '..'))
        return os.path.relpath(self.site_path, repo_path)

    @property
    def server_path(self):
        """The path of this experiment relative to the server storage directory.
        """
        try:
            expt_dir = '%0.3f' % self.get_expt_info()['__timestamp__']
        except KeyError:
            raise Exception("Directory %s index is missing __timestamp__!" % self.expt_path)
        subpath = self.site_path.split(os.path.sep)[-2:]
        return os.path.join(expt_dir, *subpath)


