import os.path, re, datetime
import numpy as np
import pyqtgraph as pg
import pyqtgraph.configfile

from aisynphys import lims, config
from aisynphys.constants import ALL_CRE_TYPES, ALL_LABELS, FLUOROPHORES, LAYERS, INJECTIONS
from aisynphys.genotypes import Genotype


class Slice():
    _cache = {}
    
    @classmethod
    def get(cls, path):
        if path not in cls._cache:
            cls._cache[path] = cls(path)
        return cls._cache[path]

    @property
    def slice_info(self):
        """Return a dictionary with meta_info about the slice."""
        raise NotImplementedError('Must be implemented in subclass.')

    @property
    def ext_id(self):
        """Return a string that identifies this slice."""
        raise NotImplementedError('Must be implemented in subclass.')

    @property
    def injections(self):
        """Return a genotype string for any viral injections or other probes added in this slice.
        """
        raise NotImplementedError('Must be implemented in subclass.')

    @property
    def genotype(self):
        """The genotype string for this specimen.
        """
        raise NotImplementedError('Must be implemented in subclass.')

    @property
    def age(self):
        """Return the age of the animal this slice came from in days at time of sectioning."""
        raise NotImplementedError('Must be implemented in subclass.')

    @property
    def sex(self):
        """Return the sex of the animal this slice came from. Options: 'M', 'F', None."""
        raise NotImplementedError('Must be implemented in subclass.')

    @property
    def species(self):
        """Return a string with the species name of the animal this slice came from. 
        Example: 'mouse' or 'human'
        """
        raise NotImplementedError('Must be implemented in subclass.')
    
    @property
    def date_of_birth(self):
        """The date of birth of the donor of this slice."""
        raise NotImplementedError('Must be implemented in subclass.')

    @property
    def orientation(self):
        """A string with the plane of section for this slice. 
        Examples: 'coronal' or 'sagital'
        """
        raise NotImplementedError('Must be implemented in subclass.')

    @property
    def surface(self):
        """A string with the surface that was exposed during the experiment. 
        Examples:'right', 'left', 'anterior', or 'posterior'
        """
        raise NotImplementedError('Must be implemented in subclass.')


    @property
    def hemisphere(self):
        """A string with the hemisphere for this slice. Options: 'left','right'
        """
        raise NotImplementedError('Must be implemented in subclass.')

    ### Not sure if we want these in the base class or not
    # @property
    # def datetime(self):
    #     ts = self.timestamp
    #     return None if ts is None else datetime.datetime.fromtimestamp(ts)

    # @property
    # def date(self):
    #     dt = self.datetime
    #     return None if dt is None else dt.date()

    @property
    def slice_time(self):
        """A datetime object representing the time of dissection of this slice, or None."""
        raise NotImplementedError('Must be implemented in subclass.')

    def __repr__(self):
        return "<%s id:%s>" % (self.__class__.__name__, self.ext_id)



class AI_Slice(Slice):
    
    def __init__(self, path):
        self.path = path
        self._slice_info = None
        self._parent_info = None
        self._genotype = None
        self._lims_record = None
        self._slice_time = None
        self._lims_specimen_name = None
        
    @property
    def slice_info(self):
        if self._slice_info is None:
            index = os.path.join(self.path, '.index')
            if not os.path.isfile(index):
                return None
            self._slice_info = pg.configfile.readConfigFile(index)['.']
        return self._slice_info
    
    @property
    def ext_id(self):
        return "%0.3f" % self.timestamp

    @property
    def lims_specimen_name(self):
        if self._lims_specimen_name is None:
            name = self.slice_info.get('specimen_ID', None)
            if name is not None:
                self._lims_specimen_name = name.strip()
            else:
                ### often opto data saves an animal_id and a slice_number -- try to reconstruct name from these
                slice_id = self.slice_info.get('slice_number', '').strip()
                if len(slice_id) > 2:
                    sid = slice_id
                else:
                    animal_id = self.parent_info.get('animal_id', '').strip()
                    if len(animal_id) == 0:
                        animal_id = self.parent_info.get('animal_ID', '').strip()
                        if len(animal_id) == 0:
                            return self._lims_specimen_name
                    if len(slice_id) == 1:
                        slice_id = '0'+slice_id
                    sid = animal_id + '.' + slice_id

                ids = lims.find_specimen_ids_matching_name(sid)
                names = []

                for i in ids:
                    name = lims.specimen_name(i)
                    m = re.match(r'(.*)(-(\d{6,7}))?(\.(\d{2}))(\.(\d{2}))$', name)
                    if m is not None:
                        names.append(name)

                if len(set(names)) == 1:
                    self._lims_specimen_name = names[0]
                elif len(set(names)) > 1:
                    raise Exception('Found more than one specimen name')
        return self._lims_specimen_name



    @property
    def parent_info(self):
        if self._parent_info is None:
            index = os.path.join(self.parent_path, '.index')
            if not os.path.isfile(index):
                raise TypeError("Cannot find index file (%s) for experiment %s" % (index, self))
            self._parent_info = pg.configfile.readConfigFile(index)['.']
        return self._parent_info

    @property
    def parent_path(self):
        return os.path.abspath(os.path.join(self.path, '..'))

    @property
    def storage_path(self):
        path = os.path.abspath(self.path)
        root = os.path.abspath(os.path.join(self.path, '..', '..'))
        return os.path.relpath(path, root)
    
    @property
    def injections(self):
        """Return a genotype string for any viral injections or other probes added in this slice.
        """
        info = self.parent_info
        inj = info.get('injections')
        if inj in (None, ''):
            return None

        # If the exact injection code is found, return the corresponding genotype description
        if inj in INJECTIONS:
            return INJECTIONS[inj]

        # injection may include multiple parts separated by " + " or " and "
        inj_parts = re.split(r'\s*\+\s*|\s+and\s+', inj)
        gtype_parts = []
        for part in inj_parts:
            if part not in INJECTIONS:
                raise KeyError("Injection part %r is unknown in constants.INJECTIONS" % part)
            gtype_parts.append(INJECTIONS[part])
        return ';'.join(gtype_parts)

    @property
    def genotype(self):
        """The genotype string for this specimen.
        """
        if self._genotype is None:
            gt_name = self.lims_record['genotype']
            inj = self.injections
            if gt_name is None:
                if inj is None:
                    return None
                else:
                    gt_parts = []
            else:
                gt_parts = gt_name.split(';')
                
            if inj is not None:
                gt_parts.extend(inj.split(';'))
                
            gt_name = ';'.join(gt_parts)
            
            self._genotype = Genotype(gt_name)
        return self._genotype

    @property
    def age(self):
        age = self.lims_record.get('age', 0)
        if self.lims_record['organism'] == 'mouse':
            if age == 0:
                raise Exception("Donor age not set in LIMS for specimen %s" % self.lims_specimen_name)
            # data not entered in to lims
            age = (self.date - self.date_of_birth).days
        return age

    @property
    def sex(self):
        return self.lims_record['sex']

    @property
    def species(self):
        return self.lims_record['organism']
    
    @property
    def date_of_birth(self):
        bd = self.lims_record['date_of_birth']
        if bd is None:
            return None
        return datetime.date(bd.year, bd.month, bd.day)

    @property
    def orientation(self):
        return self.lims_record['plane_of_section']

    @property
    def surface(self):
        return self.lims_record['exposed_surface']

    @property
    def hemisphere(self):
        return self.lims_record['hemisphere']

    @property
    def lims_record(self):
        """A dictionary of specimen information queried from LIMS.
        
        See aisynphys.lims.section_info()
        """
        if self._lims_record is None:
            if self.lims_specimen_name is not None:
                self._lims_record = lims.specimen_info(self.lims_specimen_name)
        return self._lims_record

    @property
    def timestamp(self):
        info = self.slice_info
        return None if info is None else info.get('__timestamp__', None)

    @property
    def datetime(self):
        ts = self.timestamp
        return None if ts is None else datetime.datetime.fromtimestamp(ts)

    @property
    def date(self):
        dt = self.datetime
        return None if dt is None else dt.date()

    @property
    def slice_time(self):
        if self._slice_time is None:
            slice_time = self.parent_info.get('time_of_dissection', None)
            if slice_time == '':
                slice_time = None
            if slice_time is not None:
                m = re.match(r'((20\d\d)-(\d{1,2})-(\d{1,2}) )?(\d+):(\d+)', slice_time.strip())
                if m is not None:
                    _, year, mon, day, hh, mm = m.groups()
                    if year is None:
                        date = self.date
                        slice_time = datetime.datetime(date.year, date.month, date.day, int(hh), int(mm))
                    else:
                        slice_time = datetime.datetime(int(year), int(mon), int(day), int(hh), int(mm))
            self._slice_time = slice_time
        return self._slice_time

    @property
    def quality(self):
        quality = self.slice_info.get('slice quality', None)
        try:
            quality = int(quality)
        except Exception:
            quality = None
        return quality

    @property
    def biocytin_image_url(self):
        """A LIMS URL that points to the 20x biocytin image for this specimen, or
        None if no image is found.
        """
        images = lims.specimen_images(self.lims_specimen_name)
        for image in images:
            if image['treatment'] == 'Biocytin':
                return image['url']

    @property
    def biocytin_20x_file(self):
        """File path of the 20x biocytin image for this specimen, or None if
        no image is found.
        """
        images = lims.specimen_images(self.lims_specimen_name)
        for image in images:
            if image['treatment'] == 'Biocytin':
                return image['file']

    @property
    def dapi_image_url(self):
        """A LIMS URL that points to the 20x DAPI image for this specimen, or
        None if no image is found.
        """
        images = lims.specimen_images(self.lims_specimen_name)
        for image in images:
            if image['treatment'] == 'DAPI':
                return image['url']

    @property
    def lims_drawing_tool_url(self):
        images = lims.specimen_images(self.lims_specimen_name)
        if len(images) == 0:
            return None
        else:
            return "http://lims2/drawing_tool?image_series=%d" % images[0]['image_series']
        