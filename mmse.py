
import numpy as np
import pandas as pd
import sep
import hashlib

class MultiModeSourceExtractor():
    def __init__(self, modes, bands):
        """
        Initialize MultiModeSourceExtractor instance with modes for each band
        
        Args:
            - modes (dict): dictionary with band names as keys and a list of parameters
                            for each mode as values
            - bands (dict): dictionary with band names as keys and band dictionary as values.
                            Band dictionary needs to have an image and an rms field
        """
        # check if parameters are dictionaries
        if not isinstance(modes, dict):
            raise TypeError("Parameter modes must be a dictionary")
        if not isinstance(bands, dict):
            raise TypeError("Parameter bands must be a dictionary")
        # check if name of bands in modes corresponds with those in bands
        if not set(modes.keys()).issubset(bands.keys()):
            raise ValueError("Names of bands in modes and bands dict do not match")
        # check if image and rms maps in bands have same dimension
        keys = bands.keys()
        for key in keys:
            if not bands[key]['image'].shape == bands[key]['rms'].shape:
                raise ValueError("Shape of image and rms maps in band {} do not match".format(key))
        # check that all bands have same dimensions
        image_shapes = [band['image'].shape for band in bands.values()]
        rms_shapes = [band['rms'].shape for band in bands.values()]
        if not len(set(image_shapes)) == 1:
            raise ValueError("image maps do not have same shape in all bands")
        if not len(set(rms_shapes)) == 1:
            raise ValueError("rms maps do not have same shape in all bands")
        # checks done, all good: set properties and columns
        self.modes = modes
        self.bands = bands
        self.shape = image_shapes[0]
        self.columns = ['x', 'y', 'cxx', 'cyy', 'cxy', 'a', 'b', 'theta', 'mode', 'band', 
                       'is_sub', 'main_uid', 'has_subs', 'similar_to', 'keep']

        
    def extract(self):
        """
        Run SExtractor in different modes for each band and collect detections
        for all bands in a Pandas DataFrame with unique id for each detection.
        """
        # run source extraction for every band
        data = []
        index = []
        print("Starting multi mode source extraction on {} bands".format(len(self.bands)))
        for (bid, band) in self.bands.items():
            print("- running source extraction on {} band".format(bid))
            # run source extraction for each mode in band
            for (mid, params) in self.modes[bid].items():
                objects = sep.extract(band['image'], params['sigma'], err=band['rms'], 
                                      deblend_nthresh=params['nthresh'], deblend_cont=params['mincont'],
                                      filter_kernel=params['kernel'], minarea=params['minarea'])
                for o in objects:
                    row = {}
                    index.append(self._generate_uid(o, mid, bid))
                    for c in self.columns[:8]:
                        row[c] = o[c]
                    row['mode'] = mid
                    row['band'] = bid
                    row['is_sub'] = False
                    row['main_uid'] = None
                    row['has_subs'] = False
                    row['similar_to'] = None
                    row['keep'] = None
                    data.append(row)
                print("\t + detected {} sources in {} mode".format(len(objects), mid))
        self.catalog = pd.DataFrame(data=data, columns=self.columns, index=index)
        print("Done! Detect {} objects in total across all bands.".format(self.catalog.shape[0]))
                
            
    def set_flux(self, r=1):
        """
        Commpute flux within an ellipse centered around each object
        
        Args:
            - r (float): radius of ellipse
        """
        fluxes = []
        for obj in self.catalog.itertuples(index=False):
            fluxes.append(self._get_flux(obj))
        self.catalog = self.catalog.assign(flux=fluxes)
            
        
    def find_substructures(self, main_band, main_mode, subs_band, subs_mode, percentile=0.1, r=2):
        """
        Find substructures in subs_band from subs_mode related to top percentile 
        brightest objects in main_band from main_mode
        
        Args:
            - main_band (string): name of main band
            - main_mode (string): name of main mode
            - subs_band (string): name of subs band
            - subs_mode (string): name of subs mode
            - percentile (float): top percent from which to select the objects
            - r (float): radius of ellipse
        """
        if not ('flux' in self.catalog.columns):
            raise RuntimeError("Need to set flux before searching for substrauctures")
        selection, flux_cut = self._get_flux_cut(main_band, main_mode, percentile=percentile)
        print("Looking for substructures in sources with flux count rate >= {:.2f}".format(flux_cut))
        subs = self.catalog[self.catalog['band'].str.match(subs_band) & self.catalog['mode'].str.match(subs_mode)]
        for obj in selection.itertuples():
            print("- substructures for source {}".format(obj.Index))
            ellipse = self._get_ellipse(obj, r=r)
            has_subs = False
            for sub in subs.itertuples():
                if tuple((int(sub.y), int(sub.x))) in ellipse:
                    print("\t + {}".format(sub.Index))
                    self.catalog.loc[sub.Index, 'is_sub'] = True
                    self.catalog.loc[sub.Index, 'main_uid'] = obj.Index
                    has_subs = True
            if has_subs:
                self.catalog.loc[obj.Index, 'has_subs'] = True
                
                
    def find_similar(self, band, mode1, mode2, band2=None, threshold=1e-3):
        """
        Find similar detections across two modes in the same band or in different ones
        
        Args: 
            - band (string): name of (first) band
            - mode1 (string): name of first mode
            - mode2 (string): name of second mode
            - band2 (string): name of (second) band
            - threshold (float): minimal proximity threshold
        """
        if band2==None:
            band2 = band
        objects1 = self.catalog[self.catalog['band'].str.match(band) & self.catalog['mode'].str.match(mode1)]
        objects2 = self.catalog[self.catalog['band'].str.match(band2) & self.catalog['mode'].str.match(mode2)]
        for p in objects1.itertuples():
            for q in objects2.itertuples():
                proximity = self._get_proximity(p, q)
                if proximity < threshold:
                    self.catalog.loc[q.Index, 'similar_to'] = p.Index
                    
    
    def get_positions(self, band='all', mode='all', no_duplicates=True, replace_substructures=True):
        """
        Query the catalog for the (x,y) positions of the detected objects. The parameters allow
        to choose from which band and from which mode to get the objects.
        
        Args:
            - band (string): name of band from which to get the detected objects;
                             'all' returns objects from every band
            - mode (string): name of mode from which to get the detected objects;
                             'all' returns objects from every mode
            - no_duplicates (bool): for similar objects, wheter to return both or not
            - replace_substructures (bool): wheter to replace bright objects with their
                             substructure components
                             
        Returns:
            - positions (dict or np.narray): dict with (x,y) coordinates for each band
                             or Nx2 array with (x,y) coordinates
        """
        # check if parameters are valid
        bids = self.bands.keys()
        mids = []
        for bid in bids:
            mids.extend(list(self.modes[bid].keys()))
        mids = set(mids)
        if (band!='all') and not(band in bids):
            raise ValueError("The band name specified does not exist")
        if (band!='all') and not(mode in mids):
            raise ValueError("The mode name specified does not exist")
        if band=='all' and mode=='all':
            breg = '|'.join(bids)
            mreg = '|'.join(mids)
        else:
            breg = band
            mreg = mode
        # select objects in bands and modes according to specified parameters
        print("Selecting objects from: \n\t + bands: {} \n\t + modes: {}".format(breg, mreg))
        objects = self.catalog[self.catalog['band'].str.match(breg) & self.catalog['mode'].str.match(mreg) & self.catalog['keep']]
        objects = objects.loc[:, ['x', 'y', 'band']]
        """
        # remove duplicates if requested
        if no_duplicates:
            dups = objects[objects['similar_to'].notnull()]
            objects = objects.drop(dups.index)
            print("- removing duplicates: \n\t + {} objects dropped".format(dups.shape[0]))
        # replace substructures if requested
        if replace_substructures:
            mains = objects[objects['has_subs']]
            objects = objects.drop(mains.index)
            print("- replacing substructure: \n\t + {} objects dropped".format(mains.shape[0]))
        """
        # prepare output, if all bands were considered in the previous steps, then return
        # a dictionary with the band ids as keys and the position matrix as values;
        # else simply return a Nx2 matrix with (x,y) coordinates
        if band=='all':
            positions = {}
            for bid in bids:
                positions[bid] = objects[objects['band'].str.match(bid)].as_matrix(columns=['x', 'y'])
        else:
            positions = objects.as_matrix(columns=['x', 'y'])
        return positions
    
    
    def keep(self, uids):
        """
        Set keep flag to True for objects with uid in uids
        
        Args:
            - uids (list): list of uid indices
        """
        self.catalog.loc[uids, 'keep'] = True
        
    
    def discard(self, uids):
        """
        Set keep flag to False for objects with uid in uids
        
        Args:
            - uids (list): list of uid indices
        """
        self.catalog.loc[uids, 'keep'] = False
            
            
    def _generate_uid(self, obj, mode, band):
        """
        Generate a unique id for the detected object, using its parameters as seed
        for the hash function
        
        Args:
            - obj (np.narray): detected object as structured array
            - mode (string): name of mode in which the object was detected
            - band (string): name of band in which the object was detected
        
        Returns:
            - uid (string): md5 hash
        """
        seed = ""
        for c in self.columns[:8]:
            seed += str(obj[c])
        seed += mode+band
        return hashlib.md5(seed.encode('utf-8')).hexdigest()
    
    
    def _get_ellipse(self, obj, r=1):
        """
        Find the set of pixel that satisfy the ellipse inequality:
         CXX*(x-xb)**2+CYY*(y-yb)**2+CXY*(x-xb)*(y-yb)<=r**2
         
        Args:
            - obj (np.narray): detected object as structured array
            - r (float): radius of ellipse
            
        Returns:
            - ellipse (set): set of pixel satisfying the above inequality
        """
        ellipse = []
        w1, w2 = 50, 50
        h1, h2 = 50, 50
        if (obj.x < w1):
            w1 = int(obj.x)
        if (obj.x >= self.shape[1]-w2):
            w2 = self.shape[1]-int(obj.x)
        if (obj.y < h1):
            w1 = int(obj.y)
        if (obj.y >= self.shape[0]-h2):
            w2 = self.shape[0]-int(obj.y)
        for y in range(int(obj.y)-h1, int(obj.y)+h2-1):
            for x in range(int(obj.x)-w1, int(obj.x)+w2-1):
                if (obj.cxx*(x-obj.x)**2+obj.cyy*(y-obj.y)**2+obj.cxy*(x-obj.x)*(y-obj.y)<=r**2):
                    ellipse.append((y,x))
        return set(ellipse)
    
    
    def _get_flux(self, obj, r=1):
        """
        Compute the flux of the object within the ellipse defined by:
         CXX*(x-xb)**2+CYY*(y-yb)**2+CXY*(x-xb)*(y-yb)<=r**2
         
        Args:
            - obj (np.narray): detected object as structured array
            - r (float): radius of ellipse
            
        Returns:
            - flux (float): total flux within the ellipse centered around the object
        """
        image = self.bands[obj.band]['image']
        flux = []
        w1, w2 = 50, 50
        h1, h2 = 50, 50
        if (obj.x < w1):
            w1 = int(obj.x)
        if (obj.x >= self.shape[1]-w2):
            w2 = self.shape[1]-int(obj.x)
        if (obj.y < h1):
            w1 = int(obj.y)
        if (obj.y >= self.shape[0]-h2):
            w2 = self.shape[0]-int(obj.y)
        for y in range(int(obj.y)-h1, int(obj.y)+h2-1):
            for x in range(int(obj.x)-w1, int(obj.x)+w2-1):
                if (obj.cxx*(x-obj.x)**2+obj.cyy*(y-obj.y)**2+obj.cxy*(x-obj.x)*(y-obj.y)<=r**2):
                    flux.append(image[y,x])
        flux = np.sum(flux)
        return flux
    
    
    def _get_flux_cut(self, band, mode, percentile=0.1):
        """
        Select top percentile brightest objects in band from a given mode
        
        Args:
            - band (string): name  of band
            - mode (string): name of mode
            - percentile (float): top percent from which to select the objects
            
        Returns:
            - selection (pd.DataFrame): view in to the catalog, with objects in top flux percentile
            - cut (float): minimal flux count rate for top flux percentile
        """
        selection = self.catalog[self.catalog['band'].str.match(band) & self.catalog['mode'].str.match(mode)]
        cut = selection['flux'].quantile(q=1.0-percentile, interpolation='lower')
        return selection[selection['flux']>=cut], cut
    
    def _get_proximity(self, obj1, obj2):
        """
        Compute proximity of two objects as difference between their respective barycentres,
        elliptical areas and orientations
        
        Args:
            - obj1 (np.narray): first detected object as structured array
            - obj2 (np.narray): second detected object as structured array
            
        Returns:
            - proximity (float): distance between objects in the 3D space defined by
                                 barycentre position, elliptical area and orientation
        """
        d = np.sqrt((obj1.x-obj2.x)**2 + (obj1.y-obj2.y)**2)
        a = np.abs((obj1.a*obj1.b) - (obj2.a*obj2.b))
        t = np.abs(obj1.theta - obj2.theta)
        return d + a + t