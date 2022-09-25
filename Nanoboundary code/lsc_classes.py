# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 15:48:05 2019

@author: smithd24
"""
# import common functions from default python library
import math         # import math functions
import random       # import random number functions
import numpy as np  # import numpy matrix operations
from inspect import isfunction
import pandas as pd

# import LSC specific functions
import lsc_calcs as lc       # import lsc calcs
import particle_calcs as pc  # import phosphor calcs

# import transformation functions
import rotations as rm                      # coord rotations
from coordinate_transformations import sph2cart  # sph to cart trans
from coordinate_transformations import cart2sph
from shapely.geometry import Polygon
from shapely.geometry import Point


class LSC():
    """
    A class used to represent an entire LSC (all volumes) and house
    the characteristics associated with a particular configuration (i.e. 
    geometry, relative position to the sun, etc.)

    ...

    Attributes
    ----------
    vol_list : list
        volumes within an LSC
    boundaries : list
        boundaries within an LSC 
    pv_boundaries : list
        boundaries that will convert photons to electricity
    gateways : list
        transparent boundaries
    interfaces: list
        transparent boundaries at the interfaces between volumes
    bundles: list
        List of all bundles run through simulation
    nphoton: int
        number of incident photons
    optical_efficiency: float
        incident bundles on pv boundaries divided by all incident bundles
    Isc_cell: float
        short circuit current of solar cells within the LSC
    Isc_i: float
        short circuit current of a bare solar cell with an equivalent surface
        area to the surface area of the LSC entrance boundary
    m: float
        spectral mismatch factor for the LSC. Normalized photon flux incident
        on pv boundaries divided by normalized photon flux incident on the
        LSC entrance boundary
    trials: int
        number of bundles
    norm_inc_spectrum: DataFrame
        normalized incident light spectrum on LSC
    norm_inc_spectrum_pv: DataFrame
        normalized incident light spectrum on pv boundaries

    Methods
    -------
    __init__()
        initializes LSC attributes
    __len__()
        dunder method to easily call number of volumes in an LSC
    __getitem__(position)
        dunder method to enable iteration over LSC volumes
    matching_pairs()
        groups together similar boundaries. Finds gateways (transparent
        boundaries), interfaces, and pv boundaries.
    main(trials, L, W, H, light_form, theta_o, phi_o, starting_vol,
         starting_bdy, I, emi_xenon, emi_xenon_max, particle)
        Processes bundle movement between volumes and calculates a variety of
        LSC attributes based upon bundle fates (simulation results)
    bare_photon_flux(bundle)
        finds number of converted photons for a bare solar cell that has an
        equivalent surface area to the surface area of the LSC entrance
        boundary.
    bdy_incident_spectrum()
        finds normalized incident light spectrum by bucketing bundles
    
    Notes
    -----
    The terminology "gateways" should be replaced by "transparent_bdys"
    
    The arguments that pass into "main" should pass directly into
    LSC_total()    
    """
  
    def __init__(self):
        """
        Parameters
        ----------
        None
        
        Returns
        -------
        None
            
        Notes
        -----
        This should be reformatted - inputs to "main" should be here.
        """
        
        self.vol_list = []
        self.boundaries = []
        self.pv_boundaries = []
        self.gateways = []
        self.interfaces = []
        self.bundles = []
        self.nphoton = 0 
        self.optical_efficiency = 0
        self.Isc_cell = 0
        self.Isc_i = 0
        self.m = 0
        
    def __len__(self):
        return len(self.vol_list)

    def __getitem__(self, position):
        return self.vol_list[position]

    def matching_pairs(self):
        """groups together similar boundaries. Finds gateways (transparent
        boundaries), interfaces, and pv boundaries.

        Parameters
        ----------
        
        self.boundaries : list
            Boundaries in the LSC
        self.gateways : list
            Transparent boundaries in the LSC
        self.pv_boundaries : list
            Boundaries made up of photovoltaic material
        self.interfaces : list
            Transparent boundaries lying between volumes
        
        Returns
        ------
        None
        
        Notes
        -----
        Class attributes should not be set within this method, they
        should be returned by the function and set in the __init__ method.
        """

        # separate gateways and pv boundaries into distinct lists
        
        # pv_true : int
        #     indicates if a boundary is a PV boundary (1) or not (0)
        # gateway_true : int
        #     indicates if a boundary is a Gateway boundary (1) or not (0)
        
        for vol in range(0, len(self.vol_list)):
            for bdy in range(0, len(self[vol])):
                self.boundaries.append(self[vol][bdy])
                pv_true = self[vol][bdy].class_type("PV")
                gateway_true = self[vol][bdy].class_type("Gateway")
                if gateway_true == 1:
                    self.gateways.append(self[vol][bdy])
                if pv_true == 1:
                    self.pv_boundaries.append(self[vol][bdy])

        # determine which gateways are interfaces (have matching pairs)
        
        # gateway_test : object
        #     Individual gateway used to determine if it has a matching
        #     counterpart (meaning it's an interface)        
        # gateway_potmatch : object
        #     Individual gateway used to compare with 'gateway_test'.
        #     In other words, this is a potential match (meaning the test is an
        #     interface)
        # nn_angle : float
        #     angle between normals of gateway_test and gateway_potmatch. If
        #     this is equal to pi radians, gateway_test and gateway_potmatch
        #     may be corresponding interfaces.
        # cc_line : NumPy array 
        #     Line connecting center points of gateway_test and
        #     gateway_potmatch.
        # ccn_angle : float
        #     Angle between cc_line and surface normal of gateway_potmatch. If
        #     this is equal to pi/2 radians, gateway_test and gateway_potmatch
        #     may be corresponding interfaces.
        
        for gateway in range(0, len(self.gateways)):
            gateway_test = self.gateways[gateway]
            self.gateways[gateway].matchingcenter = "None"

            for potmatch in range(0, len(self.gateways)):
                gateway_potmatch = self.gateways[potmatch]
                nn_angle = lc.incidence_angle(gateway_test.n,
                                              gateway_potmatch.n)
                cc_line = gateway_test.center - gateway_potmatch.center
                ccn_angle = lc.incidence_angle(cc_line,gateway_potmatch.n)
                
                # if gateways exist on the same plane with opposite normals,
                # check to see if they overlap
                
                # test_tilt : float
                #     Angle between the x-axis of the global coordinate system
                #     and the surface normal of gateway_test
                # potmatch_tilt : float
                #     Angle between the x-axis of the global coordinate system
                #     and the surface normal of gateway_potmatch.
                # rot_gateway_test : object
                #     gateway_test is rotated to lie within the xy plane of the
                #     global coordinate system. This is done to determine if 
                #     there is overlap between gateway_test and 
                #     gateway_potmatch.
                # rot_gateway_potmatch : object
                #     gateway_potmatch is rotated to lie within the xy plane
                #     of the global coordinate system. This is done to
                #     determine if there is overlap between gateway_test and
                #     gateway_potmatch.
                # intersection : float
                #     Area of overlap between rot_gateway_test and 
                #     rot_gateway_potmatch. If there is overlap,
                #     gateway_test is an interface.
                
                if math.isclose(abs(nn_angle), math.pi, abs_tol = 1e-5):
                    if math.isclose(abs(ccn_angle), math.pi/2, abs_tol = 1e-5):             
                        test_tilt = gateway_test.tilt
                        potmatch_tilt = gateway_potmatch.tilt
                        
                        #realign tilt angles within 0-3*pi/2
                        if test_tilt <= 0:
                            test_tilt += math.pi
                        if test_tilt >= 3*math.pi/2:
                            test_tilt -=math.pi
                        if potmatch_tilt <= 0:
                            potmatch_tilt += math.pi
                        if potmatch_tilt >= 3*math.pi/2:
                            potmatch_tilt -= math.pi
                        
                        rot_gateway_test = rm.rot_poly(
                                abs(test_tilt), gateway_test.polygon,
                                gateway_test.n)
                        rot_gateway_potmatch = rm.rot_poly(
                                abs(potmatch_tilt), gateway_potmatch.polygon,
                                gateway_potmatch.n)
                        intersection = rot_gateway_potmatch.intersection(
                                rot_gateway_test)
                        
                        # checking for overlap
                        if intersection.area > 1e-5:
                            self.gateways[gateway].interfaces.append(
                                gateway_potmatch)
                            self.interfaces.append(self.gateways[gateway])
                            self.gateways[gateway].matchingcenter = (
                                                self.gateways[gateway].center)

    def main(self, trials, L, W, H, light_form, theta_o, phi_o, tilt,
             starting_vol, starting_bdy, I, emi_xenon, emi_xenon_max, particle):
        """Processes bundle movement between volumes and calculates a variety
        of LSC attributes based upon bundle fates (simulation results)
        
        Parameters
        ----------
        trials: int
            Number of bundles
        L: float
            Length of the boundary where light enters
        W: float
            Width of the boundary where light enters
        H: float
            Height of the boundary where light enters
        light_form: string
            Determines distribution of light entering the LSC.
            'direct'  - light enters at a fixed angle
            'diffuse' - light enters with a Lambertian distribution
            'ground'  - light enters with a modified Lambertian distribution
                        due to the relative angle up from the ground
        theta_o: float
            initial polar incidence angle
        phi_o: float
            initial azimuthal incidence angle
        tilt: float
            angle relative to the xy plane, only relevant when modeling
            diffuse and ground-reflected light. For direct, orientation will
            be factored in later on.
        starting_vol: int
            index of the volume light will attempt to enter
        starting_bdy: int
            index of the boundary light will attempt to enter
        I: float
            Incident solar radiation (W)
        emi_xenon: DataFrame
            incident light spectrum
        emi_xenon_max: float
            maximum value of the incident light spectrum
        particle: object
        
        Returns
        ------
        None
        
        Notes
        -----
        These parameters should input to __init__ eventually
        
        "main" should call to two sub-functions:"locate_volume" and
        "lsc_outputs" locate_volume will handle bundle entrance into a volume
        and transition in between volumes. lsc_outputs will process
        post-simulation calculations
        
        L, W, and H should not be manual inputs, they should be determined
        automatically based on the selected boundary where light enters
        
        emi_xenon and emi_xenon_max should be renamed because they are not
        just limited to xenon.
        """
        
        self.trials = trials
        self.wave_len_log = []
        self.wave_len_log_pv = []
        
        for sim in range(0, trials):
            # initialize bundle location
            incx = random.uniform(0, W)  # begin at center of x        
            # for top
            incy = random.uniform(0, L)  # random location in y
            p_o = [incx, incy, H]        # create first point
            # for side
            # incz = random.uniform(0,H)  # random location in z
            # p_o = [incx, L, incz]       # create first point
            
            if (light_form == 'diffuse' or light_form == 'ground'):
                theta_o, phi_o = lc.incident_diffuse(tilt, light_form)
            
            bundle = Bundle(theta_o, phi_o, I, emi_xenon, emi_xenon_max,
                            p_o, self[starting_vol], particle)
        
            # capture characteristics of initial bundle
            self.wave_len_log.append(bundle.wave_len)
            #print(self.wave_len_log)
            self.bundles.append(bundle)            
            self.bare_photon_flux(bundle)                
            self.nphoton += lc.photon_generation(bundle.wave_len, bundle.energy)

            # determine direction bundle will be heading initially
            [bundle.theta, bundle.phi, bundle.reset, index] = (self[starting_vol][starting_bdy].find_direction(bundle))
        
            # vol_identifier: int
            #     indicates the index of the volume that a bundle is
            #     currently traveling through
            vol_identifier = starting_vol
        
            # run intersection calculations until bundle is converted or lost
            while bundle.reset == 0:
                [bundle.reset, index] = self[vol_identifier].bundle_fate(
                            vol_identifier, bundle)
                if vol_identifier != index:
                    vol_identifier = index
                if bundle.reset == 1:
                    break                
            
        
        # determine LSC spectral modifier factor, m
        
        # energy_i : float
        #     Total energy incident on entrace boundary
        # energy_cell : float
        #     Total energy incident on a pv boundary within an LSC
        # flux_bare_norm : float
        #     Number of photons converted by a pv boundary divided by the
        #     energy incident on the pv boundary
        # pv_inc_spectrum : DataFrame
        #     Log of bundle wavelengths incident on a pv boundary
        # total_photons_cell : float
        #     Number of photons converted by a pv boundary
        # total_photons_i : float
        #     Number of photons converted by a bare solar cell
        # cut_bins : NumPy array
        #     Determines bin edges for norm_inc_spectrum_pv

        self.norm_inc_spectrum = self.bdy_incident_spectrum(
                                            pd.DataFrame(self.wave_len_log),
                                            bundle.energy)
        energy_cell = np.zeros(len(self.pv_boundaries))
        flux_bare_norm = np.zeros(len(self.pv_boundaries))
        total_photons_cell = 0
        total_photons_i = np.zeros(len(self.pv_boundaries))
        m = np.zeros(len(self.pv_boundaries))        
        pv_inc_spectrum = pd.DataFrame()
        
        for pv in range(0, len(self.pv_boundaries)):
            
            # prepare inputs for individual and aggregate m calculation
            energy_cell[pv] = (
                bundle.energy*self.pv_boundaries[pv].bundles_absorbed)
            flux_cell_norm = self.pv_boundaries[pv].nphoton/energy_cell[pv]
            energy_i = bundle.energy*trials
            flux_bare_norm[pv] = (
                self.pv_boundaries[pv].nphoton_bare/energy_i)
            
            # find m for each individual pv boundary
            m[pv] = flux_cell_norm / flux_bare_norm[pv]
            self.pv_boundaries[pv].m = m[pv]
            # sum total photons for aggregate m calculation
            total_photons_cell += self.pv_boundaries[pv].nphoton
            total_photons_i[pv] = self.pv_boundaries[pv].nphoton_bare
            
            # find incident wavelengths for each individual pv boundary
            pv_inc_spectrum = pv_inc_spectrum.append(
                                        self.pv_boundaries[pv].wave_len_log)
        
        # determine total m and other LSC performance metrics
        total_flux_cell_norm = total_photons_cell/(sum(energy_cell))
        total_flux_bare_norm = (
            sum(flux_bare_norm*energy_cell))/(sum(energy_cell))
        self.m = total_flux_cell_norm/total_flux_bare_norm
        self.Isc_cell = total_photons_cell*1.60217662E-19
        self.Isc_i = ((sum(total_photons_i*energy_cell))
                         /(sum(energy_cell)))*1.60217662E-19
        self.optical_efficiency = sum(energy_cell)/energy_i
        
        # determine normalized spectrum incident on pv
        self.norm_inc_spectrum_pv = self.bdy_incident_spectrum(pv_inc_spectrum,
                                                               bundle.energy)

    
    def bare_photon_flux(self, bundle):
        
        """Finds number of converted photons for a bare solar cell that has an
        equivalent surface area to the surface area of the LSC entrance
        boundary. This is done for each pv boundary at the beginning of each
        simulation to determine the behavior of a bare solar cell.
        
        Parameters
        ----------
        bundle: object
            Houses bundle characteristics, this is used here to pull the
            initial bundle wavelength and number of photons in the bundle
        
        Returns
        -------
        None
        
        Notes
        -----
        Again, should configure this function to return attributes so that it
        is more readily apparent what these functions do in the main sim
        """
        
        for pv in range(0, len(self.pv_boundaries)):
            converted = 0
            
            if type(self.pv_boundaries[pv].EQE) is float:
                converted = lc.surface_efficiency(self.pv_boundaries[pv].EQE)
            else:
                probability = lc.surface_absorption(
                    bundle.wave_len,
                    self.pv_boundaries[pv].wave_len_min,
                    self.pv_boundaries[pv].wave_len_max,
                    self.pv_boundaries[pv].EQE)
                converted = lc.surface_efficiency(probability)
            # Uses EQE to determine if photons will become electricity
            if converted == 1:
                self.pv_boundaries[pv].nphoton_bare += (
                    lc.photon_generation(bundle.wave_len, bundle.energy))
    
    def bdy_incident_spectrum(self, wave_len_log, energy):
        """finds normalized incident light spectrum on LSC by bucketing bundles
        
        Parameters
        ----------
        None
        
        Returns
        -------
        bundle_counts: DataFrame
            Number of bundles per bin in the normalized spectrum
        
        Notes
        -----
        This should be moved into lsc_functions
        """
        cut_bins = np.linspace(0,1300,326)
        bundle_counts = (
            pd.cut(wave_len_log[0], bins=cut_bins).value_counts(sort = False))
        bundle_counts = bundle_counts/max(bundle_counts)

        return bundle_counts
    
    
class Volume(LSC):
    """
    A class used to represent the medium that a bundle will travel through. For
    an individual LSC, there can be a number of volumes. The geometry of an LSC
    will be constructed additively using individual volumes. Currently, each
    volume must be a convex polyhedron to process intersection calculations
    properly. Each volume will consist of a number of boundaries.
    
    The base volume class will not process bundle absorption and is a good
    candidate to model ideal scenarios when properties are unknown, or volumes
    consisting of air.

    ...

    Attributes
    ----------
    bdy_points : list
        The coordinates of each vertice of a volume.
    index : int
        A number differentiating volume objects.
    LSC : object
        Indicates the LSC_total object of which a volume belongs.
    IoR : float, DataFrame, expression
        Index of refraction. This can be independent of wavelength (float),
        defined with a user-defined equation that is a function of wavelength,
        or defined experimentally as a function of wavelength.
    Lx: float
        Length of the LSC in the x-direction.
    Ly: float
        Length of the LSC in the y-direction.
    Lz: float
        Length of the LSC in the z-direction.
    bdy_list: list
        The boundaries that border the volume.
    bundles_absorbed: int
        Count of bundles that have been absorbed by the matrix of this volume.
        For this class, there should be no bundles_absorbed. If there are, a
        bundle has been lost.
    boundary_intersect_count: int
        Number of times a boundary has been intersected in this volume (this
        can happen multiple times in a single simulation).
    center: NumPy array
        Mathematically derived center point of the volume.
    
    Methods
    -------
    __init__()
        initializes volume attributes.
    __len__()
        dunder method to easily call number of boundaries in a volume.
    __getitem__(position)
        dunder method to enable iteration over volume boundaries.
    bundle_fate(j, bundle)
        Tracks bundle trajectory/interactions until the bundle is either lost
        or absorbed by a collection surface.
    bdy_intersect(distance, i, j, bundle)
        Determines if a bundle is transmitted, reflected, or absorbed at a
        particular boundary
        
    Notes
    -----
    Get rid of the use of "i" and "j" and replace with more descriptive
    variables "bdy" and "vol".
    
    Code should raise a warning if a bundle does not find a boundary, instead
    of just printing that an intersect has not been found.
    """
    
    def __init__(self, bdy_points, index, LSC, IoR):
        
        self.bdy_points = bdy_points
        self.index = index
        self.LSC = LSC
        self.IoR = IoR
        [self.Lx, self.Ly, self.Lz] = lc.find_vol_dimensions(bdy_points)
        self.bdy_list = []
        self.bundles_absorbed = 0
        self.boundary_intersect_count = 0
        self.center = lc.find_vol_center(bdy_points)

    def __len__(self):
        return len(self.bdy_list)

    def __getitem__(self, position):
        return self.bdy_list[position]

    def bundle_fate(self, j, bundle):
                
        """Tracks bundle trajectory/interactions until the bundle is either
        lost or absorbed by a collection surface. A bundle will potentially
        interact with a number of boundaries within a volume, and this function
        probabilistically determines bundle position/direction.
        
        Parameters
        ----------
        
        j: int
            Index of the volume a bundle is currently traveling within
        bundle: object
            Houses bundle characteristics, this is used here to pull the
            initial bundle wavelength and number of photons in the bundle
        
        Returns
        -------
        reset: int
            Indicates if a particular simulation has come to an end (1) or if
            the bundle will continue to move (0)
        index: int
            Indicates the index of the volume that a bundle will enter next,
            it could remain in the same volume or go to an adjacent volume if
            it hits a TransparentBoundary.
            
        Notes
        -----
        "bundle.reset" is not the best parameter to return. First,
        return "reset" and then add to a bundle in the output.
    
        Replace reset with True/False eventually
        
        Replace index with "vol_index" to be more descriptive
        """

        for i in range(0, len(self.LSC[j])):
                index = self.index
                [intersect, distance] = self.LSC[j][i].find_intersect(bundle)

                # loop through boundaries until a valid intersection is found
                if isinstance(intersect, str) is False:
                    bundle.p_i = intersect
                    [bundle.reset] = self.LSC[j].bdy_intersect(
                                                        distance, i, j, bundle)
                    if bundle.reset == 1:
                        break
                    # intersect has been found, so update bundle direction
                    [bundle.theta, bundle.phi, bundle.reset, index] = (
                            self.LSC[j][i].find_direction(bundle))
                    bundle.p_o = bundle.p_i
                    break
                # bundle has looped through all boundaries with no intersection
                elif i == (len(self.LSC[j])-1):
                    bundle.reset = 1
                    self.LSC[j].bundles_absorbed += 1
                    print("bundle does not intersect, lost to volume")

        return [bundle.reset, index]

    def bdy_intersect(self, distance, i, j, bundle):
        
        """Determines if a bundle is transmitted, reflected, or absorbed at a
        particular boundary. This is evaluated just before a bundle leaves a
        boundary, because the bundle may not survive the "hit".
        
        Parameters
        ----------
        
        distance: float
            placeholder due to use in other functions called "bdy_intersect"
        i: int
            Index of the boundary a bundle has intersected
        j: int
            Index of the volume a bundle is currently traveling within
        bundle: object
            Houses bundle characteristics, this is used here to pull the
            initial bundle wavelength and number of photons in the bundle
        
        Returns
        -------
        reset: int
            Indicates if a bundle was absorbed at the boundary (1) or will
            continue to travel (0)
        
        Notes
        -----
        """
        
        self.LSC[j].boundary_intersect_count += 1
        bundle.boundary_intersect_count += 1
        reset = self.LSC[j][i].boundary_interaction(bundle)

        return [reset]


class AbsorbingVolume(Volume):
    """
    A class used to represent the medium that a bundle will travel through. For
    an individual LSC, there can be a number of volumes. The geometry of an LSC
    will be constructed additively using individual volumes. Currently, each
    volume must be a convex polyhedron to process intersection calculations
    properly. Each volume will consist of a number of boundaries.
    
    The AbsorbingVolue class will process bundle absorption and is a good
    candidate to model continuous media (i.e. silicone, glass, etc.)

    ...

    Attributes
    ----------
    wave_len_min: float
        Sets the minimum transmissible wavelength
    wave_len_max: float
        Sets the maximum transmissible wavelength
    absorption: float, DataFrame
        Probability of absorption, this can be independent of wavelength
        (float) or dependent upon wavelength (DataFrame)
    bdy_points : list
        The coordinates of each vertice of a volume.
    index : int
        A number differentiating volume objects.
    LSC : object
        Indicates the LSC_total object of which a volume belongs.
    IoR : float, DataFrame, expression
        Index of refraction. This can be independent of wavelength (float),
        defined with a user-defined equation that is a function of wavelength,
        or defined experimentally as a function of wavelength.
    Lx: float
        Length of the LSC in the x-direction.
    Ly: float
        Length of the LSC in the y-direction.
    Lz: float
        Length of the LSC in the z-direction.
    bdy_list: list
        The boundaries that border the volume.
    bundles_absorbed: int
        Count of bundles that have been absorbed by the matrix of this volume.
        For this class, there should be no bundles_absorbed. If there are, a
        bundle has been lost.
    boundary_intersect_count: int
        Number of times a boundary has been intersected in this volume (this
        can happen multiple times in a single simulation).
    center: NumPy array
        Mathematically derived center point of the volume.
    
    Methods
    -------
    bdy_intersect(distance, i, j, bundle)
        Determines if a bundle is transmitted, reflected, or absorbed at a
        particular boundary or if it was absorbed before reaching the boundary.
    __init__()
        initializes volume attributes.
    __len__()
        dunder method to easily call number of boundaries in a volume.
    __getitem__(position)
        dunder method to enable iteration over volume boundaries.
    bundle_fate(j, bundle)
        Tracks bundle trajectory/interactions until the bundle is either lost
        or absorbed by a collection surface.

    Notes
    -----  
    Get rid of the use of "i" and "j" and replace with more descriptive
    variables "bdy" and "vol".
    
    Code should raise a warning if a bundle does not find a boundary, instead
    of just printing that an intersect has not been found.
    """

    def __init__(self, bdy_points, index, LSC, IoR, absorption,
                 wave_len_min, wave_len_max):
        super().__init__(bdy_points, index, LSC, IoR)
        self.wave_len_min = wave_len_min
        self.wave_len_max = wave_len_max
        self.absorption = absorption

    def bdy_intersect(self, distance, i, j, bundle):
        """Determines first if a bundle was absorbed on the way to a boundary.
        If it survived the journey, this function evaluates if the bundle is 
        transmitted, reflected, or absorbed. This is evaluated just before a
        bundle leaves a boundary, because the bundle may not survive the "hit"
                
        Parameters
        ----------
        distance: float
            Distance that a bundle has traveled thus far through the LSC. If
            this is higher than the calculated pathlength, the bundle will be
            absorbed before it reaches a boundary.
        i: int
            Index of the boundary a bundle has intersected
        j: int
            Index of the volume a bundle is currently traveling within
        bundle: object
            Houses bundle characteristics, this is used here to pull the
            initial bundle wavelength and number of photons in the bundle
        
        Returns
        -------
        reset: int
            Indicates if a bundle was absorbed at the boundary (1) or will
            continue to travel (0)
        
        Notes
        -----
        Boundary intersect counts should be moved within the else statement
        """

        bundle.mag_trav += distance
        self.LSC[j].boundary_intersect_count += 1
        bundle.boundary_intersect_count += 1

        if bundle.mag_trav > bundle.pathlength:
            reset = 1
            self.LSC[j].bundles_absorbed += 1
        else:
            reset = self.LSC[j][i].boundary_interaction(bundle)

        return [reset]


class ParticleVolume(AbsorbingVolume):
    """
    A class used to represent the medium that a bundle will travel through. For
    an individual LSC, there can be a number of volumes. The geometry of an LSC
    will be constructed additively using individual volumes. Currently, each
    volume must be a convex polyhedron to process intersection calculations
    properly. Each volume will consist of a number of boundaries.
    
    The ParticleVolume class will process bundle absorption and is a good
    candidate to model continuous media (i.e. silicone, glass, etc.) with
    luminescent particles dispersed throughout.

    ...

    Attributes
    ----------
    particle: object
        The particle object assigned to this volume. Right now, only one
        particle may be assigned per volume.
    particle_intersect_count: int
        Logs all particle intersections, there can be many more than one
        intersection per simulation.
    particle_absorption_count: int
        Count of bundles absorbed by a particle, though not necessarily lost
    particle_emission_count: int
        Count of bundles that were absorbed and re-emitted by a particle
    wave_len_min: float
        Sets the minimum transmissible wavelength
    wave_len_max: float
        Sets the maximum transmissible wavelength
    absorption: float, DataFrame
        Probability of absorption, this can be independent of wavelength
        (float) or dependent upon wavelength (DataFrame)
    bdy_points : list
        The coordinates of each vertice of a volume.
    index : int
        A number differentiating volume objects.
    LSC : object
        Indicates the LSC_total object of which a volume belongs.
    IoR : float, DataFrame, expression
        Index of refraction. This can be independent of wavelength (float),
        defined with a user-defined equation that is a function of wavelength,
        or defined experimentally as a function of wavelength.
    Lx: float
        Length of the LSC in the x-direction.
    Ly: float
        Length of the LSC in the y-direction.
    Lz: float
        Length of the LSC in the z-direction.
    bdy_list: list
        The boundaries that border the volume.
    bundles_absorbed: int
        Count of bundles that have been absorbed by the matrix of this volume.
        For this class, there should be no bundles_absorbed. If there are, a
        bundle has been lost.
    boundary_intersect_count: int
        Number of times a boundary has been intersected in this volume (this
        can happen multiple times in a single simulation).
    center: NumPy array
        Mathematically derived center point of the volume.
    
    Methods
    -------
    
    bundle_fate(j, bundle)
        Tracks bundle trajectory/interactions until the bundle is either lost
        or absorbed by a collection surface. This variant of the function
        allows for luminescent particle interactions to be treated separately
        from the continous media they may be dispersed within.
    bdy_intersect(distance, i, j, bundle)
        Determines if a bundle is transmitted, reflected, or absorbed at a
        particular boundary or if it was absorbed before reaching the boundary.
    __init__()
        Initializes volume attributes.
    __len__()
        Dunder method to easily call number of boundaries in a volume.
    __getitem__(position)
        Dunder method to enable iteration over volume boundaries.


    Notes
    ----- 
    Get rid of the use of "i" and "j" and replace with more descriptive
    variables "bdy" and "vol".
    
    Code should raise a warning if a bundle does not find a boundary, instead
    of just printing that an intersect has not been found.
    """

    def __init__(self, bdy_points, index, LSC, IoR, absorption, particle,
                 wave_len_min, wave_len_max):
        super().__init__(bdy_points, index, LSC, IoR, absorption,
                         wave_len_min, wave_len_max)
        self.particle = particle
        self.particle_intersect_count = 0
        self.particle_absorption_count = 0
        self.particle_emission_count = 0

    def bundle_fate(self, j, bundle):
        """Tracks bundle trajectory/interactions until the bundle is either
        lost or absorbed by a collection surface. This variant of the function
        allows for luminescent particle interactions to be treated separately
        from the continous media they may be dispersed within.
        
        Parameters
        ----------
        
        j: int
            Index of the volume a bundle is currently traveling within
        bundle: object
            Houses bundle characteristics, this is used here to pull the
            initial bundle wavelength and number of photons in the bundle
        
        Returns
        -------
        reset: int
            Indicates if a particular simulation has come to an end (1) or if
            the bundle will continue to move (0)
        index: int
            Indicates the index of the volume that a bundle will enter next,
            it could remain in the same volume or go to an adjacent volume if
            it hits transparent boundary.
            
        Notes
        -----
        To be logically consistent, it seems like the bundle_fate() function
        should just be inherited instead of modified. The bdy_intersect()
        function would instead be expanded upon to incorporate the impact of a
        luminescent species.
        
        I'm not sure "bundle.reset" is a good parameter to return. First,
        return "reset" and then add to a bundle in the output.
    
        Replace reset with True/False eventually
        
        Replace index with "vol_index" to be more descriptive
        """
        
        for i in range(0, len(self)):
                index = self.index
                [intersect, distance] = self[i].find_intersect(bundle)
                
                # loop through boundaries until a valid intersection is found
                if isinstance(intersect, str) is False:
                    bundle.p_i = intersect

                    # bundle will hit a particle
                    if distance > (bundle.particle_pathlength -
                                   bundle.path_progress):
                        
                        # bundle is known to be moving towards a particle so
                        # find position based on trajectory and distance moved
                        bundle.path_progress = 0
                        bundle.mag_trav += bundle.particle_pathlength
                        bundle.p_i = bundle.p_o + sph2cart(
                                                    bundle.theta, bundle.phi,
                                                    bundle.particle_pathlength)
                        
                        # determine if a bundle is absorbed by the matrix
                        if bundle.mag_trav > bundle.pathlength:
                            bundle.reset = 1
                            self.bundles_absorbed += 1
                            break

                        # bundle will certainly hit a particle, determine if
                        # the bundle is absorbed by the particle
                        self.particle_intersect_count += 1
                        bundle.particle_intersect_count += 1
                        [bundle.reset, bundle.wave_len, bundle_emission,
                         absorbed] = self.particle.particle_interaction(
                                                               bundle.wave_len)
                        self.particle_absorption_count += absorbed
                        if bundle.reset == 1:
                            break
                        
                        # bundle has been absorbed, determine if it is emitted
                        bundle.mag_trav = bundle.mag_trav * bundle_emission
                        if bundle.mag_trav == 0:
                            self.particle_emission_count += 1
                            bundle.particle_emission_count += 1
                            bundle.pathlength = lc.pathlength_matrix(
                                            bundle.wave_len, self.wave_len_min,
                                            self.wave_len_max, self.absorption)
                        [bundle.theta, bundle.phi] = (
                                                self.particle.find_direction())
                        bundle.p_o = bundle.p_i
                        bundle.particle_pathlength = pc.pathlength(
                                                    self.particle.extinction)
                        break
                    
                    # bundle will not hit a particle before it hits a bdy
                    else:
                        [bundle.reset] = self.LSC[j].bdy_intersect(
                                                        distance, i, j, bundle)
                        if bundle.reset == 1:
                            break
                        # intersect has been found, so update bundle direction
                        [bundle.theta, bundle.phi, bundle.reset, index] = (
                                        self.LSC[j][i].find_direction(bundle))
                        bundle.path_progress += distance
                        bundle.p_o = bundle.p_i
                        break

                # bundle has looped through all boundaries with no intersection
                elif i == (len(self.LSC[j])-1):
                    bundle.reset = 1 
                    self.LSC[j].bundles_absorbed += 1

        return [bundle.reset, index]


class Boundary(Volume):
    """
    Boundaries create the border around volumes, giving them shape. This class
    houses the characteristics of a particular boundary, and will run methods
    to find bundle intersection, bundle direction, and bundle absorption

    ...

    Attributes
    ----------
    polygon : NumPy array
        Set of at least three vertices in the same plane
    volume : object
        Parent volume object to which a boundary belongs
    center : NumPy array
        Centerpoint of a boundary, this is used to find the normal 
    n : NumPy array
        Normal vector of a boundary, faces inward towards center of the volume
    tilt: float
        Angle relative to the xy plane
    polygon_for_check: Shapely polygon
        Shapely polygon object used for determining if an intersection is
        within the bounds of a boundary
    bundles_absorbed: int
        Count of bundles absorbed by a boundary
    bundles: list
        List of bundles absorbed by a boundary
    index: int
        Numberical identifier for a boundary
    nphoton: float
        Number of photons absorbed by a boundary

    Methods
    -------
    __init__(polygon, volume)
        initializes boundary attributes
    find_intersect(bundle)
        Determines if there is an intersection with a boundary. This is run for
        each of the boundaries in a volume to find where a bundle could
        intersect.
    boundary_interaction(bundle)
        Determines if a bundle will leave a boundary. For the base Boundary
        class, all bundles will be absorbed at a boundary.
    class_type(cls_type)
        Used alongside the matching_pairs function to distinguish this type
        of Boundary from other variants of the the Boundary class
        
    Notes
    -----
    """

    def __init__(self, polygon, volume):
        self.polygon = polygon
        self.volume = volume
        self.center = lc.find_bdy_center(polygon)
        self.n = lc.find_normal_vector(
                self.polygon, self.center, volume.center)
        self.tilt = lc.tilt_angle(self.n)
        self.polygon_for_check = rm.rot_poly(-self.tilt, self.polygon, self.n)
        self.bundles_absorbed = 0
        self.bundles = []
        self.index = volume.index
        self.nphoton = 0
        self.wave_len_log = []
                
    def find_intersect(self, bundle):
        """Determines if there is an intersection with a boundary. This is run
        for each of the boundaries in a volume to find where a bundle could
        intersect.
        
        Parameters
        ----------
        
        bundle: object
            Houses bundle characteristics, this is used here to pull bundle
            position and direction.
        
        Returns
        -------
        intersect: NumPy array, string
            Potential point of intersection assuming a bundle is not
            interrupted along the way. If this boundary does not have an
            intersect, the intersect is designated as 'No intersect'.
        mag_trav: float
            Distance a bundle will travel to get to an intersect on a boundary
            
        Notes
        -----
        More inputs could be placed into this function to be more descriptive.
        """
        
        # find possible intersect with line-plane intersection calculations
        
        # p_ray : NumPy array
        #     Bundle position (base of ray)
        # ray_direction : NumPy array
        #     Bundle direction (ray trajectory)
        # p_plane : NumPy array
        #     Point on boundary (could be any point on the plane)
        # n_plane : NumPy array
        #     Boundary surface normal (normal of the plane)
        # ray_plane_arr : NumPy array
        #     Line between point on plane and base of ray
        # t : float
        #     Pathlength between points found by solving the plane equation
        # solve_plane : float
        #     Solve plane equation using intersect that has been found. If this
        #     is zero, this may be a valid intersection point.
        
        p_ray = np.array(bundle.p_o)
        ray_direction = np.array(sph2cart(bundle.theta, bundle.phi))
        p_plane = np.array(self.polygon[0])
        n_plane = np.array(self.n)
        ray_plane_arr = p_ray - p_plane
        t = np.dot(n_plane, ray_plane_arr)/(-1*np.dot(n_plane, ray_direction))
        intersect = p_ray + ray_direction * t
        solve_plane = np.dot(n_plane, (intersect - p_plane))
        
        # find possible intersect with line-plane intersection calculations
        
        # mag_trav : float
        #     Distance travelled from initial bundle position to intersect
        # ray_check : float
        #     Angle of ray relative to intersected boundary
        # testpoint : Shapely Point
        #     2D Point projected onto xy plane to determine if it is within the
        #     boundary (or polygon)
        # testpolygon : Shapely Polygon
        #     2D Boundary projeted onto the xy plane to compare with point
        
        mag_trav = 0
        ray_check = 0
        if math.isclose(abs(solve_plane), 0, abs_tol=1e-15):
            testpoint = rm.rot_point(-self.tilt, intersect, self.n)
            testpolygon = self.polygon_for_check

            if testpoint.within(testpolygon):
                if abs(t) > 1e-15:
                    ray_check = lc.incidence_angle(ray_direction, self.n)
                    
                    # if angle is negative find equivalent positive angle
                    if ray_check < 0:
                        ray_check = 2*math.pi + ray_check
                    # if ray is not running parallel to boundary then intersect
                    if (3*math.pi/2) > ray_check > (math.pi/2):
                        mag_trav = abs(t)
                    else:
                        intersect = 'No intersect'
                else:
                    intersect = 'No intersect'
            else:
                intersect = 'No intersect'
        else:
            intersect = 'No intersect'

        return intersect, mag_trav

    def boundary_interaction(self, bundle):
        """Determines if a bundle will leave a boundary. For the base Boundary
        class, all bundles will be absorbed at a boundary.
        
        Parameters
        ----------
        
        bundle : object
            Houses bundle characteristics, this is used here to pull bundle
            wavelength.
        
        Returns
        -------
        reset : int
            For this base Boundary class, reset is always true (1)
            
        """
        
        reset = 1
        self.nphoton += lc.photon_generation(bundle.wave_len, bundle.energy)
        self.bundles_absorbed += 1
        self.bundles.append(bundle)
        
        return reset


    def class_type(self, cls_type):
        """Used alongside the matching_pairs function to distinguish this type
        of Boundary from other variants of the the Boundary class
        
        Parameters
        ----------
        cls_type: String
            String input as one of the three following options:
                (1) Boundary
                (2) Gateway
                (3) PV
            If the string input matches the defined logic then there is a match
        
        Returns
        -------
        match: int
            input cls_type matches string -> 1
            input cls_type does not match string -> 0
            
        Notes
        -----
        """

        if cls_type == "Boundary":
            match = 1
        else:
            match = 0
        return match

class NanoBoundary(Boundary):  
    
    # initalize same properties as boundary class, with reflectivity table, regression equations, and few added counters
    def __init__(self, polygon, volume, rho_forward, rho_backward, theta_equations, phi_equations):
        super().__init__(polygon, volume)
        
        self.bundles_reflected = 0
        self.bundles_refracted = 0
        self.matchingcenter = 0
        self.interfaces = []
        self.df_rho_forward = rho_forward
        self.df_rho_backward = rho_backward
        self.th_eq = theta_equations
        self.phi_eq = phi_equations
        
    def boundary_interaction(self, bundle):
        """Determines if a bundle will leave a boundary. For a 
        TransparentBoundary this is always true.
        
        Parameters
        ----------
        bundle: object
            Houses bundle characteristics, this is used here to pull bundle
            direction and wavelength.
        
        Returns
        -------
        reset: int
            Bundle is absorbed (1)
            Bundle is reflected (0)
        
        Notes
        -----
        """
        
        reset = 0
        
        return reset

    def find_direction(self, bundle):
        """Determine trajectory of a bundle as it leaves a boundary. For an
        OpaqueBoundary, a bundle can be reflected specularly or diffusely.
                
        Parameters
        ----------
        bundle : object
            Houses bundle characteristics, this is used here to pull bundle
            direction.
        
        Returns
        -------
        theta : float
            Polar incidence angle (classically defined as just the incidence
            angle) in global coordinates. This is a spherical coordinate.
        phi : float
            Azimuthal incidence angle in global coordinates. This is a
            spherical coordinate.
        reset : int
            Reset is set to zero to ensure that a bundle won't be lost
        index : int
            Reports the index of the volume a bundle will enter next. This
            could be the same volume if a bundle is reflected at a boundary,
            this could be a different volume if a bundle is transmitted.
        
        Notes
        -----
        It seems like reset and index are redundant and could be removed.
        
        More inputs could be placed into this function to be more descriptive.
        """
        reset = 0
        if type(self.volume.IoR) is float:
            self.IoR = self.volume.IoR
            """aqui"""
        elif type(self.volume.IoR) == type(pd.DataFrame()):
            self.IoR = self.volume.IoR.loc[round(float(bundle.wave_len))-400,0]
        else:
            self.IoR = self.volume.IoR.__call__(bundle.wave_len)
    
        # incident light vector is transformed into cartesian coordinates using
        # theta and phi from the global coordinates system  
        vect = sph2cart(bundle.theta, bundle.phi)
        
        # incident light vector is transformed into the local coordinate system
        # of a particular boundary 
        vect = rm.rot(self.tilt, vect, self.n)
        
        # Spherical coordinate of incident light vector
        [self.incident, self.azimuth] = cart2sph(vect)
        
        """
        Due to lack of simulation reflectivity data for incident angle more than 80 degrees,
        bundles with incident angle more than 80 degrees are skipped and bundles_reflected counter
        is raised to indicate bundle is skipped.
        """
        if self.incident*180/math.pi > 85:
            reset = 1
            index = 2
            self.bundles_reflected += 1
            #print("NO")
            return [bundle.theta, bundle.phi, reset, index]
        
        else:
            [self.refracted_vect, index, indexmatch] = self.check_interface(bundle, vect)
            
            # Select reflectivity value based on incident angle and wavelength
            self.backward_rho = self.df_rho_backward.loc[round(self.incident*180/math.pi), np.around(bundle.wave_len)]
            self.forward_rho = self.df_rho_forward.loc[round(self.incident*180/math.pi),np.around(bundle.wave_len)]
            #print(bundle.wave_len)
            #print(vect)
            # find if a bundle is refracted or reflected and resulting trajectory
            [theta, phi, bundle_reflect_stay, bundle_reflect_lost, bundle_refract_lost] = lc.sim_refraction_or_reflection(self.forward_rho, self.backward_rho, vect, self.tilt, self.refracted_vect, indexmatch, self.n, bundle.p_o, bundle.p_i, self.center)
        
            if bundle_reflect_stay == 1:  # bundle reflected back into volume
                index = self.index
            if bundle_reflect_lost == 1:  # bundle reflected out of LSC
                self.bundles_reflected += 1
                reset = 1
            if bundle_refract_lost == 1:  # bundle refracted out of LSC
                self.bundles_refracted += 1
                self.wave_len_log.append(bundle.wave_len)
                reset = 1
            #print(theta,phi)
            return [theta, phi, reset, index]        
        
        
        

    def check_interface(self, bundle, vect):
        """Determine if a bundle is heading into an adjacent volume or if it
        will remain in the same volume. Finds bundle trajectory upon entrance
        to a new boundary or reflection within the same boundary.
                
        Parameters
        ----------
        bundle: object
            Houses bundle characteristics, this is used here to pull bundle
            position, wavelength, and direction.
        
        Returns
        -------
        refracted_vect : NumPy array
            Direction bundle is heading in cartesian coordinates.
        index : int
            Reports the index of the volume a bundle will enter next. This
            could be the same volume if a bundle is reflected at a boundary,
            this could be a different volume if a bundle is transmitted.
        indexmatch : True/False
            If indexmatch is True, bundle remains in the same volume or is lost
            If indexmatch is False, bundle has the opportunity to move to an
            adjacent volume.
        
        Notes
        -----
        the line "elif distances[i] <= min(distances)" may be unnecessary
        """
       
        # Nanoboundary has no associated interface, it will always process sim_refracted_vect_flip
        if self.interfaces == []:
            refracted_vect = lc.sim_refracted_vect_flip(self.n, self.center, bundle.p_o, bundle.p_i, bundle.wave_len, self.incident, self.th_eq, self.phi_eq)
            index = self.index  # do not enable transfer to adjacent volume
            indexmatch = True
        # if there is an associated interface, select the first one
        else:  
            adj_interface = self.interfaces[0]
            # if there is more than one interface, determine which interface
            # a bundle should enter
            if len(self.interfaces) > 1:
                distances = np.ones(len(self.interfaces))
                testpoint = rm.rot_point(self.tilt, bundle.p_i, self.n)
                for i in range(0,len(self.interfaces)):   
                    testpolygon = rm.rot_poly(self.tilt,
                                              self.interfaces[i].polygon,
                                              self.n)
                    distances[i] = testpoint.distance(testpolygon)                    
                    # if point is within polygon, the interface has been found
                    if distances[i] == 0:
                        adj_interface = self.interfaces[i]
                        break
                    # if point is closer to an interface than the last point,
                    # update the current interface to be this one
                    elif distances[i] <= min(distances):
                        adj_interface = self.interfaces[i]     
            
            # determine index of refraction for the adjacent volume upon enter
            if type(adj_interface.volume.IoR) is float:
                adj_interface.IoR = adj_interface.volume.IoR
                '''aqui'''
            elif type(adj_interface.volume.IoR) == type(pd.DataFrame()):
                adj_interface.IoR = adj_interface.volume.IoR.loc[round(float(bundle.wave_len))-400,0]
            else:
                adj_interface.IoR = adj_interface.volume.IoR.__call__(
                                                            bundle.wave_len)
            IoR_adj = adj_interface.IoR
            refracted_vect = lc.refracted_vect(self.IoR, IoR_adj, vect)
            index = adj_interface.index
            indexmatch = False
        #print(refracted_vect)
        return refracted_vect, index, indexmatch

    def class_type(self, cls_type):
        """Used alongside the matching_pairs function to distinguish this type
        of Boundary from other variants of the the Boundary class
        
        Parameters
        ----------
        cls_type: String
            String input as one of the three following options:
                (1) Boundary
                (2) Gateway
                (3) PV
            If the string input matches the defined logic then there is a match
        
        Returns
        -------
        match: int
            input cls_type matches string -> 1
            input cls_type does not match string -> 0
            
        Notes
        -----
        This feels like a rudimentary way to do this. There might be something
        more elegant to try here.
        """
        
        if cls_type == "Gateway":
            match = 1
        else:
            match = 0
        return match    

class TransparentBoundary(Boundary):
    """
    Boundaries create the border around volumes, giving them shape. This class
    houses the characteristics of a particular boundary, and will run methods
    to find bundle intersection and bundle direction.
    
    The TransparentBoundary class is unique in that it will allow for bundle
    transmittance. The entrance boundary of an LSC must be a
    TransparentBoundary if a bundle is expected to move beyond the first
    boundary, and any interfaces between volumes must be TransparentBoundarys.

    ...

    Attributes
    ----------
    bundles_reflected : int
        Number of bundles that are reflected at the entrance boundary of the
        LSC, and will therefore never have a chance to be absorbed.
    bundles_refracted : int
        Number of bundles that are refracted out of the entrance boundary of
        the LSC, and will therefore never have a chance to be absorbed
    matchingcenter : NumPy array
        Old method of determining an interface. This will be removed.
    interfaces : list
        Interfaces to which a boundary is associated. This can be more than one
        if one boundary has multiple smaller boundaries that overlap.
    polygon : NumPy array
        Set of at least three vertices in the same plane
    volume : object
        Parent volume object to which a boundary belongs
    center : NumPy array
        Centerpoint of a boundary, this is used to find the normal 
    n : NumPy array
        Normal vector of a boundary, faces inward towards center of the volume
    tilt: float
        Angle relative to the xy plane
    polygon_for_check: Shapely polygon
        Shapely polygon object used for determining if an intersection is
        within the bounds of a boundary
    bundles_absorbed: int
        Count of bundles absorbed by a boundary
    bundles: list
        List of bundles absorbed by a boundary
    index: int
        Numberical identifier for a boundary
    nphoton: float
        Number of photons absorbed by a boundary

    Methods
    -------
    find_direction(bundle)
        Determine trajectory of a bundle as it leaves a boundary
    check_interface(bundle, vect)
        Determine if a bundle is heading into an adjacent volume or if it will
        remain in the same volume. Finds bundle trajectory upon entrance to a
        new boundary or reflection within the same boundary.
    boundary_interaction(bundle)
        Determines if a bundle will leave a boundary. For a TransparentBoundary
        this is always true.
    class_type(cls_type)
        Used alongside the matching_pairs function to distinguish this type
        of Boundary from other variants of the the Boundary class
    __init__(polygon, volume)
        initializes boundary attributes
    find_intersect(bundle)
        Determines if there is an intersection with a boundary. This is run for
        each of the boundaries in a volume to find where a bundle could
        intersect.    
    
    Notes
    ----- 
    matchingcenter attribute should be removed throughout the code.
    """
    
    # initalize same properties as boundary class, with a few added counters
    def __init__(self, polygon, volume):
        super().__init__(polygon, volume)
        
        self.bundles_reflected = 0
        self.bundles_refracted = 0
        self.matchingcenter = 0
        self.interfaces = []
        
    def boundary_interaction(self, bundle):
        """Determines if a bundle will leave a boundary. For a 
        TransparentBoundary this is always true.
        
        Parameters
        ----------
        bundle: object
            Houses bundle characteristics, this is used here to pull bundle
            direction and wavelength.
        
        Returns
        -------
        reset: int
            Bundle is absorbed (1)
            Bundle is reflected (0)
        
        Notes
        -----
        """
        
        reset = 0

        return reset

    def find_direction(self, bundle):
        """Determine trajectory of a bundle as it leaves a boundary. For an
        OpaqueBoundary, a bundle can be reflected specularly or diffusely.
                
        Parameters
        ----------
        bundle : object
            Houses bundle characteristics, this is used here to pull bundle
            direction.
        
        Returns
        -------
        theta : float
            Polar incidence angle (classically defined as just the incidence
            angle) in global coordinates. This is a spherical coordinate.
        phi : float
            Azimuthal incidence angle in global coordinates. This is a
            spherical coordinate.
        reset : int
            Reset is set to zero to ensure that a bundle won't be lost
        index : int
            Reports the index of the volume a bundle will enter next. This
            could be the same volume if a bundle is reflected at a boundary,
            this could be a different volume if a bundle is transmitted.
        
        Notes
        -----
        It seems like reset and index are redundant and could be removed.
        
        More inputs could be placed into this function to be more descriptive.
        """
        reset = 0  # initialize reset if lost through a surface of LSC
        if type(self.volume.IoR) is float:
            self.IoR = self.volume.IoR
            '''aqui'''
        elif type(self.volume.IoR) == type(pd.DataFrame()):
            self.IoR = self.volume.IoR.loc[round(float(bundle.wave_len))-400,0]
        else:
            #print(type(self.volume.IoR))
            self.IoR = self.volume.IoR.__call__(bundle.wave_len)

        # incident light vector is transformed into cartesian coordinates using
        # theta and phi from the global coordinates system  
        vect = sph2cart(bundle.theta, bundle.phi)
       # print(vect)
        # incident light vector is transformed into the local coordinate system
        # of a particular boundary 
        vect = rm.rot(self.tilt, vect, self.n)
        
        [self.refracted_vect, index, indexmatch] = self.check_interface(
                                                            bundle, vect)
        self.rho = lc.reflectivity(vect, self.refracted_vect)

        # find if a bundle is refracted or reflected and resulting trajectory
        [theta, phi, bundle_reflect_stay,
         bundle_reflect_lost, bundle_refract_lost] = (
         lc.refraction_or_reflection(self.rho, vect,
                                      self.tilt, self.refracted_vect,
                                      indexmatch, self.n))
        if bundle_reflect_stay == 1:  # bundle reflected back into volume
            index = self.index
        if bundle_reflect_lost == 1:  # bundle reflected out of LSC
            self.bundles_reflected += 1
            reset = 1
        if bundle_refract_lost == 1:  # bundle refracted out of LSC
            self.bundles_refracted += 1
            self.wave_len_log.append(bundle.wave_len)
            reset = 1
        #print(theta,phi)  
        return [theta, phi, reset, index]

    def check_interface(self, bundle, vect):
        """Determine if a bundle is heading into an adjacent volume or if it
        will remain in the same volume. Finds bundle trajectory upon entrance
        to a new boundary or reflection within the same boundary.
                
        Parameters
        ----------
        bundle: object
            Houses bundle characteristics, this is used here to pull bundle
            position, wavelength, and direction.
        
        Returns
        -------
        refracted_vect : NumPy array
            Direction bundle is heading in cartesian coordinates.
        index : int
            Reports the index of the volume a bundle will enter next. This
            could be the same volume if a bundle is reflected at a boundary,
            this could be a different volume if a bundle is transmitted.
        indexmatch : True/False
            If indexmatch is True, bundle remains in the same volume or is lost
            If indexmatch is False, bundle has the opportunity to move to an
            adjacent volume.
        
        Notes
        -----
        the line "elif distances[i] <= min(distances)" may be unnecessary
        """

        # if this boundary has no associated interface, it is leaving the LSC
        if self.interfaces == []:
            refracted_vect = lc.refracted_vect_flip(self.n, self.center,
                                                     self.IoR, bundle.p_o,
                                                     bundle.p_i, vect)
            index = self.index  # do not enable transfer to adjacent volume
            indexmatch = True
        # if there is an associated interface, select the first one
        else: 
            adj_interface = self.interfaces[0]
            # if there is more than one interface, determine which interface
            # a bundle should enter
            if len(self.interfaces) > 1:
                distances = np.ones(len(self.interfaces))
                testpoint = rm.rot_point(self.tilt, bundle.p_i, self.n)
                for i in range(0,len(self.interfaces)):   
                    testpolygon = rm.rot_poly(self.tilt,
                                              self.interfaces[i].polygon,
                                              self.n)
                    distances[i] = testpoint.distance(testpolygon)                    
                    # if point is within polygon, the interface has been found
                    if distances[i] == 0:
                        adj_interface = self.interfaces[i]
                        break
                    # if point is closer to an interface than the last point,
                    # update the current interface to be this one
                    elif distances[i] <= min(distances):
                        adj_interface = self.interfaces[i]     
            
            # determine index of refraction for the adjacent volume upon enter
            if type(adj_interface.volume.IoR) is float:
                adj_interface.IoR = adj_interface.volume.IoR
            else:
                adj_interface.IoR = adj_interface.volume.IoR.__call__(
                                                            bundle.wave_len)
            IoR_adj = adj_interface.IoR
            refracted_vect = lc.refracted_vect(self.IoR, IoR_adj, vect)
            index = adj_interface.index
            indexmatch = False
        return refracted_vect, index, indexmatch

    def class_type(self, cls_type):
        """Used alongside the matching_pairs function to distinguish this type
        of Boundary from other variants of the the Boundary class
        
        Parameters
        ----------
        cls_type: String
            String input as one of the three following options:
                (1) Boundary
                (2) Gateway
                (3) PV
            If the string input matches the defined logic then there is a match
        
        Returns
        -------
        match: int
            input cls_type matches string -> 1
            input cls_type does not match string -> 0
            
        Notes
        -----
        This feels like a rudimentary way to do this. There might be something
        more elegant to try here.
        """
        
        if cls_type == "Gateway":
            match = 1
        else:
            match = 0
        return match


class OpaqueBoundary(Boundary):
    """
    Boundaries create the border around volumes, giving them shape. This class
    houses the characteristics of a particular boundary, and will run methods
    to find bundle intersection, bundle direction, and bundle absorption

    ...

    Attributes
    ----------
    surface_type : string
        Indicates the behavior of reflected bundles, options are:
            (1) Specular
            (2) Diffuse
    efficiency : float, DataFrame
        Likelihood of bundle absorption. This can be independent of wavelength
        and given as a float number or as a function of wavelength if defined
    wave_len_min : float
        Minimum wavelength that can be absorbed
    wave_len_max : float
        Maximum wavelength bundle that can be absorbed
    polygon : NumPy array
        Set of at least three vertices in the same plane
    volume : object
        Parent volume object to which a boundary belongs
    center : NumPy array
        Centerpoint of a boundary, this is used to find the normal 
    n : NumPy array
        Normal vector of a boundary, faces inward towards center of the volume
    tilt: float
        Angle relative to the xy plane
    polygon_for_check: Shapely polygon
        Shapely polygon object used for determining if an intersection is
        within the bounds of a boundary
    bundles_absorbed: int
        Count of bundles absorbed by a boundary
    bundles: list
        List of bundles absorbed by a boundary
    index: int
        Numberical identifier for a boundary
    nphoton: float
        Number of photons absorbed by a boundary

    Methods
    -------
    __init__(polygon, volume, surface_type, efficiency,
             wave_len_min, wave_len_max)
        Initializes boundary attributes
    boundary_interaction(bundle)
        Determines if a bundle will leave a boundary. For an OpaqueBoundary,
        bundles can be absorbed or reflected.
    find_direction(bundle)
        Determine trajectory of a bundle as it leaves a boundary
    class_type(cls_type)
        Used alongside the matching_pairs function to distinguish this type
        of Boundary from other variants of the the Boundary class
    find_intersect(bundle)
        Determines if there is an intersection with a boundary. This is run for
        each of the boundaries in a volume to find where a bundle could
        intersect.
        
    Notes
    -----
    """
    
    def __init__(self, polygon, volume, surface_type,
                 efficiency, wave_len_min=0, wave_len_max=0):
        super().__init__(polygon, volume)

        self.surface_type = surface_type
        self.efficiency = efficiency
        self.wave_len_min = wave_len_min
        self.wave_len_max = wave_len_max

    def boundary_interaction(self, bundle):
        """Determine trajectory of a bundle as it leaves a boundary. For an
        OpaqueBoundary, a bundle can be reflected specularly or diffusely.
        
        Parameters
        ----------
        bundle: object
            Houses bundle characteristics, this is used here to pull bundle
            direction and wavelength.
        
        Returns
        -------
        reset: int
            Bundle is absorbed (1)
            Bundle is reflected (0)
        
        Notes
        -----
        More inputs could be placed into this function to be more descriptive.
        """
        
        reset = 0
        if type(self.efficiency) is float:
            reset = lc.surface_efficiency(self.efficiency)
        else:
            probability = lc.surface_absorption(
                    bundle.wave_len, self.wave_len_min,
                    self.wave_len_max, self.efficiency)
            reset = lc.surface_efficiency(probability)

        if reset == 1:
            self.wave_len_log.append(bundle.wave_len)
            self.nphoton += lc.photon_generation(bundle.wave_len,
                                                  bundle.energy)
            self.bundles_absorbed += 1
            self.bundles.append(bundle)
            
        return reset

    def find_direction(self, bundle):
        """Determine trajectory of a bundle as it leaves a boundary. For an
        OpaqueBoundary, a bundle can be reflected specularly or diffusely.
                
        Parameters
        ----------
        bundle: object
            Houses bundle characteristics, this is used here to pull bundle
            direction.
        
        Returns
        -------
        theta: float
            Polar incidence angle (classically defined as just the incidence
            angle) in global coordinates. This is a spherical coordinate.
        phi: float
            Azimuthal incidence angle in global coordinates. This is a
            spherical coordinate.
        reset: int
            Reset is set to zero to ensure that a bundle won't be lost
        index: int
            Only reflection is processed in this class, so the boundary index
            is refreshed to ensure the bundle stays within the same volume.
        
        Notes
        -----    
        More inputs could be placed into this function to be more descriptive.
        """
        
        reset = 0
        index = self.index
        # incident light vector is transformed into cartesian coordinates using
        # theta and phi from the global coordinates system
        vect = sph2cart(bundle.theta, bundle.phi)
        # incident light vector is transformed into the local coordinate system
        # of a particular boundary
        vect = rm.rot(self.tilt, vect, self.n)
        # incident light vector is reflected (output is in global coordinates)
        [theta, phi] = lc.bundle_reflection(self.surface_type, vect,
                                             self.tilt, self.n)
        return [theta, phi, reset, index]

    def class_type(self, cls_type):
        """Used alongside the matching_pairs function to distinguish this type
        of Boundary from other variants of the the Boundary class
        
        Parameters
        ----------
        cls_type: String
            String input as one of the three following options:
                (1) Boundary
                (2) Gateway
                (3) PV
            If the string input matches the defined logic then there is a match
        
        Returns
        -------
        match: int
            input cls_type matches string -> 1
            input cls_type does not match string -> 0
            
        Notes
        -----
        """

        if cls_type == "Boundary":
            match = 1
        else:
            match = 0
        return match
        

class PVBoundary(OpaqueBoundary):
    """
    Boundaries create the border around volumes, giving them shape. This class
    houses the characteristics of a particular boundary, and will run methods
    to find bundle intersection, bundle direction, and bundle absorption

    ...

    Attributes
    ----------
    EQE : float, DataFrame
        The external quantum efficiency (EQE) of a solar cell. This is the
        ratio of the number of charge carriers to the number of incident
        photons. This can be given as independent of wavelength (float) or as a
        function of wavelength (DataFrame).
    IQE : float, DataFrame
        The internal quantum efficiency (IQE) of a solar cell. This is the
        ratio of the number of charge carriers to the number of absorbed
        photons. This number is always higher than the EQE, and, if IQE is not
        specified, it is defaulted to be the EQE.This can be given as
        independent of wavelength (float) or as a function of wavelength
        (DataFrame).
    bundles_converted : int
        Number of bundles that are absorbed AND converted into electricity
    nphoton_bare : 
        Number of photons that would be absorbed if this boundary was placed
        over the entrance boundary of the LSC. This is used to determine the
        relative benefit of the spectral shift within the LSC.
    m : float
        spectral mismatch factor for this boundary. Normalized photon flux
        incident on a boundary divided by normalized photon flux incident on
        the LSC entrance boundary        
    surface_type : string
        Indicates the behavior of reflected bundles, options are:
            (1) Specular
            (2) Diffuse
    wave_len_min : float
        Minimum wavelength that can be absorbed
    wave_len_max : float
        Maximum wavelength bundle that can be absorbed
    polygon : NumPy array
        Set of at least three vertices in the same plane
    volume : object
        Parent volume object to which a boundary belongs
    center : NumPy array
        Centerpoint of a boundary, this is used to find the normal 
    n : NumPy array
        Normal vector of a boundary, faces inward towards center of the volume
    tilt: float
        Angle relative to the xy plane
    polygon_for_check: Shapely polygon
        Shapely polygon object used for determining if an intersection is
        within the bounds of a boundary
    bundles_absorbed: int
        Count of bundles absorbed by a boundary
    bundles: list
        List of bundles absorbed by a boundary
    index: int
        Numberical identifier for a boundary
    nphoton: float
        Number of photons absorbed by a boundary

    Methods
    -------
    __init__(polygon, volume, surface_type, efficiency,
             wave_len_min, wave_len_max)
        Initializes boundary attributes
    boundary_interaction(bundle)
        Determines if a bundle will leave a boundary. For an OpaqueBoundary,
        bundles can be absorbed or reflected.
    class_type(cls_type)
        Used alongside the matching_pairs function to distinguish this type
        of Boundary from other variants of the the Boundary class

    find_intersect(bundle)
        Determines if there is an intersection with a boundary. This is run for
        each of the boundaries in a volume to find where a bundle could
        intersect.
        
    Notes
    -----
    """
    
    def __init__(self, polygon, volume, surface_type, EQE, IQE=0,
                 wave_len_min=0, wave_len_max=0):
        super().__init__(polygon, volume, surface_type, EQE)
 
        self.surface_type = surface_type
        self.wave_len_min = wave_len_min
        self.wave_len_max = wave_len_max        
        self.EQE = EQE
        if IQE == 0:
            self.IQE = EQE
        else:
            self.IQE = IQE
        self.bundles_converted = 0
        self.nphoton_bare = 0
        self.m = 0


    def boundary_interaction(self, bundle):
        """Determine trajectory of a bundle as it leaves a boundary. For a
        PVBoundary, a bundle can be reflected specularly or diffusely.
        
        Parameters
        ----------
        bundle: object
            Houses bundle characteristics, this is used here to pull bundle
            direction and wavelength.
        
        Returns
        -------
        reset: int
            Bundle is absorbed (1)
            Bundle is reflected (0)
        
        Notes
        -----
        More inputs could be placed into this function to be more descriptive.
        """
        
        reset = 0
        converted = 0
        self.wave_len_log.append(bundle.wave_len)
        if type(self.EQE) is float:
            efficiency = self.EQE/self.IQE
            reset = lc.surface_efficiency(efficiency)
        else:
            # call EQE at a particular wavelength
            probability_EQE = lc.surface_absorption(bundle.wave_len,
                                                     self.wave_len_min,
                                                     self.wave_len_max,
                                                     self.EQE)
            # call IQE at a particular wavelength
            probability_IQE = lc.surface_absorption(bundle.wave_len,
                                                     self.wave_len_min,
                                                     self.wave_len_max,
                                                     self.IQE)
            # calculate transmittance at a particular wavelength
            if probability_EQE == 0:  # avoid divide by zero error
                probability = 0
            else:
                probability = probability_EQE/probability_IQE     
            reset = lc.surface_efficiency(probability)
        
        # if bundle is transmitted, evaluate against IQE
        if reset == 1:
            self.bundles_absorbed += 1
            self.bundles.append(bundle)
            if type(self.IQE) is float:
                converted = lc.surface_efficiency(self.IQE)
            else:
                probability = lc.surface_absorption(bundle.wave_len,
                                                     self.wave_len_min,
                                                     self.wave_len_max,
                                                     self.IQE)
                converted = lc.surface_efficiency(probability)
            if converted == 1:
                self.nphoton += lc.photon_generation(bundle.wave_len,
                                                      bundle.energy)
                self.bundles_converted += 1

        return reset

    def class_type(self, cls_type):
        """Used alongside the matching_pairs function to distinguish this type
        of Boundary from other variants of the the Boundary class
        
        Parameters
        ----------
        cls_type: String
            String input as one of the three following options:
                (1) Boundary
                (2) Gateway
                (3) PV
            If the string input matches the defined logic then there is a match
        
        Returns
        -------
        match: int
            input cls_type matches string -> 1
            input cls_type does not match string -> 0
            
        Notes
        -----
        """
    
        if cls_type == "PV":
            match = 1
        else:
            match = 0
        return match
    
    
class Bundle():
    """
    A bundle represents a "packet" of photons, and each bundle has a fixed
    amount of energy. The energy that reaches pv collection surfaces will be
    determined by the fate of each bundle in the simulation. This class logs 
    characteristics of a bundle as it travels through an LSC. "Bundle" houses
    bundle position, direction, wavelength, and distance traveled until
    absorption or particle interaction. A new bundle object is created at the
    start of every simulation.
    ...

    Attributes
    ----------
    theta : float
        Polar incidence angle (classically defined as just the incidence angle)
        in global coordinates. This is a spherical coordinate.
    phi : float
        Azimuthal incidence angle in global coordinates. This is a spherical 
        coordinate.
    I : float
        Incident solar radiation (W)
    spectrum : 
        Incident light spectrum
    spectrum_max : float
        Maximum wavelength in incident light spectrum
    p_o : NumPy array
        Confirmed bundle position. This is where the bundle currently resides
        within the LSC.
    p_i : NumPy array
        Potential bundle position. This is where a bundle would intersect if
        it is not interrupted on the way (i.e. absorbed by the matrix or bundle
                                          intersects a particle).
    volume : object
        The volume that a bundle is currently traveling within.
    particle : object
        The particle object that a bundle may intersect.
    mag_trav : float
        Distance travelled by the bundle. When mag_trav > pathlength the bundle
        is absorbed by the matrix.
    path_progress : float
        Distance travelled by the bundle toward a particle. Once a bundle
        interacts with a particle, this counter is reset.
    wave_len : float
        Current bundle wavelength
    pathlength : float
        Distance bundle can travel before it is absorbed by the matrix
    particle_pathlength : float
        Distance bundle can travel before it interacts with a particle
    reset : int
        Simulation is finished (1)
        Simulation continues (0)
    boundary_intersect_count : int
        Number of boundaries intersected by a bundle (there can be more than 1
                                                      in a single simulation)
    particle_intersect_count : int
        Number of particles intersected by a bundle (there can be more than 1
                                                     in a single simulation)
    particle_emission_count : int
        Number of particles emitted from a particle (there can be more than 1
                                                     in a single simulation)
    energy : float
        Incident solar radiation divided by number of trials. This is the
        energy of a single bundle.
    nphoton : float
        Number of photons in a bundle

    Methods
    -------
    None
    
    Notes
    -----
    """
    
    def __init__(self, theta, phi, I, spectrum, spectrum_max, p_o,
                 volume, particle, mag_trav=0, path_progress=0):

        self.theta = theta
        self.phi = phi
        self.I = I
        self.spectrum = spectrum
        self.spectrum_max = spectrum_max
        self.p_o = p_o
        self.volume = volume
        self.particle = particle
        self.mag_trav = mag_trav
        self.path_progress = path_progress
        self.p_i = p_o
        self.wave_len = lc.xenon_spectrum(self.spectrum, self.spectrum_max)
        self.pathlength = lc.pathlength_matrix(self.wave_len,
                                                self.volume.wave_len_min,
                                                self.volume.wave_len_max,
                                                self.volume.absorption)
        self.particle_pathlength = pc.pathlength(particle.extinction)
        self.reset = 0
        self.boundary_intersect_count = 0
        self.particle_intersect_count = 0
        self.particle_emission_count = 0
        self.energy = self.I/self.volume.LSC.trials
        self.nphoton = lc.photon_generation(self.wave_len, self.energy)


class Particle():
    """
    Houses characteristics of a luminescent particle while additionally
    modeling the behavior of a particle luminescent particle. This method of
    luminescent particle characterization applies to particles that are much
    larger than the wavelengths of light, and has been detailed experimentally
    using luminescent phosphors.
    ...

    Attributes
    ----------
    poa : float
        Used to determine if scattering or absorption will take place. This is
        the maximum probability of absorption, and is adjusted using the
        normalized absorption spectra of the luminescent species.
    extinction : float
        Used to determine the distance a bundle travels before reaching a
        particle. This is determined experimentally and is analogous to the
        traditional extinction coefficient of a continous material
    wave_len_min : float
        Minimum incident wavelength that can be absorbed by a particle
    wave_len_max : float
        Maximum incident wavelength that can be absorbed by a particle
    qe : float
        Quantum efficiency of a particle. Number of emitted photons divided
        by the number of absorbed photons.
    absorption_spect : DataFrame
        Absorption spectrum of the luminescent species associated with particle
    emission_spect : DataFrame
        Emission spectrum of the luminescent species associated with particle
    emission_spect_max : float
        Maximum wavelength of the emission spectrum
    bundles_absorbed: int
        Number of bundles absorbed (and not emitted) by a particle
        
    Methods
    -------
    particle_interaction(wave_len)
        Determines the fate of a bundle incident on a particle. A bundle may
        be scattered, absorbed, and/or emitted upon particle interaction.
    find_direction()
        Models isotropic emission/scattering of a bundle 

    Notes
    -----
    """
    
    def __init__(self, poa, extinction, wave_len_min, wave_len_max, qe,
                 absorption_spectrum, emission_spectrum,
                 emission_spectrum_max):
        self.poa = poa
        self.extinction = extinction
        self.wave_len_min = wave_len_min
        self.wave_len_max = wave_len_max
        self.qe = qe
        self.absorption_spect = absorption_spectrum
        self.emission_spect = emission_spectrum
        self.emission_spect_max = emission_spectrum_max
        self.bundles_absorbed = 0

    def particle_interaction(self, wave_len):
        """Determines the fate of a bundle incident on a particle. A bundle may
        be scattered, absorbed, and/or emitted upon particle interaction.
        
        Parameters
        ----------
        wave_len: float
            Incident wavelength of a bundle when it interacts with a particle
        
        Returns
        -------
        reset : int
            Bundle was absorbed and lost at a particle (1)
            Bundle continues through simulation (0)
        wave_len : float
            Modified wavelength after particle interaction. If a bundle is
            scattered there is no change in wavelength. If a bundle is emitted,
            the red-shift is applied.
        bundle_emission : int 
            Bundle has been emitted (0)
            Bundle has not been emitted (1)
        absorbed : 
            Bundle has been absorbed or emitted (1)
            Bundle has been scattered
            
        Notes
        -----
        bundle_emission is a bit confusing. Instead of this being 1 or 0, it
        should be True/False and the bundle_fate function in the scattering 
        volume class should be modified accordingly.  
        """

        reset = 0
        bundle_emission = 1
        wave_len_e = 0
        absorbed = 0
        if self.wave_len_min < wave_len < self.wave_len_max:
            # find probability that a bundle is absorbed based upon normalized
            # absorption spectrum of the phosphor
            probability = pc.absorption_spectrum(self.absorption_spect,
                                                 wave_len,
                                                 self.wave_len_min,
                                                 self.wave_len_max)
            # adjust poa with spectrum value to evaluate bundle absorption
            if random.uniform(0, 1) <= self.poa*probability:
                absorbed = 1
                wave_len_e = pc.emission_spectrum(self.emission_spect,
                                                  self.emission_spect_max)
                reset = pc.quantum_efficiency(wave_len, wave_len_e, self.qe)
                if reset == 1:
                    self.bundles_absorbed += 1
                else:
                    # bundle is emitted
                    bundle_emission = 0
                    wave_len = wave_len_e
        
        return reset, wave_len, bundle_emission, absorbed

    def find_direction(self):
        """Isotropic emission when bundle is emitted or scattered by a particle
        
        Parameters
        ----------
        None
        
        Returns
        -------
        theta : float
            Polar incidence angle (classically defined as just the incidence
            angle in global coordinates). This is a spherical coordinate. Theta
            is randomly generated to reflect isotropic emission.
        phi : float
            Azimuthal incidence angle in global coordinates. This is a
            spherical coordinate. Phi is randomly generated to reflect
            isotropic emission.
        
        Notes
        -----
        """

        theta = math.acos(1-2*random.uniform(0, 1))
        phi = 2*math.pi*random.uniform(0, 1)

        return theta, phi
