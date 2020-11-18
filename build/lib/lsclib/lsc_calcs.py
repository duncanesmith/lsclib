# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 10:33:32 2019

@author: smithd24
"""

# import common functions from default python library
import math         # import math functions
import random       # import random number functions
import numpy as np  # import numpy matrix operations
import time

# import transformation functions
import rotations as rm
from coordinate_transformations import sph2cart
from coordinate_transformations import cart2sph


def xenon_spectrum(spectrum, spectrum_max):
    """From the normalized incident light spectrum on an LSC, probabilistically
    determine an individual wavelength.
    
    Parameters
    ----------
    spectrum : DataFrame, float
        Normalized incident light spectrum or individual wavelength
    spectrum_max : float
        Maximum normalized probability
    
    Returns
    -------
    wave_len : float
        Initial wavelength of a bundle
    
    Notes
    -----
    This function should be renamed - it can handle other incident spectra, not
    just xenon.
    
    spectrum_max should be automated here, should not be an input
    """

    wave_len = 0
    if type(spectrum) is float:
        wave_len = spectrum
    else:
        wave_len = spectrum.__call__(random.uniform(0,spectrum_max))
    
    return wave_len


def pathlength_matrix(wave_len, wave_len_min, wave_len_max, absorption):
    """Determines pathlength through a volume as a function of wavelength using
    the Bouger's law.
    
    Parameters
    ----------
    wave_len : float
        Bundle wavelength
    wave_len_min : float
        Minimum wavelength absorbed by matrix
    wave_len_max : float
        Maximum wavelength absorbed by matrix
    absorption : float, DataFrame
        Float value for probability of absorption by the matrix, or the
        normalized absorption spectrum for the matrix.
    
    Returns
    -------
    matrix_path : float
        Distance a bundle can travel before it is absorbed.
    
    Notes
    -----
    """

    if type(absorption) is float:
        matrix_abco = absorption  # matrix_abco = matrix absorption coefficient
    else:
        matrix_abco = 0
        if wave_len < wave_len_min or wave_len > wave_len_max:
            matrix_abco = 10000
        if wave_len >= wave_len_min and wave_len <= wave_len_max:
            matrix_abco = absorption.__call__(wave_len)
    
    # calculate pathlength using Bouger's law
    matrix_path = ((-1 / matrix_abco) * math.log(random.uniform(0, 1)))

    return matrix_path


def surface_absorption(wave_len, wave_len_min, wave_len_max, abs_surface):
    """Determines probability of absorption as a function of wavelength at a
    particular boundary.
    
    Parameters
    ----------
    wave_len : float
        Bundle wavelength
    wave_len_min : float
        Minimum wavelength absorbed by matrix
    wave_len_max : float
        Maximum wavelength absorbed by matrix
    abs_surface : DataFrame
        Probability of absorption of a surface as a function of wavelength
    
    Returns
    -------
    probability : float
        Probability that a bundle is absorbed
    
    Notes
    -----
    surface_absorption should be renamed to boundary_absorption
    """

    probability = 0
    if wave_len >= wave_len_min and wave_len <= wave_len_max:
        probability = abs_surface.__call__(wave_len)

    return probability


def refracted_vect(n_i, n_f, vect):
    """Determines refracted angle based on incoming vector and index of
    refraction of 1st/2nd mediums
    
    Parameters
    ----------
    n_i : float
        Index of refraction of the current medium
    n_f : float
        Index of refraction of the medium a bundle is attempting to enter
    vect : NumPy array
        Current bundle trajectory in local coordinates
    
    Returns
    -------
    refr_vect : NumPy array, string
        Refracted vector - bundle trajectory after refraction. If no refraction
        occurred, the output is a string.
    
    Notes
    -----
    """

    [theta, phi] = cart2sph(vect)
    # solve Snell's law to evaluate critical angle
    reflect_test = math.sin(theta) * n_i / n_f
    # test for critical angle
    if (n_f < n_i) and (reflect_test > 1):
        refr_vect = "no refraction"
    # use Snell's law to solve for the refracted vector
    else:
        refr_angle = math.asin(reflect_test)
        # realign refr_angle within one quadrant
        if theta > (math.pi/2):  

            refr_angle = math.pi - refr_angle
        refr_vect = sph2cart(refr_angle, phi)

    return refr_vect


def reflectivity(vect, refracted_vect):
    """Calculate reflectivity at an interface using Fresnel's equations. Light
    is assumed to be unpolarized.
    
    Parameters
    ----------
    vect : NumPy array
        Current bundle trajectory in local coordinates
    refr_vect : NumPy array, string
        Refracted vector - bundle trajectory after refraction. If refraction
        was impossible, this is a string.

    Returns
    -------
    rho : float
        Reflectivity at an interface calculated using Fresnel's equations
    
    Notes
    -----
    """

    if isinstance(refracted_vect, str) is False:
        [xi, phi] = cart2sph(refracted_vect)
        [theta, phi] = cart2sph(vect)
        # apply Fresnel's equations for reflectivity to determine reflectivity
        rho = 0.5 * ((math.tan(theta - xi))**2 / (math.tan(theta + xi))**2 +
                     (math.sin(theta - xi))**2/(math.sin(theta + xi))**2)
    else:
        rho = 1  # critical angle was achieved ensuring reflection

    return rho


def refraction_or_reflection(rho, vect, tilt, refracted_vect, indexmatch,
                             n, surface_type = 'specular'):
    """Determine if a bundle is refracted or reflected at an interface.
    
    Parameters
    ----------
    rho : float
        Reflectivity at an interface calculated using Fresnel's equations
    vect : NumPy array
        Current bundle trajectory in local coordinates
    tilt : float
        Boundary angle relative to the xy plane
    refracted_vect : NumPy array, string
        Refracted vector - bundle trajectory after refraction. If refraction
        was impossible, this is a string.
    indexmatch : int
        If indexmatch is 1, bundle remains in the same volume or is lost
        If indexmatch is 0, bundle has the opportunity to move to an
        adjacent volume.
    n : NumPy array
        Normal vector of a boundary, faces inward towards center of the volume
    surface_type : string, optional
        Indicates the behavior of reflected bundles, options are:
            (1) Specular
            (2) Diffuse
            
    Returns
    -------
    theta : float
        Polar incidence angle (classically defined as just the incidence
        angle) in global coordinates. This is a spherical coordinate.
    phi : float
        Azimuthal incidence angle in global coordinates. This is a
        spherical coordinate.
    bundle_reflect_stay : int
        Bundle is reflected but remains within the LSC (1) otherwise (0)
    bundle_reflect_lost : int
        Bundle is reflected and is lost (1) otherwise (0)
    bundle_refract_lost : int
        Bundle is refracted out of the LSC (1) otherwise (0)
    Notes
    -----
    There is an excess rotation to produce "ray_vector", this could be included
    within bundle_reflection and/or bundle_reflection and processed there
    instead.
    """

    bundle_reflect_stay = 0
    bundle_reflect_lost = 0
    bundle_refract_lost = 0
    
    # ray_vector : NumPy array
    #   reflected/refracted vector in local coordinates
    # ray_normal_angle : NumPy array
    #     Angle between surface normal and reflected/refracted vector
    
    if random.uniform(0, 1) < rho:
        [theta, phi] = bundle_reflection(surface_type, vect, tilt, n)
        ray_vector = sph2cart(theta, phi)
        ray_vector = rm.rot(tilt, ray_vector, n)  # rotate into local coords
        ray_normal_angle = tilt_angle(ray_vector)
        if ray_normal_angle < 0:
            ray_normal_angle = 2*math.pi + ray_normal_angle

        # if outgoing angle will cause bundle to move in opposite direction of
        # normal, otherwise bundle stays in LSC
        if ((3*math.pi/2) > ray_normal_angle > (math.pi/2)):
            bundle_reflect_lost = 1
        else: 
            bundle_reflect_stay = 1
            
    else:
        [theta, phi] = bundle_refraction(surface_type, refracted_vect, tilt, n)
        ray_vector = sph2cart(theta, phi)
        ray_vector = rm.rot(tilt, ray_vector, n)  # rotate into local coords
        ray_normal_angle = tilt_angle(ray_vector)
        if ray_normal_angle < 0:
            ray_normal_angle = 2 * math.pi + ray_normal_angle

        # if outgoing angle will cause bundle to move in opposite direction of
        # the normal and bundle will not enter a new volume then the bundle is
        # lost. Otherwise, the bundle stays in the LSC
        if (((3 * math.pi / 2) > ray_normal_angle > (math.pi / 2))
            and (indexmatch == 1)):
            bundle_refract_lost = 1
  
    return [theta, phi, bundle_reflect_stay,
            bundle_reflect_lost, bundle_refract_lost]


def bundle_reflection(surface_type, vect, tilt, n):
    """Determine bundle trajectory upon reflection. Currently, bundles can be
    reflected specularly or diffusely.
    
    Parameters
    ----------
    vect : NumPy array
        Current bundle trajectory in local coordinates
    tilt : float
        Boundary angle relative to the xy plane
    n : NumPy array
        Normal vector of a boundary, faces inward towards center of the volume
    surface_type : string, optional
        Indicates the behavior of reflected bundles, options are:
            (1) Specular
            (2) Diffuse
            
    Returns
    -------
    theta : float
        Polar incidence angle (classically defined as just the incidence
        angle) in global coordinates. This is a spherical coordinate.
    phi : float
        Azimuthal incidence angle in global coordinates. This is a
        spherical coordinate.

    Notes
    -----
    """

    if surface_type == 'specular':
        vect = np.array(vect) * -1  # flip direction of vector
        vect = rm.z(math.pi, vect)  # rotate 180 around normal
        vect = rm.rot(-tilt, vect, n)  # rotate back to global coords
        [theta, phi] = cart2sph(vect)

        return [theta, phi]

    elif surface_type == 'diffuse':
        theta_rand = math.asin(math.sqrt(random.uniform(0, 1)))
        phi_rand = 2 * math.pi * random.uniform(0, 1)
        vect = sph2cart(theta_rand, phi_rand)
        vect = rm.rot(-tilt, vect, n)  # rotate back to global coords
        [theta, phi] = cart2sph(vect)

        return[theta, phi]

    else:
        print("The type of surface you have selected is not available")


def bundle_refraction(surface_type, refracted_vect, tilt, n):
    """Determine bundle trajectory upon refraction. Currently, bundles can be
    refracted specularly or diffusely.
    
    Parameters
    ----------
    vect : NumPy array
        Current bundle trajectory in local coordinates
    tilt : float
        Boundary angle relative to the xy plane
    n : NumPy array
        Normal vector of a boundary, faces inward towards center of the volume
    surface_type : string, optional
        Indicates the behavior of reflected bundles, options are:
            (1) Specular
            (2) Diffuse
            
    Returns
    -------
    theta : float
        Polar incidence angle (classically defined as just the incidence
        angle) in global coordinates. This is a spherical coordinate.
    phi : float
        Azimuthal incidence angle in global coordinates. This is a
        spherical coordinate.

    Notes
    -----
    """


    if surface_type == 'specular':
        vect = rm.rot(-tilt, refracted_vect, n)  # rotate back to global coords
        [theta, phi] = cart2sph(vect)

        return [theta, phi]

    elif surface_type == 'diffuse':
        theta_rand = math.asin(math.sqrt(random.uniform(0, 1)))
        phi_rand = 2 * math.pi * random.uniform(0, 1)
        vect = sph2cart(theta_rand, phi_rand)
        vect = rm.rot(-tilt, vect, n)  # rotate back to global coords
        [theta, phi] = cart2sph(vect)
        
        return[theta, phi]

    else:
        print("The type of surface you have selected is not available")


def surface_efficiency(efficiency):
    """Determine if a bundle is absorbed or not based on a float number.
    
    Parameters
    ----------
    efficiency : float
        efficiency of a surface (perfect mirror has an efficiency of 0, a
                                 perfect absorbed has an efficiency of 1)

    Returns
    -------
    reset : int
        Bundle is absorbed (1), bundle continues (2)

    Notes
    -----
    This should be combined with surface absorption.
    """

    if random.uniform(0, 1) <= efficiency:
        reset = 1
    else:
        reset = 0

    return reset


def photon_generation(wave_len, energy):
    """Determine number of photons in a bundle
    
    Parameters
    ----------
    wave_len : float
        Wavelength associated with a particular bundle
    energy : float
        Energy in a bundle

    Returns
    -------
    nphoton : float
        Number of photons in a bundle as a function of wavelength 

    Notes
    -----
    """
    
    nphoton = 0
    h = 6.626e-34  # [joule*s] planck's constant
    c = 2.998e8    # [m/s] speed of light
    ephoton = (h*c)/(wave_len*1e-9)  # find energy in a photon
    if wave_len <= 1107:
        nphoton = energy/ephoton
    
    return nphoton


def refracted_vect_flip(n, bdy_center, n_max, p_o, p_i, vect):
    """Determines if a bundle is attempting to leave an LSC or if it is
    attempting to enter the LSC. Then, finds the refracted vector.
    
    Parameters
    ----------
    n : NumPy array
        Normal vector of a boundary, faces inward towards center of the volume
    bdy_center : NumPy array
        Center point of a boundary in 3D
    n_max : float
        Index of refraction of matrix
    p_o : NumPy array
        Confirmed bundle position. This is where the bundle currently resides
        within the LSC.
    p_i : NumPy array
        Potential bundle position. This is where a bundle would intersect if
        it is not interrupted on the way (i.e. absorbed by the matrix or bundle
                                          intersects a particle).
    vect : NumPy array
        Current bundle trajectory in global coordinates
    
    Returns
    -------
    refracted_vect: NumPy array

    Notes
    -----
    This is set up with the assumption that the exterior is air, this is fine,
    but not the most elegant
    """
    
    # find point on the normal that is very far away
    exterior_point = bdy_center + np.array(100000 * n)
    # find distance between distant point and initial position
    dist_ext_point_p_o = np.linalg.norm(exterior_point - np.array(p_o))
    # find distance between distant point and potential intersect
    dist_ext_point_p_i = np.linalg.norm(exterior_point - np.array(p_i))

    # if distance to initial point is less than distance to potential
    # intersect then bundle is attempting to leave the volume, otherwise it
    # must be attempting to enter the volume.
    refract_vect = 0
    if dist_ext_point_p_o < dist_ext_point_p_i:
        refract_vect = refracted_vect(n_max, 1.0, vect)
    else:
        refract_vect = refracted_vect(1.0, n_max, vect)

    return refract_vect


def incidence_angle(n, n_o):
    """Calculate incidence angle between two vectors
    
    Parameters
    ----------
    n : NumPy array
        normal vector of second boundary
    n_o : NumPy array
        normal vector of first boundary or the xy plane (z-axis)

    Returns
    -------
    theta_i: float
        Polar incidence angle between two normals

    Notes
    -----
    This needs to be improved if more complex geometries are to be evaluated.
    Currently, boundaries must either run parallel to the x-axis or their 
    surface normal must run parallel to the x-axis.
    """    
    
    n = np.array(n)
    n_o = np.array(n_o)
    dot = np.dot(n, n_o)
    
    # calculate incidence angle in radians
    if math.isclose(dot,-1,abs_tol = 1e-7):
        theta_i = math.pi
    elif math.isclose(dot, 1, abs_tol = 1e-7):
        theta_i = 0
    else:
        theta_i = math.acos(dot)
    # flip direction of angle if cross product between two vectors is negative
    cross = np.cross(n, n_o)
    
    # if x or y coordinates are negative then flip sign of incidence angle
    if list(abs(n_o)) == [1,0,0]:
        if cross[1] < 0:
            theta_i = -theta_i
    else:
        if cross[0] < 0:
            theta_i = -theta_i
    
    return theta_i


def tilt_angle(n):
    """Calculate incidence angle between normal vector of a boundary and the
    z-axis of the global coordinate system.
    
    Parameters
    ----------
    n : NumPy array
        Normal vector of a boundary, faces inward towards center of the volume

    Returns
    -------
    theta_i: float
        Polar incidence angle between boundary normal and global z-axis

    Notes
    -----
    This should be consolidated with incidence_angle function, they are doing
    pretty much the same thing.
    
    This needs to be improved if more complex geometries are to be evaluated.
    Currently, boundaries must either run parallel to the x-axis or their 
    surface normal must run parallel to the x-axis.
    """
    
    n = np.array(n)
    n_glob = np.array([0,0,1])
    dot = np.dot(n, n_glob)  # calculate dot product of two surface normals
    
    # calculate incidence angle in radians
    if math.isclose(dot,-1,abs_tol = 1e-7):
        theta_i = math.pi
    elif math.isclose(dot, 1, abs_tol = 1e-7):
        theta_i = 0
    else:
        theta_i = math.acos(dot)
    # flip direction of angle if cross product between two vectors is negative
    cross = np.cross(n, n_glob)
    
    # if x or y coordinates are negative then flip sign of incidence angle
    if list(abs(n)) == [1,0,0]:
        if cross[1] < 0:
            theta_i = -theta_i
    else:
        if cross[0] < 0:
            theta_i = -theta_i
        
    return theta_i


def find_bdy_center(polygon):
    """Find center point of a boundary
    
    Parameters
    ----------
    polygon : NumPy array
        Set of at least three vertices in the same plane

    Returns
    -------
    center : NumPy array
        Center point of a boundary

    Notes
    -----
    """
    
    polygon = np.array(polygon)
    center = sum(polygon)/len(polygon)

    return center


def find_normal_vector(polygon, bdy_center, vol_center):
    """Given a polygon (four points), and the center of a volume, reports the
    normal vector of the corresponding plane facing inward towards the center.
    This is necessary because it will eventually help determine bundle
    intersection based on relative direction of a bundle to the direction of
    the normal vector (i.e. the way to tell if a bundle is moving upwards or
    downwards is to check trajectory against the normal)
    
    Parameters
    ----------
    polygon : list, NumPy array
        Set of at least three vertices in the same plane
    bdy_center : NumPy array
        Center point of a boundary
    vol_center : list, NumPy array
        Center point of a volume

    Returns
    -------
    unit_vector : NumPy array
        unit normal vector of a boundary facing inward
        
    Notes
    -----
    """
    
    polygon = np.array(polygon)
    vol_center = np.array(vol_center)

    # Use three points to establish two vectors in plane
    v1 = polygon[2] - polygon[0]
    v2 = polygon[1] - polygon[0]

    # the cross product of these two vectors is a vector normal to the plane
    cp = np.cross(v1, v2)
    a, b, c = cp

    # check to see if unit normal vector is facing the correct direction by
    # checking distance to tip of unit normal vector against volume center
    unit_vector = cp / np.linalg.norm(cp)
    dist_center = np.linalg.norm(bdy_center - vol_center)
    dist_with_normal = np.linalg.norm(bdy_center + unit_vector
                                      * dist_center - vol_center)
    # if distance to the center is lower than distance to tip of normal vector
    # then flip direction of unit vector to ensure it is pointing inward
    if(dist_center < dist_with_normal):
        unit_vector *= -1

    return unit_vector


def find_vol_center(bdy_points):
    """Find center point of a volume
    
    Parameters
    ----------
    bdy_points : list
        The coordinates of each vertice of a volume.

    Returns
    -------
    center : NumPy array
        Center point of a volume
        
    Notes
    -----
    """
    
    all_points = []  # initalize array to house all individual points
    # for every boundary extract points
    for j in range(0, len(bdy_points)):
        points_from_polygon = np.array(bdy_points[j])
        for k in range(0, len(points_from_polygon)):
            point_from_polygon = points_from_polygon[k]
            all_points.append(point_from_polygon)
    # eliminate all duplicate points
    unique_points = np.unique(all_points, axis=0)
    center = sum(unique_points) / len(unique_points)

    return center

def find_vol_dimensions(bdy_points):
    """Find center point of a volume
    
    Parameters
    ----------
    bdy_points : list
        The coordinates of each vertice of a volume.

    Returns
    -------
    Lx: float
        Length of the LSC in the x-direction.
    Ly: float
        Length of the LSC in the y-direction.
    Lz: float
        Length of the LSC in the z-direction.

    Notes
    -----
    """
    
    # create list of individual volume vertices
    vertices = []
    for bdy in bdy_points:
        vertices += bdy
    
    # create lists of all coordinate values
    vertexcount = len(vertices)
    x = np.zeros(vertexcount)
    y = np.zeros(vertexcount)
    z = np.zeros(vertexcount)
    for vertex in range(0, vertexcount):
        x[vertex] = vertices[vertex][0]
        y[vertex] = vertices[vertex][1]
        z[vertex] = vertices[vertex][2]
    
    # subtract highest coordinate from lowest coordinate to get lengths
    Lx = max(x) - min(x)
    Ly = max(y) - min(y)
    Lz = max(z) - min(z)
    
    return Lx, Ly, Lz
