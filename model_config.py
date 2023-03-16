from thetis import *
import thetis.coordsys as coordsys
import thetis.forcing as forcing
import csv
from netCDF4 import Dataset
import os
import scipy.interpolate as si
import numpy

#Setup timezones
sim_tz = timezone.pytz.utc
coord_system = coordsys.UTMCoordinateSystem(utm_zone=30)

"""
Function here to read station data at specific points

def readStationData():
    return

"""

def interpolate_bathymetry(bathymetry_2d, dataset="etopo1", min_depth=10.0):
    """
    Interpolate a bathymetry field from some data set
    
    :arg bathymetry_2d: :class: Firedrake 'Function' to store data in
    :kwarg dataset: the data set name, which defines the NetCDF file name
    :kwarg min_depth: minimum value to cap the bathymetry at in the shallows
    """

    if min_depth <= 0.0:
        raise NotImplementedError(
            "Minimum bathymetric depth must be postive because"
            " wetting and drying is not modelled"
        )
    
    #Make the mesh from the bathymetry function explicit
    mesh = bathymetry_2d.function_space().mesh()
    
    #Read bathymetry data from the netcdf file
    with Dataset(f"{dataset}.nc", "r") as nc:
        interp = si.RectBivariateSpline(
            nc.variables["lat"][:],
            nc.variables["lon"][:],
            nc.variables["Band1"][:, :]
        )

    #Interpolate at mesh vertices
    lonlat_func = coord_system.get_mesh_lonlat_function(mesh)
    lon, lat = lonlat_func.dat.data_ro.T
    bathymetry_2d.dat.data[:] = numpy.maximum(-interp(lat, lon), min_depth)