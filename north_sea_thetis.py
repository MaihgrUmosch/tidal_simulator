#Casts bathymetry onto mesh.

import os
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import datetime
import gmsh
import time as time_mod
os.environ["OMP_NUM_THREADS"] = "1"

from thetis import *
from thetis.forcing import *

sim_tz = timezone.pytz.utc
coord_system = coordsys.UTMCoordinateSystem(utm_zone=30)

#Bathymetry
def readBathymetry(filename):
    """
    Returns interpolator for specified bathymetry
    """
    b = Dataset(filename)
    lat = b.variables['lat'][:] #y
    lon = b.variables['lon'][:] #x
    elv = b.variables['elevation'][:]
    bath = -elv
    bath = bath.filled(5) #Replaced masked data with 5
    
    return lat, lon, bath


########################################################################


lat, lon, bath = readBathymetry('/home/epipremnum/Documents/tidal_dam_model/north_sea_example/bathymetric_data/gebco_2022_n64.0723_s47.0215_w-13.9746_e12.7441.nc')
interpolator = RegularGridInterpolator((lon, lat), bath.T)
"""
gmsh.initialize()
gmsh.open('dungeness.msh')
ocean_nodes = gmsh.model.mesh.get_nodes_for_physical_group(1,2)[0]

gmsh.finalize()
"""

mesh = Mesh('/home/epipremnum/Documents/tidal_dam_model/dungeness_contour/dungeness_wall_A_v2.msh', dim = 2)
H = get_functionspace(mesh, 'Lagrange', 1)
V = get_functionspace(mesh, 'Lagrange', 1, vector=True)
elev_field  = Function(H, name='elevation')
uv_field    = Function(V, name='transport')
bathymetry  = Function(H, name = 'bathymetry')

#Check mesh sits within bounds of bathymetry file
mesh_coords = mesh.coordinates
mesh_lon = []
mesh_lat = []
for coord in mesh_coords.dat.data[:]:
    mesh_lon.append(coord[0])
    mesh_lat.append(coord[1])

b_max_lat = max(lat)
b_max_lon = max(lon)
b_min_lat = min(lat)
b_min_lon = min(lon)
m_max_lat = max(mesh_lat)
m_max_lon = max(mesh_lon)
m_min_lat = min(mesh_lat)
m_min_lon = min(mesh_lon)

if max(lat) < max(mesh_lat):
    print("ERROR 1A: max bathymetry latitude is smaller than the max mesh latitude. Use taller bathymetry file.")
if max(lon) < max(mesh_lon):
    print("ERROR 1B: max bathymetry longitude is smaller than the max mesh longitude. Use wider bathymetry file.")
if min(lat) > min(mesh_lat):
    print("ERROR 1C: min bathymetry latitude is greater than the min mesh latitude. Use taller bathymetry file.")
if min(lon) > min(mesh_lon):
    print("ERROR 1D: min bathymetry longitude is greater than the min mesh longitude. Use wider bathymetry file.")

#Interpolate bathymetry onto mesh
bathymetry.dat.data[:] = numpy.maximum(interpolator(mesh_coords.dat.data[:]),1) 
lat, lon = coord_system.get_mesh_lonlat_function(mesh)


########################################################################


def constructSolver(spinup=False, store_station_time_series=True, **model_options):
    """
    Construct a :class: Firedrake 'FlowSolver2d' instance for inverse modelling

    :kwarg spinup: Boolean to start a spin-up run or a simulation
    :kwarg store_station_time_series: Boolean to indicate whether gauge measurements should be stored
    :return: :class: 'FlowSolver2d' instance, the start date for the simulation and a function for updating forcings
    """

    #Setup Manning friction
    manning_2d = Function(H, name="Manning coefficient")
    manning_2d.assign(3.0e-2)

    #Setup Coriolis forcing
    omega = 7.292e-5
    coriolis_2d = Function(H, name = "Coriolis forcing")
    coriolis_2d.interpolate(2*omega*sin(lat * pi / 180.0))

    #Setup temporal discretisation
    default_start_date = datetime.datetime(2022, 1, 1, tzinfo=sim_tz)
    default_end_date = datetime.datetime(2022, 1, 2, tzinfo=sim_tz)
    start_date = model_options.pop("start_date", default_start_date)
    end_date = model_options.pop("end_date", default_end_date)
    dt = 3600.0         # one hour
    t_export = 3600.0   # export time of every hour
    t_end = (end_date - start_date).total_seconds()

    if os.getenv("THETIS_REGRESSION_TEST") is not None:
        # if test, only run five hours
        t_end = 5*t_export

    #Create solver
    solver_obj = solver2d.FlowSolver2d(mesh, bathymetry)

    options = solver_obj.options
    options.element_family = "dg-dg"
    options.polynomial_degree = 1
    options.coriolis_frequency =  coriolis_2d
    options.manning_drag_coefficient = manning_2d
    options.horizontal_velocity_scale = Constant(1.5)
    options.use_lax_friedrichs_velocity = True
    options.simulation_initial_date = start_date
    options.simulation_end_date = end_date
    options.simulation_export_time = t_export
    options.swe_timestepper_type = "DIRK22"
    options.swe_timestepper_options.use_semi_implicit_linearization = True
    options.timestep = dt
    options.fields_to_export = ["elev_2d", "uv_2d"]
    #options.fields_to_export_hdf5 = []

    #Coarse mesh, solve using a full LU decomposition as a preconditioner (???)

    options.swe_timestepper_options.solver_parameters = {
        "snes_type": "newtonls",
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "snes_monitor": None,
        #"snes_view": None,
        #"ksp_monitor_true_residual": None,
        #"snes_converged_reason": None,
        #"ksp_converged_reason": None
    }

    options.update(model_options)
    print_output(f"Exporting to {options.output_directory}")
    solver_obj.create_equations()

    #Forcings
    data_dir = "/home/epipremnum/Programs/OTPSnc/DATA/TPXO9v5a_nc/DATA"
    if not os.path.exists(data_dir):
        raise IOError(f"Data directory {data_dir} does not exist")
    forcing_constituents = ["Q1", "O1", "P1", "K1", "N2", "M2", "S2", "K2"]
    elev_tide_2d = Function(solver_obj.function_spaces.P1_2d, name = "Tidal elevation")
    tidecls = TPXOTidalBoundaryForcing
    tbnd = tidecls(
        elev_tide_2d, 
        start_date, 
        coord_system, 
        #data_dir = data_dir,
        uv_field=uv_field,
        data_dir='/home/epipremnum/Programs/OTPSnc/DATA/TPXO9v5a_nc/DATA',
        constituents=forcing_constituents,
        boundary_ids=[8] # Printed out at the end of generic_mesh_prep.py run
    )

    #Set time to zero for the tidal forcings
    tbnd.set_tidal_field(0.0)
    
    #Account for spinup
    bnd_time = Constant(0.0)
    if spinup:
        ramp_t = t_end
        elev_ramp = conditional(bnd_time < ramp_t, bnd_time / ramp_t, 1.0)
    else:
        elev_ramp = Constant(1.0)
    tide_elev_expr_2d = elev_ramp * elev_tide_2d

    #Boundary conditions for open ocean segments
    solver_obj.bnd_functions["shallow_water"] = {
        2: {"elev": tide_elev_expr_2d, "uv": Constant(as_vector([0,0]))},
    }

    def update_forcings(t):
        bnd_time.assign(t)
        tbnd.set_tidal_field(t)

    return solver_obj, start_date, update_forcings

########################################################################

# Setup solver
solver_obj, start_time, update_forcings = constructSolver(
    output_directory="runs/wall_A/outputs_spinup",
    spinup=True,
    start_date=datetime.datetime(2022, 1, 1, tzinfo=sim_tz),
    end_date=datetime.datetime(2022, 1, 15, tzinfo=sim_tz),
    fields_to_export=["elev_2d", "uv_2d"],
    fields_to_export_hdf5=["elev_2d", "uv_2d"],
    simulation_export_time=24 * 3600.0,
)

output_dir = solver_obj.options.output_directory
mesh2d = solver_obj.mesh2d
solver_obj.assign_initial_conditions()
update_forcings(0.0)

# Time integrate
tic = time_mod.perf_counter()
solver_obj.iterate(update_forcings=update_forcings)
toc = time_mod.perf_counter()
print_output(f"Total duration: {toc-tic:.2f} seconds")

from model_config import *
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
"""
with CheckpointFile("/home/epipremnum/Documents/tidal_dam_model/dungeness_contour/outputs_spinup/hdf5/Elevation2d_00014.h5", "r") as f:
    m = f.load_mesh("firedrake_default")
    elev_2d = f.load_function(m, "elev_2d")
fig, axes = plt.subplots(figsize=(6.4, 4.8))
triplot(mesh2d, axes=axes, boundary_kw={"linewidth": 0.5, "edgecolor": "k"})
norm = mcolors.CenteredNorm()
tc = tripcolor(elev_2d, axes=axes, cmap="RdBu_r", norm=norm)
cb = fig.colorbar(tc, ax=axes)
cb.set_label("Initial elevation (m)")
axes.axis(False)
plt.tight_layout()
imgfile = "north_sea_init.png"
print_output(f'Saving {imgfile}')
plt.savefig(imgfile, dpi=300)
"""

"""
#Tidal forcing
tidecls = TPXOTidalBoundaryForcing
tbnd = tidecls(
    elev_field, init_date, to_latlon, COORDSYS,
    uv_field=uv_field,
    data_dir='/home/epipremnum/Programs/OTPSnc/DATA/TPXO9v5a_nc/DATA',
    #constituents=[],
    boundary_ids=[2] # Can be 2 or 4, in this instance. This is from the mesh's unique ids
)
tbnd.set_tidal_field(0.0)

elev_out_filename = 'tmp/tidal_elev.pvd'
uv_out_filename = 'tmp/tidal_uv.pvd'
print('Saving to {:} {:}'.format(elev_out_filename, uv_out_filename))
elev_out = File(elev_out_filename)
uv_out = File(uv_out_filename)
for t in numpy.linspace(0, 12*3600., 49):
    tbnd.set_tidal_field(t)
    if elev_field.function_space().mesh().comm.rank == 0:
        print('t={:7.1f} elev: {:7.1f} uv: {:7.1f}'.format(t, norm(elev_field), norm(uv_field)))
    elev_out.write(elev_field)
    uv_out.write(uv_field)
"""

"""
tnci = uptide.tidal_netcdf.FESTidalInterpolator(b'/home/epipremnum/Programs/FES/aviso-fes/data/fes2014/ocean_tide.ini')
t = np.arange(0, 30*24*3600, 600)
tnci.set_time(t)
eta = tnci.get_val((mesh_lat,mesh_lon))
print(eta)
"""

#File('bathymetry_dungeness.pvd').write(bathymetry_2d)
