import numpy as np
import gmsh
import json
import sys

#Read data
metadata = np.loadtxt('dungeness.asc', max_rows = 6, dtype = 'str')
shape = json.load(open('dungeness_coast_v2.geojson'))
features = shape['features']

data_properties = {}
for i in metadata:
    data_properties['%s' % (i[0])] = i[1]

#Memorise bathymetry metadata
ncols       = np.int(data_properties['ncols'])
nrows       = np.int(data_properties['nrows'])
xllcorner   = np.float(data_properties['xllcorner'])
yllcorner   = np.float(data_properties['yllcorner'])
cellsize    = np.float(data_properties['cellsize'])

#print("Cellsize: %s" % cellsize)

#Memorise Coast metadata
# note: only works for one coast and one ocean, i.e. a curve representing an island would break this functionality
for feature in features:
    if feature['properties']['Type'] == "ocean":
        ocean_geo = np.array(feature['geometry']['coordinates'][0])
    elif feature['properties']['Type'] == "coast":
        coast_geo = np.array(feature['geometry']['coordinates'][0])
    else:
        print("Unknown feature found in geojson")

#//Test for close points. 
#This in preparation for the entire curve that should start at the boundary of the ocean and move through the ocean points to the coast point, then move through the coast points and finish at the original start point
#It may be useful to turn this part into a function that works with more than two features.

oceanLastPoint_lat  = round(ocean_geo[-1][0],7)
oceanLastPoint_long = round(ocean_geo[-1][1],7)
oceanFrstPoint_lat  = round(ocean_geo[0][0],7)
oceanFrstPoint_long = round(ocean_geo[0][1],7)
coastFrstPoint_lat  = round(coast_geo[0][0],7)
coastFrstPoint_long = round(coast_geo[0][1],7)
coastLastPoint_lat  = round(coast_geo[-1][0],7)
coastLastPoint_long = round(coast_geo[-1][1],7)

#Test A is to see if the last ocean point is on top or very near the first coast point and, if it is, remove it.
test_A = oceanLastPoint_lat == coastFrstPoint_lat and oceanLastPoint_long == coastFrstPoint_long
if test_A:
    average_x = (ocean_geo[-1][0]+coast_geo[0][0])/2
    average_y = (ocean_geo[-1][1]+coast_geo[0][1])/2
    ocean_geo       = ocean_geo[:-1]
    coast_geo[0]    = [average_x, average_y]
#Test B is to see if the first ocean point is on top or very near the last coast point. If it is, remove the coastal point.
test_B = oceanFrstPoint_lat == coastLastPoint_lat and oceanFrstPoint_long == coastLastPoint_long
if test_B:
    average_x = (ocean_geo[0][0]+coast_geo[-1][0])/2
    average_y = (ocean_geo[0][1]+coast_geo[-1][1])/2
    ocean_geo[0]    = [average_x, average_y]
    coast_geo       = coast_geo[:-1]
#End of test for close points//

#Join features together
total_geo = np.concatenate((ocean_geo,coast_geo),axis=0)
#Calculate feature lengths
ocean_len = len(ocean_geo)
coast_len = len(coast_geo)
boundary_len = len(total_geo)

#//Define GMSH points

gmsh.initialize()
gmsh.option.setNumber('General.NumThreads', 20)

gmsh.model.add('coast_geometry')

points = []
lines = []

#Create points and lines in preparation for a loop. Note gmsh tags entities automatically from 1.
for i in range(boundary_len):
    points.append(gmsh.model.geo.addPoint(total_geo[i][0],total_geo[i][1],0,cellsize))
    #From the second point onwards, keep adding lines
    if i != 0:
        lines.append(gmsh.model.geo.addLine(i+1,i))
        if i == boundary_len-1:
            lines.append(gmsh.model.geo.addLine(1,i+1))
    
loop = [gmsh.model.geo.addCurveLoop(lines)]
surface = gmsh.model.geo.addPlaneSurface(loop)
gmsh.model.geo.synchronize()

#Define groups. 
#Ocean points does not include the second node shared with the coast, only the first; likewise, the coast points do not include the second node with the ocean, only the first

ocean_points    = gmsh.model.addPhysicalGroup(0,points[:ocean_len])
gmsh.model.setPhysicalName(0,ocean_points,"Ocean Points")
ocean_lines     = gmsh.model.addPhysicalGroup(1,lines[:ocean_len])
gmsh.model.setPhysicalName(1,ocean_lines,"Ocean Curve")
coast_points    = gmsh.model.addPhysicalGroup(0,points[ocean_len:])
gmsh.model.setPhysicalName(0,coast_points,"Coast points")
coast_lines     = gmsh.model.addPhysicalGroup(1,lines[ocean_len:])
gmsh.model.setPhysicalName(1,coast_lines,"Coast curve")

#Define fields.

distance = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(distance, "PointsList", gmsh.model.getEntitiesForPhysicalGroup(0,coast_points))
gmsh.model.mesh.field.setNumbers(distance, "CurvesList", gmsh.model.getEntitiesForPhysicalGroup(1,coast_lines))
gmsh.model.mesh.field.setNumber(distance, "Sampling", 600)

threshold = gmsh.model.mesh.field.add('Threshold')
gmsh.model.mesh.field.setNumber(threshold, "DistMax", 0.1)
gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.005)
gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
gmsh.model.mesh.field.setNumber(threshold, "SizeMin", 0.002)
gmsh.model.mesh.field.setNumber(threshold, "SizeMax", 0.004)

structure = gmsh.model.mesh.field.add("Structured")
gmsh.model.mesh.field.setString(structure, "FileName", r"dungeness_structured_bath_v2A.txt")
gmsh.model.mesh.field.setNumber(structure, "TextFormat", 1)

bath_threshold = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(bath_threshold, "DistMax", 66)
gmsh.model.mesh.field.setNumber(bath_threshold, "DistMin", 0)
gmsh.model.mesh.field.setNumber(bath_threshold, "InField", structure)
gmsh.model.mesh.field.setNumber(bath_threshold, "SizeMin", 0)
gmsh.model.mesh.field.setNumber(bath_threshold, "SizeMax", 66)

bath_math = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(bath_math, "F", r"(0.5 + 2*Sqrt(F4/66))*F2") #66 is the max figure in the bathymetry file; this normalizes the effect of the depth on the result.

gmsh.option.setNumber('General.NumThreads', 20)

gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

gmsh.model.mesh.field.setAsBackgroundMesh(bath_math)   
#Generate mesh
gmsh.model.mesh.generate(2)

#Save mesh
gmsh.option.setNumber("Mesh.SaveAll", 1)
gmsh.write('dungeness.msh')

#Show grid in gmsh popup
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#close gmsh
gmsh.finalize()
