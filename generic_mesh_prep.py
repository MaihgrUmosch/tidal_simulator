from multiprocessing.resource_sharer import stop
import numpy as np
import gmsh
import json
import sys

from sympy import ordered

#Read data
metadata = np.loadtxt(
    '/home/epipremnum/Documents/tidal_dam_model/dungeness_contour/dungeness.asc', 
    max_rows = 6, 
    dtype = 'str'
    )

bath = np.loadtxt(
    '/home/epipremnum/Documents/tidal_dam_model/dungeness_contour/dungeness.asc',
    skiprows = 6,
    dtype = 'int'
    )
max_depth = -np.min(bath)
shape = json.load(
    open(
        '/home/epipremnum/Documents/tidal_dam_model/dungeness_contour/dungeness_coast_wall_A_v3.geojson'
        )
        )
features = shape['features']

meshOptions_file = np.loadtxt('/home/epipremnum/Documents/tidal_dam_model/dungeness_contour/meshOptions_dungeness_wall.txt', dtype = 'str')

data_properties = {}
for i in metadata:
    data_properties['%s' % (i[0])] = i[1]

meshOptions = {}
for i in meshOptions_file:
    meshOptions['%s' % (i[0])] = float(i[1])

print(meshOptions)

#Memorise bathymetry metadata
ncols       = int(data_properties['ncols'])
nrows       = int(data_properties['nrows'])
xllcorner   = float(data_properties['xllcorner'])
yllcorner   = float(data_properties['yllcorner'])
cellsize    = float(data_properties['cellsize'])

#print("Cellsize: %s" % cellsize)

#Memorise Coast metadata
# note: only works for one coast and one ocean, i.e. a curve representing an island would break this functionality
ocean_geo = []
coast_geo = []
islnd_geo = []
for feature in features:
    if feature['properties']['Type'] == 200 or feature['properties']['Type'] == '200':
        ocean_geo.append(feature['geometry']['coordinates'][0])
    elif feature['properties']['Type'] == 100 or feature['properties']['Type'] == '100':
        lines = feature['geometry']['coordinates']
    
        for line in lines:
                coast_geo.append(line)
        else:
            coast_geo.append(line)
    elif feature['properties']['Type'] == 300 or feature['properties']['Type'] == '300':
        if len(feature['geometry']['coordinates']) != 0:
            islnd_geo.append(feature['geometry']['coordinates'][0])
    else:
        print("Unknown feature found in geojson")

# Arrays need sorting. Imported coast features must be ordered as merges can cause
# lines to go out of order. Start and end points of arrays should not be duplicated. 
# fThis should be optional.
# Island features are OK. 
def remove_duplicate_points(geo_array):
    #removes adjacent double up points in array of geo coords
    geo_array_len = len(geo_array)
    output_geo_array = []
    for i in range(0,geo_array_len-1):
        if i == 0:
            continue
        current_coords = np.round(geo_array[i],6)
        prv_coords = np.round(geo_array[i-1],6)
        nxt_coords = np.round(geo_array[i+1],6)
        prv_coords_test = prv_coords == current_coords
        nxt_coords_test = nxt_coords == current_coords
        if i > 0:
            if prv_coords_test.all() and nxt_coords_test.all():
                geo_array.pop(i)
                geo_array.pop(i+1)
            elif prv_coords_test.all():
                geo_array.pop(i)
            elif nxt_coords_test.all():
                geo_array.pop(i+1)
        if i >= len(geo_array)-2:
            return geo_array
    
    
            

def remove_duplicates(array_of_lines):
    array = []
    num_arrays = len(array_of_lines)
    popnum = 0
    for i in range(num_arrays):
        if i != num_arrays-popnum:
            ith_array = np.round(array_of_lines[i],6)
            for j in range(1,num_arrays-i-popnum,1):
                if i+j < num_arrays-popnum:
                    jth_array = np.round(array_of_lines[i+j],6)
                    #NumPy behaves strangely when testing truths in two arrays. Sometimes output is bool, sometimes it is an (NumPy) array
                    ts = jth_array == ith_array
                    if hasattr(ts,'all'):
                        if ts.all():
                            array_of_lines.pop(i+j)
                            popnum += 1
                    else:
                        if ts:
                            array_of_lines.pop(i+j)
                            popnum += 1
        else:
            break
    return array_of_lines

def order_lines(array_of_lines):
    num_lines = len(array_of_lines)
    ordered_lines = []
    ordered_indexes = []
    for i in range(num_lines):
        if len(ordered_lines) == 0:
            line = array_of_lines[i]
            ordered_lines = [line]
            ordered_indexes = [i]
        elif len(ordered_lines) == num_lines:
            break
        for j in range(num_lines):
            if i == j:
                continue
            elif len(ordered_lines) == num_lines:
                break
            else:
                i_frst_coord = np.round(ordered_lines[0][0],4)
                i_last_coord = np.round(ordered_lines[-1][-1],4)
                j_frst_coord = np.round(array_of_lines[j][0],4)
                j_last_coord = np.round(array_of_lines[j][-1],4)
                i_then_j = i_last_coord == j_frst_coord
                j_then_i = j_last_coord == i_frst_coord
                if i_then_j.all():
                    ordered_lines.append(array_of_lines[j][1:])
                    ordered_indexes.append(j)
                elif j_then_i.all():
                    ordered_lines.reverse()
                    ordered_lines.append(array_of_lines[j][:-1])
                    ordered_lines.reverse()
                    ordered_indexes.reverse()
                    ordered_indexes.append(-j) #minus indicates line has been reversed
                    ordered_indexes.reverse()                       
        """num_ordered_lines = len(ordered_lines)
        lines_left = num_lines - num_ordered_lines
        last_coord_lat, last_coord_lon = ordered_lines[-1][-1]
        for j in range(lines_left):
            if i+j >= lines_left-1:
                break
            else:
                line = array_of_lines[i+j+1]
                if np.round(last_coord_lat,6) == np.round(line[0][0],6) and np.round(last_coord_lon,6) == np.round(line[0][1],6):
                    #Tolerance too high or points too far apart.
                    ordered_lines.append(line[1:])
                    ordered_indexes.append(i+j+1)
                    break
    """
    unordered_lines = []
    for i in range(num_lines):
        if i in ordered_indexes:
            continue
        else:
            unordered_lines.append(array_of_lines[i])
    return ordered_lines, ordered_indexes

#Remove possible coast duplicates
coast_geo = remove_duplicates(coast_geo)

#Set up entire border geometry
border_geo = coast_geo
border_tag = []
for i in range(len(coast_geo)):
    border_tag.append('100')
for ocean_line in ocean_geo:
    border_geo.append(ocean_line)
    border_tag.append('200')

ordered_border_geo, indexes = order_lines(border_geo)


def get_tag(index):
    if index < 0:
        return border_tag[-index]
    else:
        return border_tag[index]


# Once arrays are imported, the coast and ocean arrays must also be ordered.
# Tests to ensure that start and end points are not doubled.

# Requirement: ocean and coast lines are identifiable within the array.
# This could be done with the length of the joined up arrays, assuming that coasts and ocean boundaries have been kept separately

def flatten_array(nested_array):
    array = []
    for sub_array in nested_array:
        for element in sub_array:
            array.append(element)
    return array

def flatten_tagged_array(nested_array, indexes):
    array = []
    tag_array = []
    len_input_array = len(nested_array)
    for i in range(len_input_array):
        for element in nested_array[i]:
            tag_array.append(get_tag(indexes[i]))
            array.append(element)
    return array, tag_array

def check_for_line_intersections(array_of_lines):
    #https://en.wikipedia.org/wiki/Line-line_intersection
    a = array_of_lines
    num = len(a)
    for i in range(0,num):
        x1 = a[i][0][0]
        y1 = a[i][0][1]
        x2 = a[i][1][0]
        y2 = a[i][1][1]
        for j in range(i+1,num):
            x3 = a[j][0][0]
            y3 = a[j][0][1]
            x4 = a[j][1][0]
            y4 = a[j][1][1]
            #Determinant calculation of intersection
            denom = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)
            if denom != 0:
                t_num = (x1-x3)*(y3-y4)-(y1-y3)*(x3-x4)
                u_num = (x1-x3)*(y1-y2)-(y1-y3)*(x1-x2)
                t_test = 0 <= t_num <= denom
                u_test = 0 <= u_num <= denom
                if t_test and u_test:
                    t = t_num/denom
                    u = u_num/denom
                    Px,Py = [x1+t*(x2-x1),y1+t*(y2-y1)]
                    line_end_test_1 = np.round([Px,Py],6) != np.round(a[i][0],6) 
                    line_end_test_2 = np.round([Px,Py],6) != np.round(a[i][1],6)

                    if line_end_test_1.all() and line_end_test_2.all():
                        print("Intersection found between lines: " + str(a[i]) + " and " + str(a[j]) + ". The intersection P is " + str([Px,Py]))
        
    return

ordered_border_geo, tag_array = flatten_tagged_array(ordered_border_geo, indexes)
test_overlap = np.round(ordered_border_geo[-1],3) == np.round(ordered_border_geo[0],3)
if test_overlap.all():
    ordered_border_geo = ordered_border_geo[:-1]
    tag_array = tag_array[:-1]
lengths = {}
piece_num = 0
for i in range(len(tag_array)-1):
    if i == 0:
        key = tag_array[i]+"_"+str(piece_num)
        lengths[key] = 0
    if tag_array[i] == tag_array[i+1]:
        key = tag_array[i]+"_"+str(piece_num)
        lengths[key] += 1
    else:
        key = tag_array[i]+"_"+str(piece_num)
        lengths[key] += 1
        piece_num += 1
        next_key = tag_array[i+1]+"_"+str(piece_num)
        lengths[next_key] = 0

#//Test for close points. 
#This in preparation for the entire curve that should start at the boundary of the ocean and move through the ocean points to the coast point, then move through the coast points and finish at the original start point
#It may be useful to turn this part into a function that works with more than two features.

"""

oceanLastPoint_lat  = np.round(ocean_geo[-1][0],7)
oceanLastPoint_long = np.round(ocean_geo[-1][1],7)
oceanFrstPoint_lat  = np.round(ocean_geo[0][0],7)
oceanFrstPoint_long = np.round(ocean_geo[0][1],7)
coastFrstPoint_lat  = np.round(coast_geo[0][0],7)
coastFrstPoint_long = np.round(coast_geo[0][1],7)
coastLastPoint_lat  = np.round(coast_geo[-1][0],7)
coastLastPoint_long = np.round(coast_geo[-1][1],7)


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
"""
boundary_len = len(ordered_border_geo)


#//Define GMSH points

gmsh.initialize()
gmsh.option.setNumber('General.NumThreads', 20)

gmsh.model.add("Sea_mesh")

#Add islands
num_islnds = len(islnd_geo)
all_islnd_points = []
all_islnd_lines = []
all_islnd_curve_loops = []
_all_islnd_lines = []
point_index = 0
for islnd_index in range(0,num_islnds):
    single_islnd_geo = remove_duplicate_points(islnd_geo[islnd_index])
    islnd_len = len(single_islnd_geo)
    islnd_points = []
    islnd_lines = []
    _islnd_lines = []
    loop_test = np.round(single_islnd_geo[0],6) != np.round(single_islnd_geo[-1],6) #Island doesn't form complete loop

    for i in range(0,islnd_len):
        islnd_geo_i = single_islnd_geo[i]
        if i == 0:
            islnd_points.append(gmsh.model.geo.addPoint(islnd_geo_i[0],islnd_geo_i[1],0,cellsize))
            first_tag = islnd_points[0]
        elif 0 < i < islnd_len-1:
            islnd_points.append(gmsh.model.geo.addPoint(islnd_geo_i[0],islnd_geo_i[1],0,cellsize))
            islnd_lines.append(gmsh.model.geo.addLine(islnd_points[i],islnd_points[i-1]))
            _islnd_lines.append([islnd_geo_i,single_islnd_geo[i-1]])
        elif i == islnd_len-1:
            if loop_test.any():
                islnd_points.append(gmsh.model.geo.addPoint(islnd_geo_i[0],islnd_geo_i[1],0,cellsize))
                islnd_lines.append(gmsh.model.geo.addLine(islnd_points[i],islnd_points[i-1]))
                _islnd_lines.append([single_islnd_geo[i-1],islnd_geo_i])
                islnd_lines.append(gmsh.model.geo.addLine(first_tag,islnd_points[i]))
                _islnd_lines.append([islnd_geo_i,single_islnd_geo[0]])
            else:
                islnd_lines.append(gmsh.model.geo.addLine(first_tag,islnd_points[i-1]))
                _islnd_lines.append([single_islnd_geo[i-1],single_islnd_geo[0]])
         
    all_islnd_points.append(islnd_points)
    all_islnd_lines.append(islnd_lines)
    _all_islnd_lines.append(_islnd_lines)

    loop = gmsh.model.geo.addCurveLoop(islnd_lines)
    all_islnd_curve_loops.append(loop)
    
    #surface = gmsh.model.geo.addPlaneSurface([loop])

gmsh.model.geo.synchronize()

all_islnd_points_group = []
all_islnd_lines_group = []

#Add borders
#gmsh.model.add('border_geo')

border_points = []
border_lines = []
_border_lines = []
border_loop_test = np.round(ordered_border_geo[0],6) != np.round(ordered_border_geo[-1],6)
#Create points and lines in preparation for a loop. Note gmsh tags entities automatically from 1.
#ordered_border_geo = np.round(ordered_border_geo,6)
for i in range(boundary_len):
    if i == 0:
        border_points.append(gmsh.model.geo.addPoint(ordered_border_geo[i][0],ordered_border_geo[i][1],0,cellsize))
    #From the second point onwards, keep adding lines
    if 0 < i < boundary_len-1:
        border_points.append(gmsh.model.geo.addPoint(ordered_border_geo[i][0],ordered_border_geo[i][1],0,cellsize))
        border_lines.append(gmsh.model.geo.addLine(border_points[i],border_points[i-1]))
        _border_lines.append([ordered_border_geo[i],ordered_border_geo[i-1]])
    elif i == boundary_len-1:
        if border_loop_test.any():
            border_points.append(gmsh.model.geo.addPoint(ordered_border_geo[i][0],ordered_border_geo[i][1],0,cellsize))
            border_lines.append(gmsh.model.geo.addLine(border_points[i],border_points[i-1]))
            _border_lines.append([ordered_border_geo[i],ordered_border_geo[i-1]])
            border_lines.append(gmsh.model.geo.addLine(border_points[0],border_points[i]))
            _border_lines.append([ordered_border_geo[0],ordered_border_geo[i]])
        else:
            border_lines.append(gmsh.model.geo.addLine(border_points[0],border_points[i-1]))
            _border_lines.append([ordered_border_geo[0],ordered_border_geo[i-1]])

#gmsh.model.geo.synchronize()

#Test for intersecting lines:
_flat_all_islnd_lines = flatten_array(_all_islnd_lines)
for i in range(0,num_islnds):
    _border_lines.append(_flat_all_islnd_lines[i])
#print("Checking for line intersections...")
#check_for_line_intersections(_border_lines) #time-consuming
#print("Done")

border_loop = [gmsh.model.geo.addCurveLoop(border_lines)]

surface_loops = [border_loop[0]]

for i in range(0,num_islnds):
    surface_loops.append(all_islnd_curve_loops[i])
#surface = gmsh.model.geo.addPlaneSurface(surface_loops)             
surface = gmsh.model.geo.addPlaneSurface(surface_loops)
gmsh.model.geo.synchronize()

#Add physical groups
for islnd_index in range(0,num_islnds):
    new_points_group = gmsh.model.addPhysicalGroup(0,all_islnd_points[islnd_index])
    new_lines_group = gmsh.model.addPhysicalGroup(1,all_islnd_lines[islnd_index])
    gmsh.model.setPhysicalName(0,new_points_group,'300'+str(islnd_index))
    gmsh.model.setPhysicalName(1,new_lines_group,'300'+str(islnd_index))
    all_islnd_points_group.append(new_points_group)
    all_islnd_lines_group.append(new_lines_group)
    
gmsh.model.geo.synchronize()

#Define groups. 
#Ocean points does not include the second node shared with the coast, only the first; likewise, the coast points do not include the second node with the ocean, only the first

length_covered = 0
coast_points = []
coast_lines = []
ocean_points = []
ocean_lines = []
for piece in lengths:
    tag = piece[0:3]
    length = lengths[piece]

    new_points = gmsh.model.addPhysicalGroup(0,border_points[length_covered:length_covered+length])
    new_lines = gmsh.model.addPhysicalGroup(1,border_lines[length_covered:length_covered+length])
    gmsh.model.setPhysicalName(0,new_points,piece)
    gmsh.model.setPhysicalName(1,new_lines,piece)    
    if tag == '100':
        coast_points.append(new_points)
        coast_lines.append(new_lines)
    elif tag == '200':
        ocean_points.append(new_points)
        ocean_lines.append(new_lines)
    length_covered += length

gmsh.model.geo.synchronize()



    
"""
ocean_points    = gmsh.model.addPhysicalGroup(0,points[:ocean_len])
gmsh.model.setPhysicalName(0,ocean_points,"Ocean Points")
ocean_lines     = gmsh.model.addPhysicalGroup(1,lines[:ocean_len])
gmsh.model.setPhysicalName(1,ocean_lines,"Ocean Curve")
coast_points    = gmsh.model.addPhysicalGroup(0,points[ocean_len:])
gmsh.model.setPhysicalName(0,coast_points,"Coast points")
coast_lines     = gmsh.model.addPhysicalGroup(1,lines[ocean_len:])
gmsh.model.setPhysicalName(1,coast_lines,"Coast curve")
"""
#Define fields.

distance = gmsh.model.mesh.field.add("Distance")
entities_points = []
entities_lines = []
num_coasts = len(coast_points)
for i in range(0,num_coasts+num_islnds):
    if i < num_coasts:
        relevant_points = coast_points[i]
        relevant_lines = coast_lines[i]
    else:
        relevant_points = all_islnd_points_group[i-num_coasts]
        relevant_lines = all_islnd_lines_group[i-num_coasts]
    ith_points = gmsh.model.getEntitiesForPhysicalGroup(0,relevant_points)
    ith_lines = gmsh.model.getEntitiesForPhysicalGroup(1,relevant_lines)
    for point_index in range(0,len(ith_points)):
        entities_points.append(ith_points[point_index])
    for line_index in range(0,len(ith_lines)):
        entities_lines.append(ith_lines[line_index])


gmsh.model.mesh.field.setNumbers(distance, "PointsList", entities_points)
gmsh.model.mesh.field.setNumbers(distance, "CurvesList", entities_lines)
gmsh.model.mesh.field.setNumber(distance, "Sampling", 500)

threshold = gmsh.model.mesh.field.add('Threshold')
#Set distance from coast that there are small mesh sizes
gmsh.model.mesh.field.setNumber(threshold, "DistMax", meshOptions['distMax']) #After this distance, mesh size will be normal
gmsh.model.mesh.field.setNumber(threshold, "DistMin", meshOptions['distMin']) #Up to this distance, mesh size will be SizeMin, after this and up to DistMax, mesh size will be graded to SizeMax
gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
gmsh.model.mesh.field.setNumber(threshold, "SizeMin", meshOptions['sizeMin'])
gmsh.model.mesh.field.setNumber(threshold, "SizeMax", meshOptions['sizeMax'])


#Then, flex the mesh to use a finer mesh in shallower depths
structure = gmsh.model.mesh.field.add("Structured")
gmsh.model.mesh.field.setString(structure, "FileName", r"/home/epipremnum/Documents/tidal_dam_model/north_sea_example/north_sea_structured_bath.txt")
gmsh.model.mesh.field.setNumber(structure, "TextFormat", 1000)

bath_threshold = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(bath_threshold, "DistMax", max_depth)
gmsh.model.mesh.field.setNumber(bath_threshold, "DistMin", 0)
gmsh.model.mesh.field.setNumber(bath_threshold, "InField", structure)
gmsh.model.mesh.field.setNumber(bath_threshold, "SizeMin", 0)
gmsh.model.mesh.field.setNumber(bath_threshold, "SizeMax", max_depth)

bath_math = gmsh.model.mesh.field.add("MathEval")
expr = r"(0.5 + 2*Sqrt(F4/"+str(max_depth)+"))*F2"
gmsh.model.mesh.field.setString(bath_math, "F", expr) #this normalizes the effect of the depth on the result.

gmsh.option.setNumber('General.NumThreads', 20)

gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

gmsh.model.mesh.field.setAsBackgroundMesh(bath_math)   
#Generate mesh"""
gmsh.model.mesh.generate(2)

#Save mesh
gmsh.option.setNumber("Mesh.SaveAll", 1)
gmsh.write('/home/epipremnum/Documents/tidal_dam_model/dungeness_contour/dungeness_wall_A_v2.msh')
#print(expr)

#Show grid in gmsh popup
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

#close gmsh
gmsh.finalize()

print("Ocean line ids: %s" % ocean_lines)
