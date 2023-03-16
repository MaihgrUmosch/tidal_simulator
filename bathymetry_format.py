def toBinary(a):
  l,m=[],[]
  for i in a:
    l.append(ord(i))
  for i in l:
    m.append(int(bin(i)[2:]))
  return m

import numpy as np

ascii_file = np.loadtxt('/home/epipremnum/Documents/tidal_dam_model/north_sea_example/GEBCO_10_Sep_2022_71c6879208fe/gebco_2022_n64.0723_s47.0215_w-13.9746_e12.7441.asc', skiprows=6)
_ascii_metadata = np.genfromtxt('/home/epipremnum/Documents/tidal_dam_model/north_sea_example/GEBCO_10_Sep_2022_71c6879208fe/gebco_2022_n64.0723_s47.0215_w-13.9746_e12.7441.asc',max_rows=6,dtype=[('cust_str','S10'), ('cust_flt','f16')])

ascii_metadata = {}

for data in _ascii_metadata:
    if type(data[1]) == np.float128:
        ascii_metadata[data[0]] = data[1]
    else:
        ascii_metadata[data[0]] = toBinary(data[1])
#print(_ascii_metadata)
#print(len(ascii_file[0]))


# Generate arrays of metadata and put into array. These are saved at the top of the file
coordinates_O   = [ascii_metadata[b'xllcorner'],ascii_metadata[b'yllcorner'],0]
distance_D      = [ascii_metadata[b'cellsize'],ascii_metadata[b'cellsize'],1]
number_n        = [ascii_metadata[b'ncols'],ascii_metadata[b'nrows'],1]
structured_metadata = [coordinates_O, distance_D, number_n]

structured_array = []

for i in range(0,int(ascii_metadata[b'ncols'])):
    for j in range(0,int(ascii_metadata[b'nrows'])):
        structured_array.append(ascii_file[int(ascii_metadata[b'nrows'])-j-1,i])

f = open('north_sea_structured_bath.txt', 'wb')
for i in range(0,3):
    if i<2: 
        f.write(bytes('%f %f %f\n' % (structured_metadata[i][0],structured_metadata[i][1],structured_metadata[i][2]),encoding='utf8'))
    else: 
        f.write(bytes('%i %i %i\n' % (structured_metadata[i][0],structured_metadata[i][1],structured_metadata[i][2]),encoding='utf8'))
for i in range(len(structured_array)):
    f.write(bytes('%i\n' % -structured_array[i],encoding='utf8'))

f.close()