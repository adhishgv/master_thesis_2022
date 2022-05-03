# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 10:14:48 2022

"""

import numpy as np      # for operations and data types
from scipy.interpolate import RegularGridInterpolator   # for 3D interpolation
import meshio   # read/write meshes
import random
import csv


inputfile="mwe_input.vtu"   # generateStructuredMesh -e hex --lx 1 --ly 1 --lz 1 --nx 30 --ny 30 --nz 30 -o mwe_input.vtu
                            # generateStructuredMesh -e quad --lx 1 --ly 1 --nx 30 --ny 30  -o mwe_input.vtu


#outputfile="mwe_output.vtu"
n=30
nnn = np.array([n, n, n])  # divisions of data mesh (to be appended to input mesh) in x-, y-, z-direction
xyz = np.array([0,1, 0,1, 0,1])    # start end (x0,x1, y0,y1, z0,z1)
k = np.linspace(1,1000,1000)          #Permeability values
my_data = np.random.standard_normal(size=(nnn[0], nnn[1], nnn[2]))


H = 1
density = 1000
beta = 0.01
dT = 1
porosity = 0.1
g = 9.8
viscosity = 10e-3
FTC = 0.65
FSH = 4200

#Create csv file to store k values
with open('key.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Block 1', 'Block 2', 'Block 3', 'Block 4', 'Block 5', 'Block 6', 'Block 7', 'Block 8', 'Block 9'])
    
#Loop for generating output vtu files
for num in range(0,2):
    key = np.zeros([3,3])
    outputfile="output_%s.vtu"%(num+1)
    #Initialize and assign k values
    k_data = np.zeros(np.shape(my_data))
    k_indices = random.sample(range(len(k)),9) #indices of random k values
    key = k[k_indices]
    #Initialize coordinates for blocks
    x=np.array([0,10,20,30])
    y=np.array([0,10,20,30])
    #x1 = 0
    #x2=1

    z1=0
    #Nested loops for assigning k values to cells
    for zy in range(0,3):
        y1=zy
        y2=y1+1
        for zx in range(0,3):
            x1=zx
            x2=x1+1
            for i in range(x[x1],x[x2]):
                for j in range(y[y1],y[y2]):
                    k_data[i,j,0]=k[k_indices[z1]]
            
            z1=z1+1
            
        
    '''     
    z=0 
    for i in range(0,n):
        for j in range(0,n):
            for l in range(0,n):
                k_data[i,j,l] = k[k_indices[z]]
            z=z+1
            
    #print(k_data)        
    '''

    # vm data in domain: [X1R,X2R] [Y1R,Y2R] [ZBAS,Z2R]
    NXW = nnn[0]
    NYW = nnn[1]
    NZW = nnn[2]
    x = np.linspace(0.0, xyz[1]-xyz[0], NXW) # shifted to x=0 as mesh does
    y = np.linspace(0.0, xyz[3]-xyz[2], NYW) # shifted to y=0 as mesh does
    z = np.linspace(0.0, xyz[5]-xyz[4], NZW) # shifted to z=0 as mesh does
    #print("Data grid nxw, nyw, nzw: ", NXW, NYW, NZW)
    #print("Data grid x1,x2, y1,y2, zb,zr: ", xyz)

    # READ IN MESH AND ADD INTERPOLATED DATA
    mesh = meshio.read(inputfile)
    #print(f'Input mesh: {len(mesh.points)} points') # array
    #print(f'Input mesh: {len(mesh.cells[0].data)} cells')   # dictionary
    #print(mesh)

    # add interpolated data mesh to input mesh ([0,LX]  [0,LY] [0,LZ])
    interpolating_function = RegularGridInterpolator((x, y, z),
                                                     k_data[0:NXW, 0:NYW, 0:NZW],
                                                     method='nearest',
                                                     bounds_error=False,
                                                     fill_value=None)   # method='nearest', 'linear'

    # TODO there may be multiple data on multiple cell blocks
    cells_count = len(mesh.cells[0].data)   # MeshIO < 5.0   cells_count = len(mesh.cells[0][1])
    centerpoints=np.zeros((cells_count, 3))   # evaluate interpolation at center points
    for cell_index, cellblock_cell in enumerate(mesh.cells[0].data):    # MeshIO < 5.0   for cell_index, cellblock_cell in enumerate(mesh.cells[0][1]):
        centerpoints[cell_index] = (np.sum(mesh.points[cellblock_cell], axis=0) / len(cellblock_cell))

    my_data_interpolation = interpolating_function(centerpoints)

    my_data_min = np.min(my_data_interpolation)
    my_data_max = np.max(my_data_interpolation)
    print("Data appended to input mesh are within range {}...{}".format(my_data_min, my_data_max))

    mesh.cell_data['k_data'] = [np.array(my_data_interpolation)]   # append data

    # write to file
    meshio.write(outputfile, mesh)
    
    #write k values to csv file
    with open('key.csv', 'a', newline='') as file:
         writer = csv.writer(file)
         writer.writerow(key)
    
    #Delete used k values to avoid repetition
    k = np.delete(k,k_indices)
    
