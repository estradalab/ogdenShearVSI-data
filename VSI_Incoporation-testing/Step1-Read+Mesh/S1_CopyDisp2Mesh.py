import os
from ufl import *
from dolfin import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

temp = 0
filename_disp = 'disp_field_%i'%temp
filename_specimen = 'specimen_%i'%temp

#meshdir='../mesh/gmsh/holes_5mm_L%ia'%temp
#meshfilename = '/test%i.xdmf'%temp
meshdir = '.'
meshfilename = '/sq-8mm_sin-per-4_sin-amp-2mm_tet.xdmf' # reading in .xdmf mesh file created in Step1_msh2xdmf.py which describes the original MATLAB .inp mesh in a mesh.io python mesh

datadir = '.'
TimeStep = [2]

mesh_py=Mesh() # create empty FEniCS mesh object

with XDMFFile(meshdir+meshfilename) as infile: # create 'infile', an XDMFFile object imported from 'meshfilename'
    infile.read(mesh_py) # read in data from .xdmf file to mesh object called infile (used to access geometry/connectivity info about the mesh) 

V = VectorFunctionSpace(mesh_py, "CG", 1) # create vector function space for 'mesh' composed of linear piece-wise continuous Lagrangian functions
W = FunctionSpace(mesh_py, "Lagrange", 1) # Create a scalar function space for 'mesh'composed of linear piece-wise continuous Lagrangian functions

# Define functions
u = Function(V) 

coordinate_elem_ABAQUS=pd.read_csv(datadir + '/X_061623_validationtest2.csv', header=None).to_numpy()[:,0:3] # make coordinate_data variable which contains node positions from ABAQUS .txt data
print(coordinate_elem_ABAQUS.shape)
check_pos = mesh_py.coordinates()
check_disp = u.compute_vertex_values(mesh_py)
filename_geom = "GeometryCheck.png"
if os.path.exists(filename_geom):
    os.remove(filename_geom)
print("x from u",max(check_pos[:,0]))
print("y from u",max(check_pos[:,1]))
print("z from u",max(check_pos[:,2]))
print("x from .txt",max(coordinate_elem_ABAQUS[:,0])) # check max amplidue 
print("y from .txt",max(coordinate_elem_ABAQUS[:,1])) # check max amplidue
print("z from .txt",max(coordinate_elem_ABAQUS[:,2])) # check max amplidue
plt.xlabel("x (mm)")
plt.ylabel("y (mm)")
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
axs[0].scatter(coordinate_elem_ABAQUS[:,0],coordinate_elem_ABAQUS[:,1],marker=",",alpha=0.1,color='blue',label='From .txt ABAQUS results file')
axs[0].set_title('From .txt ABAQUS results file')
axs[1].scatter(check_pos[:,0],check_pos[:,1],marker=",",alpha=0.1,color='red',label='From .inp file')
axs[1].set_title('From .inp file')
# plt.legend()
plt.savefig(filename_geom)
plt.close()

displacement_elem_ABAQUS = [0]*len(TimeStep) # initialize 'displacement_data' as a zeros vector of dimension 1 by len(NumStepList) (currently yields a 1 by 1).

for count,tid in enumerate(TimeStep):
	# displacement_elem_ABAQUS[count]=pd.read_csv(datadir + '/disp_'+str(tid)+'.txt', header=None).to_numpy()[:,0:3]
  displacement_elem_ABAQUS[count]=pd.read_csv(datadir + '/U_061623_validationtest_2.csv', header=None).to_numpy()[:,0:3]

DOF2VertexMap = pd.read_csv('dof2vertex.txt', header=None).to_numpy()[:,0:3]

coordinates_nodes = mesh_py.coordinates()

num_elements_data,_= coordinate_elem_ABAQUS.shape
num_nodes_data,_=coordinates_nodes.shape
print('num_elements =',num_elements_data)
print('num_nodes =',num_nodes_data)

u_array=np.zeros((num_nodes_data*3, len(TimeStep))) # create array (1 x num_coordinate in mesh) of 0's to preallocate
for i in range(num_nodes_data): # looping over all nodes in data 
  # Find nearest MATLAB node for each Python node
  X=coordinates_nodes[i] # import coordinates for node i in 'mesh'
  index=-1
  _dist = coordinate_elem_ABAQUS[:,:]-X*np.ones([num_elements_data, 1]) 
  _dist = np.linalg.norm(_dist,axis = 1) # calculating distance from X[i] to each node 
  index = np.argmin(_dist) # index of element located closest to the ith node
  
  # Write displacement data from nearest node to Python node
  for count,tid in enumerate(TimeStep):
    u_array[i*3:i*3+3, count] = np.reshape(displacement_elem_ABAQUS[count][index,:],(-1))

# checking u in MATLAB
check_disp = u.compute_vertex_values(mesh_py)
check_loc = mesh_py.coordinates()
print(check_disp.shape,check_loc.shape)
np.savetxt("disp_u_61923.txt", check_disp, delimiter=',', fmt='%f')    
np.savetxt("nod_u_61923.txt", check_loc, delimiter=',', fmt='%f')

# Write .xdmf/.h5 files with 
file_1 = XDMFFile(mesh_py.mpi_comm(),datadir + '/'+filename_disp+'.xdmf')
file_2 = HDF5File(MPI.comm_world, datadir + '/'+filename_specimen+'.h5', 'w')
file_2.write(u,'/mesh')

for count,tid in enumerate(TimeStep):
  # Adding displacements to 'u'
  print(DOF2VertexMap.shape)
  u.vector().set_local(u_array[DOF2VertexMap, count])
  u.vector().apply("insert")

  # Save solution to .xdmf/.h5
  file_1.write(u, tid)
  file_2.write(u,'/disp%i'%tid)

# Check u
print("Joseph")
print(u)
print("Beckett")

file_1.close()
file_2.close()