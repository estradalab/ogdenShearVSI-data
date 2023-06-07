from ufl import *
from dolfin import *
import numpy as np
import pandas as pd

temp = 0
filename_disp = 'disp_field_%i'%temp
filename_specimen = 'specimen_%i'%temp

#meshdir='../mesh/gmsh/holes_5mm_L%ia'%temp
#meshfilename = '/test%i.xdmf'%temp
# meshdir='../mesh/gmsh/noholes_L%ia'%temp
# meshfilename = '/testa.xdmf'
meshdir = '.'
meshfilename = '/ShearWavy_6.25MMDisp_Amp_0.6_tet.xdmf'


# datadir='../result/filtered/20220211-no_holes/raw_data'
datadir = '.'
NumStepList = [1]#,2,3]

mesh=Mesh() # Creating empty mesh object

with XDMFFile(meshdir+meshfilename) as infile: # create 'infile', an XDMFFile object imported from the file named in 'meshfilename'
    infile.read(mesh) # read in data from .xdmf file to mesh object called infile (used to access geometry/connectivity info about the mesh)


V = VectorFunctionSpace(mesh, "CG", 1) # create vector function space for 'mesh' composed of linear piece-wise continuous Lagrangian functions
W = FunctionSpace(mesh, "Lagrange", 1) # Create a scalar function space for 'mesh'composed of linear piece-wise continuous Lagrangian functions

# Define functions
u = Function(V) 

coordinate_data=pd.read_csv(datadir + '/node_pos.txt', header=None).to_numpy()[:,0:3] # make coordinate_data variable which contains node positions from node_pos.txt ()
# print(coordinate_data.shape)

displacement_data = [0]*len(NumStepList) # initialize 'displacement_data' as a zeros vector of dimension 1 by len(NumStepList) (currently yields a 1 by 1). should this have the same dimensions as coordinate_data.shape?

for count,tid in enumerate(NumStepList):
	displacement_data[count]=pd.read_csv(datadir + '/disp_'+str(tid)+'.txt', header=None).to_numpy()[:,0:3]
# print(count)
# print(tid)
# print(displacement_data)

# print(dof_to_vertex_map(V))
# print(dof_to_vertex_map(W))
# d2v_vector = dof_to_vertex_map(V)

# random_Class = V.dofmap()
# print(dir(random_Class))
a = V.dofmap().tabulate_local_to_global_dofs() # currently this is just print an array going from 0 to 3*NumberOfNodes in increments of 1 (in order) (i.e., simple mapping between local and global DOF)
# print(mesh)
# print(a[198065])
# print(len(a))
# exit()
coordinates_nodes = mesh.coordinates()

num_elements_data,_= coordinate_data.shape
num_nodes_data,_=coordinates_nodes.shape
print('num_elements =',num_elements_data)
print('num_nodes =',num_nodes_data)
# num_displacement_file,_=displacement_data.shape
# print('num_displacement in file =',num_displacement_file)

u_array=np.zeros((num_nodes_data*3, len(NumStepList))) # create array (1 x num_coordinate in mesh) of 0's to preallocate

for i in range(num_nodes_data): # looping over all nodes in data and finding near 
  X=coordinates_nodes[i] # import coordinates for node i in 'mesh'
  index=-1
  # print(coordinate_data[:,:].shape)
  # test = X*np.ones([num_elements_data, 1])
  # print(test.shape)
  # print(X.shape)
  _dist = coordinate_data[:,:]-X*np.ones([num_elements_data, 1]) 
  _dist = np.linalg.norm(_dist,axis = 1)
  index = np.argmin(_dist)
  print(index)
  # print(index)
  # print(_dist[index])
  # print('i=',i,' index=',index,' X=',X,' coordinate_data=',coordinate_data[index,:])
  
  for count,tid in enumerate(NumStepList):
    u_array[i*3:i*3+3, count]=np.reshape(displacement_data[count][index,:],(-1))
    # u_array[i*3:i*3+3, count]=np.reshape([1,2,3],(-1))



file_1 = XDMFFile(mesh.mpi_comm(),datadir + '/'+filename_disp+'.xdmf')
file_2 = HDF5File(MPI.comm_world, datadir + '/'+filename_specimen+'.h5', 'w')
file_2.write(u,'/mesh')


for count,tid in enumerate(NumStepList):
  u.vector().set_local(u_array[a, count])
  u.vector().apply("insert")



  # Save solution 
  file_1.write(u, tid)
  # #
  file_2.write(u,'/disp%i'%tid)

file_1.close()
file_2.close()

