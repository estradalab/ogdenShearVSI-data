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
# meshdir='../mesh/gmsh/noholes_L%ia'%temp
# meshfilename = '/testa.xdmf'
meshdir = '.'
meshfilename = '/ShearWavy_6.25MMDisp_Amp_0.6_tet.xdmf' # reading in .xdmf mesh file created in Step1_msh2xdmf.py which describes the mesh stored in the .inp file


# datadir='../result/filtered/20220211-no_holes/raw_data'
datadir = '.'
NumStepList = [3]#,2,3]

mesh=Mesh() # create empty FEniCS mesh object

with XDMFFile(meshdir+meshfilename) as infile: # create 'infile', an XDMFFile object imported from the file named in 'meshfilename'
    infile.read(mesh) # read in data from .xdmf file to mesh object called infile (used to access geometry/connectivity info about the mesh)
    

V = VectorFunctionSpace(mesh, "CG", 1) # create vector function space for 'mesh' composed of linear piece-wise continuous Lagrangian functions
W = FunctionSpace(mesh, "Lagrange", 1) # Create a scalar function space for 'mesh'composed of linear piece-wise continuous Lagrangian functions

# Define functions
u = Function(V) 

coordinate_elem_TXTdata=pd.read_csv(datadir + '/nod_u_GOOD.txt', header=None).to_numpy()[:,0:3] # make coordinate_data variable which contains node positions from node_pos.txt ()
# print(coordinate_data.shape)
# print("FROM THE MESH")
# print(mesh.coordinates()[:]) # This may very well be the issue
# print(mesh.coordinates()[:].shape)
# print("FROM SUS TXT FILE")
# print(coordinate_elem_TXTdata)
# print(coordinate_elem_TXTdata.shape)


displacement_elem_TXTdata = [0]*len(NumStepList) # initialize 'displacement_data' as a zeros vector of dimension 1 by len(NumStepList) (currently yields a 1 by 1). should this have the same dimensions as coordinate_data.shape?

for count,tid in enumerate(NumStepList):
	displacement_elem_TXTdata[count]=pd.read_csv(datadir + '/disp_'+str(tid)+'.txt', header=None).to_numpy()[:,0:3]
      
print(type(displacement_elem_TXTdata))
testing123 = np.concatenate(displacement_elem_TXTdata, axis=0)
first_array = np.array(displacement_elem_TXTdata[0])  # Select the first array
if isinstance(first_array,np.ndarray):
  print("is array!! woot woot oh yeah")
np.savetxt("nodeTXT.txt", coordinate_elem_TXTdata, delimiter=',', fmt='%f')    
np.savetxt("dispTXT.txt", testing123, delimiter=',', fmt='%f')
# print(count)
# print(tid)
# print(displacement_TXTdata)


# print(dof_to_vertex_map(V))
# print(dof_to_vertex_map(W))
# d2v_vector = dof_to_vertex_map(V) # Sid is going to try to find
# a = dof_to_vertex_map(V)
a = pd.read_csv('dof2vertex.txt', header=None).to_numpy()[:,0:3]

# random_Class = V.dofmap()
# print(dir(random_Class))
# a = V.dofmap().tabulate_local_to_global_dofs() # currently this is just print an array going from 0 to 3*NumberOfNodes in increments of 1 (in order) (i.e., simple mapping between local and global DOF)
# print(mesh)
# print(a[198065])
# print(len(a))
#print(a.shape)
# exit()
coordinates_nodes = mesh.coordinates()





# print(coordinates_nodes[0])
# print(coordinates_nodes[1])
# exit()

num_elements_data,_= coordinate_elem_TXTdata.shape
num_nodes_data,_=coordinates_nodes.shape
print('num_elements =',num_elements_data)
print('num_nodes =',num_nodes_data)

# num_displacement_file,_=displacement_data.shape
# print('num_displacement in file =',num_displacement_file)

u_array=np.zeros((num_nodes_data*3, len(NumStepList))) # create array (1 x num_coordinate in mesh) of 0's to preallocate
test_disp = np.zeros((num_nodes_data, len(NumStepList)))
test_index = np.zeros((num_nodes_data, len(NumStepList)))
test_count = np.zeros((num_nodes_data, len(NumStepList)))
for i in range(num_nodes_data): # looping over all nodes in data and finding nearest element
  X=coordinates_nodes[i] # import coordinates for node i in 'mesh'
  index=-1
  # print(coordinate_data[:,:].shape)
  # test = X*np.ones([num_elements_data, 1])
  # print(test.shape)
  # print(X.shape)
  _dist = coordinate_elem_TXTdata[:,:]-X*np.ones([num_elements_data, 1]) 
  _dist = np.linalg.norm(_dist,axis = 1) # calculating distance from X[i] to each node 
  index = np.argmin(_dist) # index of element located closest to the ith node
  test_disp[i] = np.min(_dist)
  test_index[i] = index
  # print(f'loop: {i}, index: {index}\n')
  # print(index)
  # print(index)
  # print(_dist[index])
  # print('i=',i,' index=',index,' X=',X,' coordinate_data=',coordinate_data[index,:])
  for count,tid in enumerate(NumStepList):
  #  print(index)
  #  print(coordinate_elem_TXTdata[index,:])
  #  print(displacement_elem_TXTdata[count][index,:],(-1))
  #  input("next")
    test_count[i] = count
    #pos_i = np.reshape(X,(-1))
    #u_array[i*3:i*3+3, count] = pos_i
    u_array[i*3:i*3+3, count] = np.reshape(displacement_elem_TXTdata[count][index,:],(-1))
  #  print(u_array[0:3])
   # print(u_array.shape)
   # print(u_array[i*3:i*3+3, count][:])
    #u_array[i*3:i*3+3, count]=np.reshape([1,2,3],(-1))

plt.hist(test_index, bins=60)
plt.show()
plt.savefig('histogramINDEX.png')

print("test_count")
print(sum(test_count))

# print("test")
# print(u_array[0:6])
# print(u_array.shape)

# checking u
check_disp = u.compute_vertex_values(mesh)
check_loc = mesh.coordinates()
np.savetxt("disp_u.txt", check_disp, delimiter=',', fmt='%f')    
np.savetxt("nod_u.txt", check_loc, delimiter=',', fmt='%f')



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

print(u)


file_1.close()
file_2.close()