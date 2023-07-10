# Overview: read quad tet mesh defined by <filename>.inp and write simplified linear tet mesh to <filename>.xdmf/.h5

import meshio
from dolfin import *
# from dolfin_adjoint import *
import os
import numpy as np

def gmsh2vtk_hex(P):
    return [P[0], P[3], P[2], P[1], P[4], P[5], P[6], P[7]]
	
def vtk2dolfin_hex(P):
    return [P[0], P[1], P[3], P[2], P[4], P[5], P[7], P[6]]

def dolfin2vtk_hex(P):
    return [P[0], P[1], P[3], P[2], P[4], P[5], P[7], P[6]]

def vtk2gmsh_hex(P):
    return [P[0], P[3], P[2], P[1], P[4], P[5], P[6], P[7]]
	
def ufl_simplicial_order(P):
    return np.sort(np.array(P), axis=None).tolist()	 # sorts vector P in ascending order

# loadPath = "./"
# loadPath = "./sq-8mm_sin-per-4_sin-amp-2mm_tet/"  
# savePath = "."
# filename = "sq-8mm_sin-per-4_sin-amp-2mm_tet"

filenames = ["sq-8mm_sin-per-2_sin-amp-1mm_tet","sq-8mm_sin-per-2_sin-amp-2mm_tet","sq-8mm_sin-per-2_sin-amp-3mm_tet","sq-8mm_sin-per-4_sin-amp-1mm_tet","sq-8mm_sin-per-4_sin-amp-2mm_tet","sq-8mm_sin-per-4_sin-amp-3mm_tet"]
loadPath = "/home/fenics/shared/ogdenShearVSI-data/sensitivity_data/"
savePath = "/home/fenics/shared/ogdenShearVSI-data/sensitivity_data/"

# Convert quad tet mesh (.inp) to a linear tet python mesh (.xdmf/.h5)
def inp2xdmf(loadPath,savePath,filename):

    # Write quad tet data to meshio.Mesh object
    mesh_matlab = meshio.read(os.path.join(loadPath, filename + '.inp'))

    # Find quad tet mesh size params
    num_vertex = mesh_matlab.points.shape[0] 
    num_cell = mesh_matlab.cells[0].data.shape[0]

    # Initialize FEniCS mesh 
    mesh_py = Mesh() # create empty dofin.Mesh object 
    editor = MeshEditor() # iniate MeshEditor object
    editor.open(mesh_py, type="tetrahedron",tdim=3, gdim=3) # define tetrahedron (or "hexahedron") mesh with 3 topologic and geometric dimensions
    # editor.init_vertices(num_vertex) # initialize cells in dolfin mesh - the number of verticies is equal to the number of linear nodes on the quad tet mesh 
    editor.init_vertices(np.max(mesh_matlab.cells[0].data[:,:4])+1) # initialize cells in dolfin mesh - the number of verticies is equal to the number of linear nodes on the quad tet mesh 
    editor.init_cells(num_cell) # initialize cells in python mesh - the number of elements is equal to the quad tet mesh

    # for i in range(num_vertex): 
    for i in range(np.max(mesh_matlab.cells[0].data[:,:4])+1): 
        editor.add_vertex(i, Point(mesh_matlab.points[i,:])) # defining coordinates of ith vertex from the nodes included on the linear portion of the quadratic tets
    for i in range(num_cell):
        try:		
            _cell_dolfin = ufl_simplicial_order(mesh_matlab.cells[0].data[i,:]) # defining verticies of the ith cell in ascending order 
            editor.add_cell(i,_cell_dolfin[0:4]) # add first 4 nodes (linear portion of tet) to ith connectivity matrix  to MeshEditor object. Effectively converts quad --> linear tet mesh.
        except RuntimeError:
            print("Error in cell index %i"%i)
            print(mesh_matlab.cells[0].data[i,:4])
            raise
    # Call out relevant information about the transitioin to linear tets from quadratic tets
    print("linear tet uses nodes between",np.min(mesh_matlab.cells[0].data[:,:4]),"to",np.max(mesh_matlab.cells[0].data[:,:4])) 
    print("quad tet uses nodes between",np.min(mesh_matlab.cells[0].data[:,:]),"to",np.max(mesh_matlab.cells[0].data[:,:]))
    print("quad tet uses linear nodes and",np.min(mesh_matlab.cells[0].data[:,4:]),"to",np.max(mesh_matlab.cells[0].data[:,4:]))
    print("num_vertex variable is equal to",num_vertex)
    try: 
        print('==============================')
        editor.close() # close MeshEditor object
        print("Editor closed")
    except RuntimeError as err:
        print("Error in closing the editor")
        print(err)
        pass	
    try:
        print('==============================')
        mesh_py.order() # Finding the highest order of basis function used in the mesh
        print('Mesh ordered: ', mesh_py.ordered()) # checking that the mesh is well-defined in terms of the ordering of its verticies and cells
    except RuntimeError:
        print("Error ordering the mesh and generating mesh function")
        raise
    # Save python mesh to .xdmf/.h5 files
    try:
        print('==============================')
        with XDMFFile(MPI.comm_world,os.path.join(savePath,filename+'.xdmf')) as file:
            file.write(mesh_py)
        print("Mesh saved")
    except RuntimeError:
        print("Error in saving file")
        raise

    # Attempt to open mesh
    mesh_py = Mesh()

    # Final check on .xdmf form of Python mesh
    try:
        print('==============================')
        with XDMFFile(MPI.comm_world,os.path.join(savePath,filename+'.xdmf')) as file:
            file.read(mesh_py)
            V = VectorFunctionSpace(mesh_py, "CG", 1)
            a = dof_to_vertex_map(V)
            np.savetxt(os.path.join(savePath, 'dof2vertex.txt'),a, fmt='%i') # save DOF to vertex map for linear tet mesh
            print('Mesh is readable')	
    except RuntimeError as err:
        print("Error in reading the file")
        print(err)

# LOOPING CODE: comment out regular inp2xdmf function call at the end of the code and uncomment the code below to run the script in a loop 
# filenames = ["sq-8mm_sin-per-2_sin-amp-1mm_tet","sq-8mm_sin-per-2_sin-amp-2mm_tet","sq-8mm_sin-per-2_sin-amp-3mm_tet","sq-8mm_sin-per-4_sin-amp-1mm_tet","sq-8mm_sin-per-4_sin-amp-2mm_tet","sq-8mm_sin-per-4_sin-amp-3mm_tet"]
# for i in range(len(filenames)):
#     print("Running file",i+1,"of",len(filenames))
#     loadPath_loop = "/home/fenics/shared/ogdenShearVSI-data/sensitivity_data/" + str(filenames[i]) + "/"
#     savePath_loop = loadPath_loop + "VSI/"
#     inp2xdmf(loadPath_loop,savePath_loop,str(filenames[i]))

# SINGLE CODE: comment out if running batch code
# inp2xdmf(loadPath,savePath,filename)

# BATCH CODE:
for i in range(len(filenames)):
    inp2xdmf(loadPath+filenames[i]+'/',loadPath+filenames[i]+'/VSI/',filenames[i])