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

# filename = "STA26_27"
# loadPath = "./InputFiles/22-1215-Wavy_Sweep_v2/"
loadPath = "./"
savePath = "."
filename = "ShearWavy_6.25MMDisp_Amp_0.6_tet"

# Convert MATLAB mesh (.inp) to a python mesh (.xdmf/.h5)
def msh2xdmf(loadPath,savePath,filename):

    msh_matlab = meshio.read(os.path.join(loadPath, filename + '.inp')) # read .inp mesh file from MATLAB and create meshio.Mesh object to store information describing the mesh

    #Find mesh size params
    num_vertex = msh_matlab.points.shape[0] 
    num_cell = msh_matlab.cells[0].data.shape[0]

    #Initialize FEniCS mesh 
    mesh_py = Mesh() # create empty dofin.Mesh object 

    editor = MeshEditor() # iniate MeshEditor object
    editor.open(mesh_py, type="tetrahedron",tdim=3, gdim=3) # define tetrahedron (or "hexahedron") mesh with 3 topologic and geometric dimensions
    editor.init_vertices(num_vertex) # defining verticies
    editor.init_cells(num_cell) # initialize cells

    for i in range(num_vertex):
        editor.add_vertex(i, Point(msh_matlab.points[i,:])) # defining coordinates of ith vertex
    for i in range(num_cell):
        try:		
            _cell_dolfin = ufl_simplicial_order(msh_matlab.cells[0].data[i,:]) # defining verticies of the ith cell in ascending order 
            editor.add_cell(i,_cell_dolfin) # add ith cell to MeshEditor object 
        except RuntimeError:
            print("Error in cell index %i"%i)
            print(msh_matlab.cells[0].data[i,:4])
            raise
    try: 
        print('==============================')
        editor.close()
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
    
    # Saving python mesh to .xdmf file
    try:
        print('==============================')
        with XDMFFile(MPI.comm_world,os.path.join(savePath,filename+'.xdmf')) as file:
            file.write(mesh_py)
        print("Mesh saved")
    except RuntimeError:
        print("Error in saving file")
        raise

    # Now try opening the file
    mesh_py = Mesh()

    try:
        print('==============================')
        with XDMFFile(MPI.comm_world,filename+'.xdmf') as file:
            file.read(mesh_py)
            V = VectorFunctionSpace(mesh_py, "CG", 1)
            a = dof_to_vertex_map(V)
            # np.savetxt('dof2vertex.txt',a, fmt='%i') # check DOF to vertex map
            print('Mesh is readable')	
    except RuntimeError as err:
        print("Error in reading the file")
        print(err)
msh2xdmf(loadPath,savePath,filename)