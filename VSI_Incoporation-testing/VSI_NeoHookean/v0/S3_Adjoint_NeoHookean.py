from ufl import *
from dolfin import *
import numpy as np
#import h5py as h5
from dolfin_adjoint import *
set_log_active(True)
import os
import sys
import time as timelib

#a = Constant(1.1)
#b = Control(a)
#print(b)
#print(ass)
#exit()

def iter_cb(m):
  print ("m = ", m)

def eval_cb(j, m):
  print ("j = %f, m = %f." % (j, m))

def derivative_cb(j, dj, m):
  print("j = %f, dj = %s, m = %s." %(j, [float(k) for k in dj], [float(k) for k in m]))


################################ Main Code ####################################
###############################################################################


################################ Problem Setup ################################




## Create Mesh and Define Function Space
#datadir='../result/unfiltered/20220211-holes/raw_data'	
result_dir = './'
in_mesh="./sq-8mm_sin-per-4_sin-amp-2mm_tet.xdmf"
	
data_tag = 'specimen_0'	
mesh_tag = 'test0'
plot_dir = result_dir + 'plots/'
solution_tag = 'adj_disp'

data_dir = result_dir

#In the data, force is in N. 
#Units of stress is taken to be N/mm^2
#Denis's inital values are in kPa. x kPa = x * 1e-3 N/mm^2

# Initial values for 2 parameter model
bounds=np.zeros((2,2))
bounds[0,0]= 1e-4	#Positive E
bounds[1,0]= np.inf	
bounds[0,1]= -1. #Range for poisson ratio = [-1,0.5]
bounds[1,1]= 0.5 

#set initial values and control

theta = [Constant(6.0e-1), Constant(0.11)]  # [E(MPa), nu(dimensionless)]
control_index = [0,1]

# # Initial values for parameter for model version 2
# bounds=np.zeros((2,3))
# bounds[0,:]= 0	#Positive E
# bounds[1,:]= np.inf	

# #set initial values and control

# theta = [Constant(1e3), Constant(2e3), Constant(3e3)]  # [E(MPa), nu(dimensionless)]
# control_index = [0,1,2]


# control_index=[]
# theta=[]
# for i in range(len(gamma)):
#   theta.append(Constant(gamma[i]))  
#   control_index.extend([i])			#To make sure only non-zero gamma are fine tuned


#Set forward

start_forward = timelib.perf_counter()	

saveFig_flag = True
NumStep = 2 #Corresponds to time-step in Step2_addU2xdmf.py
adjoint_dict = {'estimate_loss_flag':True,
'NumStep': NumStep,  
'data_dir':result_dir,
'data_tag':data_tag, 
'force_val':38.4139,
'dispU_val':3.0,
}


#u, strain, stress, loss = forward_model(in_mesh, result_dir, mesh_tag, solution_tag, saveFig_flag, adjoint_dict, theta)

#===========================================================================================
#============================ Forward model starts here ====================================
#===========================================================================================
q_degree = 3
dx_ = dx(metadata={'quadrature_degree': q_degree})
zeros = Constant((0.0, 0.0, 0.0))
tol = 1E-14

mesh=Mesh()

print(in_mesh)

with XDMFFile(MPI.comm_world,in_mesh) as infile:
	infile.read(mesh)


#mesh=Mesh()
#with XDMFFile(result_dir + mesh_tag + '.h5') as infile:
#	infile.read(mesh)

V_vector = VectorFunctionSpace(mesh, 'Lagrange', 1)
u = Function(V_vector, name='u')
x=SpatialCoordinate(mesh)
dof_coordinates = V_vector.tabulate_dof_coordinates()                    
dof_coordinates.resize((V_vector.dim(), mesh.geometry().dim()))                           
dof_x = dof_coordinates[:, 0] 
dof_y = dof_coordinates[:, 1]
dof_z = dof_coordinates[:, 2]	

#x-(up direction), y-(out of page), z-(right direction) 
# bottom =  CompiledSubDomain("near(x[1] , (-4 + 2*sin(2*pi*x[0]/5))) && on_boundary")
# top =  CompiledSubDomain("near(x[1] , (4 - 2*sin(2*pi*x[0]/5))) && on_boundary")
bottom =  CompiledSubDomain("x[1]< (-4 + 2*sin(pi*x[0]/5)) + 1e-1 && on_boundary")
top =  CompiledSubDomain("x[1]> (4 - 2*sin(pi*x[0]/5)) - 1e-1 && on_boundary")


# def bottom(x, on_boundary):
#     return on_boundary and x[1]  < (-4 + 2*sin(2*pi*x[0]/5)) + _tol
# def top(x, on_boundary):
#     return on_boundary and x[1]  > (4 - 2*sin(2*pi*x[0]/5)) -  _tol



back =  CompiledSubDomain("near(x[2], side) && on_boundary", side = dof_z.min())
front =  CompiledSubDomain("near(x[2], side) && on_boundary", side = dof_z.max())
left =  CompiledSubDomain("near(x[0], side) && on_boundary", side = dof_x.min())
right =  CompiledSubDomain("near(x[0], side) && on_boundary", side = dof_x.max())	

#Top surface normal
#[cos(atan(dy/dx)), sin(atan(dy/dx)), 0]
#top_surf_normal = Expression(("1/(pow(1+pow(pi*cos(pi*x[0]/5)/5,2),0.5))","x[0]/(pow(1+pow(pi*cos(pi*x[0]/5)/5,2),0.5))","0"), degree=2) 
n = FacetNormal(mesh)
	
#print(dof_x.min())
#print(dof_x.max())
#exit()
#surface =  CompiledSubDomain('on_boundary')
#bcl_V = DirichletBC(V_vector, up, facetfct,0)
#bcl_S = DirichletBC(V_vector, zeros, surface)

#Set dispalcements here
top_disp = Constant((adjoint_dict['dispU_val'],0.,0.))
bottom_disp = Constant((0.0, 0.0, 0.0))



#bcl_bottom = DirichletBC(V_vector.sub(2), Constant(0.0), left)
#bcl_top = DirichletBC(V_vector, Constant((0.0, 0.0, 6.0)), right)	
bcl_top = DirichletBC(V_vector, top_disp, top)
bcl_bottom = DirichletBC(V_vector, bottom_disp, bottom)

bcs = [bcl_top, bcl_bottom]

v = TestFunction(V_vector) 

## Elasticity Parameters
d = len(u)
Zero = zero((d,d))             # Identity Tensor
I = Identity(d)             # Identity Tensor
F = I + grad(u)             # Deformation Gradient
C = F.T*F                   # Right Cauchy-Green Tensor

# #===============================Neo Hookean model begins here===============================

## Elasticity Parameters
#E, nu = 5e3, 0.49
mu, lmbda = theta[0]/(2*(1 + theta[1])), theta[0]*theta[1]/((1 + theta[1])*(1 - 2*theta[1]))

## Invariants of Deformation Tensors
Ic = tr(C)
Jc  = det(C)
invC = inv(C)

## Compressible Neo-Hookean Model
psi = 0.25*lmbda*(Jc-1)-0.5*(lmbda/2+mu)*ln(Jc)+0.5*mu*(inner(F,F)-3)
Pi = psi*dx_

S = 0.5*lmbda*Jc*invC-(0.5*lmbda+mu)*invC+mu*I
P = F*S


# #===============================NeoHookean model ends here ===============================

# #=============================== Test Neo-Hookean version 2 ===============================

# Ic = tr(C)
# Jc  = det(C)
# invC=inv(C)
# psi=0.25*theta[0]*(Jc-1)-0.5*theta[1]*ln(Jc)+0.5*theta[2]*(inner(F,F)-3)
# Pi = psi*dx 

# S=0.5*theta[0]*Jc*invC -(0.5*theta[1])*invC + theta[2]*I
# P=F*S
# #=============================== Test Neo-Hookean ends ===============================


file = HDF5File(MPI.comm_world, result_dir + mesh_tag + '.h5', 'w')
file.write(mesh, '/mesh')

R = inner(P, grad(v))*dx_
J = derivative(R, u)

problem = NonlinearVariationalProblem(R, u, bcs, J)
solver = NonlinearVariationalSolver(problem)

prm = solver.parameters
prm['newton_solver']['absolute_tolerance'] = 1E-8
prm['newton_solver']['relative_tolerance'] = 1E-9
prm['newton_solver']['maximum_iterations'] = 25
prm['newton_solver']["error_on_nonconvergence"] = True
# prm['nonlinear_solver']='snes'
# prm['snes_solver']['line_search']='bt'
# prm['snes_solver']['linear_solver']='lu'
# prm['snes_solver']['maximum_iterations'] = 100
# prm['snes_solver']["error_on_nonconvergence"] = True
solver.solve()

if saveFig_flag:
	print('forward: total energy E =', assemble(Pi))

	file.write(u,'/displacement')

	## Save to .xdmf
	file = XDMFFile(mesh.mpi_comm(),result_dir + solution_tag + '_plot.xdmf')
	file.write(u)
	file.close()

loss = 0
if adjoint_dict['estimate_loss_flag']: 
	print('Reading data for NumStep %i'%adjoint_dict['NumStep'])
	u_true = Function(V_vector, name='u_true')
	with HDF5File(MPI.comm_world, adjoint_dict['data_dir']+adjoint_dict['data_tag'] + '.h5','r') as infile:
		infile.read(u_true, 'disp%i'%adjoint_dict['NumStep'])	

	data_Term = (inner(u - u_true, u - u_true)) *dx
	#Estimate net traction on top surface
	
	traction = dot(P, n)
	#Traction term on front surface
	boundary_subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
	boundary_subdomains.set_all(0)
	top.mark(boundary_subdomains,1)
	ds_top = Measure("ds", domain=mesh, subdomain_data=boundary_subdomains, subdomain_id=1, metadata={'quadrature_degree': q_degree})
	area_top = assemble(1*ds_top)
	
	x_dir = Constant((1.0, 0.0, 0.))
	y_dir = Constant((0.0, 1.0, 0.))
	z_dir = Constant((0.0, 0.0, 1.))
	load_Term = dot(traction , x_dir) *ds_top

	#print(area_top)
	#print(assemble( load_Term ))
	
	loss1 = assemble( data_Term )
	loss2 = (assemble( load_Term ) - adjoint_dict['force_val'])**2
	
	print('Displacement Loss: %f'%loss1)
	print('Load Loss        : %f'%loss2)
	loss_factor = 1e0
	loss = loss1+loss_factor*loss2

#=========================================================================================
#============================ Forward model ends here ====================================
#=========================================================================================

end_forward = timelib.perf_counter()	

#do adjoint
control_parameter=[Control(theta[i]) for i in control_index]  

	


start_adjoint = timelib.perf_counter()	
maxiter=50
method = "L-BFGS-B"#"SLSQP"#"L-BFGS-B"
print('Estimating reduced functional')
reduced_functional = ReducedFunctional(loss, control_parameter, derivative_cb_post=derivative_cb)
print('Starting minimization process')
results_opt = minimize(reduced_functional, method = method, bounds=bounds, tol=1.0e-9, options = {'ftol':1.0e-11,'maxiter':maxiter, 'disp': True},callback = iter_cb)
#results_opt = minimize(reduced_functional, method = method, bounds=bounds,tol=1.0e-9, options = {'ftol':1.0e-11,'maxiter':maxiter, 'disp': True},callback = iter_cb)
#converge_flag=True 


end_adjoint = timelib.perf_counter()	


gamma_matrix_all = results_opt#[tem for i in range(total_num_step)]



np.savetxt(result_dir+ 'Adjoint_sol.txt',gamma_matrix_all)

print('Forward solve time: %f seconds'%(float(end_forward)-float(start_forward)))
print('Adjoint solve time: %f seconds'%(float(end_adjoint)-float(start_adjoint)))

