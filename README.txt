Dataset Title: MRI Raw Experimental Data Pipeline (as of 10/06/2022)
Dataset Creators: D.P. Nikolov, S. Srivastava, B.A. Abeid, U.M. Scheven, E.M. Arruda, K. Garikipati, J.B. Estrada
Dataset Contact: J.B. Estrada jbestrad@umich.edu
Funding: 1729166 (NSF)
-------------------------------------------
Research Overview:
Contemporary material characterisation techniques that leverage deformation fields and the weak form of the equilibrium equations face challenges in the numerical solution procedure of the inverse characterisation problem.  As material models and descriptions differ, so too must the approaches for identifying parameters and their corresponding mechanisms. The widely-used Ogden material model can be comprised of a chosen number of terms of the same mathematical form, which presents challenges of parsimonious representation, interpretability, and stability.  Robust techniques for system identification of any material model are important to assess and improve experimental design, in addition to their centrality to forward computations. Using fully 3D displacement fields acquired in silicone elastomers with our recently-developed magnetic resonance cartography (MR-u) technique on the order of ~20,000 points per sample, we leverage PDE-constrained optimisation as the basis of variational system identification of our material parameters. We incorporate the statistical F-test to maintain parsimony of representation.  Using a new, local deformation decomposition locally into mixtures of biaxial and uniaxial tensile states, we evaluate experiments based on an analytical sensitivity metric, and discuss the implications for experimental design.  

Methodology:
This repository contains the acquired kinematic data and MRI processing code used in this work.

Instruments:
	- 7T small-animal MRI system
	- Distant captive linear actuator (L5918S2008-T10X2-A50, Nanotec Electronic GmbH and Co. KG, Germany)
	- Load cell (LCM300, Futek Advanced Sensor Technology Inc., Irvine, CA)
	- Silicone samples (Ecoflex OO-20 formulation; Dragon Skin, Smooth-On Inc., Macungie, PA)

Software: MathWorks MATLAB v. 2020b or later (Natick, MA); Abaqus FEA (Providence, RI)
-------------------------------------------
Notes: 		Archive has all datasets and files from prior research and aren't intended for the main investigation of this paper
      		It's intended for future investigations of the authors of this paper, however the resources are publically available for miscellaneous investigation
		The primary files for this paper are in the following parent directories located in the main directory with this readme file:
			- UMItools
			- Ogden_RawMRData
			- MR_Processing
			- Abaqus2Matlab
			- Decomposition_Sensitivity
			- Visualization

File Inventory:	For further explanation of the datasets in each of these folders, refer to Data_Inventory.xlsx
		For each run/main function of the dataset, you may utilize the following function to see dependencies and all other functions used in the primary functions:
			- [fList,pList] = matlab.codetools.requiredFilesAndProducts('[functions_name].m');
			- fList contains all the dependent functions
-------------------------------------------
Use and Access:

1. Open MATLAB and add all folders and subfolders to path
2. Processing raw data
	- Open Ogden_RawMRData directory
	- See table below to enter the directory for a given sample (i.e. for Solid_Rectangular (7 mm), enter 20211012-ogdenss-0020/silicone0020)
	- Run the following inline function:
		og=ogden_deste_3d(PE,SL,RO,Ref,'blurvec',[0.8,0.8,0.8]);
	- A file named 'DESTE_strains_RO_PE_SL_Ref_blurvec0.80.80.8.mat' is produced in directory that contains the og structure with the required data for processing

	Exp Date	Sample Type		Prescribed displacement	     Timesteps:	PE    SL    RO    Ref
	------------------------------------------------------------------------------------------------------
	211012		Solid_Rectangular 	7 mm					1958, 1904, 1931, 1836
	220124		Holes_Rectangle 	5 mm					0924, 0950, 1016, 0852
	220124		Holes_Rectangle 	2.5 mm					1139, 1204, 1113, 0852
	220209		Solid_Rectangular 	5 mm					1336, 1402, 1429, 1554
	220209		Solid_Rectangular 	2.5 mm					1456, 1521, 1546, 1554

3. Numerical differentiation for kinematic quantities
	- Open MR_Processing directory
	- Open runComplexFilter_Ogden.m
	- Uncomment the following for 'sampleName' variable to run for a particular sample:
		20211012-ogdenss			Solid_Rectangular	7 mm
		20220209-ogdenss_5MMH_apod_64_16	Holes_Rectangular	2.5 mm & 5 mm
		20220209-ogdenss_2moreloadsteps		Solid_Rectangular	2.5 mm & 5 mm
	- Due to size limitations, you'll have to visit https://drive.google.com/file/d/1sl29qNNzVtvLLefp_xoE7m4ru78bWVhV/view?usp=sharing and download the zip file.
	- Populate the DESTE_strains_RO_PE_SL_Ref_blurvec0.80.80.8.mat files into their respective experiment subdirectory (labeling in the zip file should help out)
	- Re-enter the parent directory of MR_Processing (or you may comment out 'cd(sampleName)' in the code)
	- Run runComplexFilter_Ogden
	- A file named 'disp_data_sampleName.mat' is produced in the respective sampleName subdirectory that contains the following information:
		F_t{t}{i,j}(x,y,z) 	- Deformation gradient tensor
					- t: displacement step (1 - 2.5 mm, 2 - 5 mm, 3 - 7 mm)
					- i,j: component of deformation gradient
					- x,y,z: voxel mesh grid indices
		E_t{t}{i,j}(x,y,z)	- Lagrange strain tensor
		U_t{t}{i}(x,y,z)	- Displacement vector
		mask(x,y,z)		- Mask of material
		prescribedU		- prescribedU values
		osc			- voxel resolution (mm/px)
	- A file named 'refpositions_sampleName.mat' is produced in the respective sampleName subdirectory that contains the following information:
		X{j}(x,y,z)		- Reference configuration positions
					- j: 1,2 or 3 coordinate axis direction

4. Simulation deformation gradients:
	- Open Abaqus2Matlab directory
	- Open DefGrad_sim_v8
	- Uncomment the following for 'curdir' variable to run for a particular sample:
		22-0201-5MM_Holes_SS	Holes_Rectangular	2.5 mm & 5 mm
		22-0301-NoHoles_SS	Solid_Rectangular	2.5 mm & 5 mm & 7 mm
		22-0325-Uniaxial	Uniaxial		7 mm
	- A file named 'MRI-3Ddefs_SimpleShear_' curdir '.mat' is produced in the respective 'Simulations_tet\curdir' subdirectory and contains the deformation information on the deformation gradient and displacement of each element
	- A file named 'refpositions.mat' is produced in the respective 'Simulations_tet\curdir' subdirectory and contains information about the referenece configuration positions
	- Note: Abaqus simulations don't run by default, since the 'Data' files are already produced. If you'd like to inspect how Abaqus runs in MatLab, simply delete the folder with the respective data in the 'Data' subdirectory
	- Feel free to use your own input files from Abaqus to investigate DefGrad_sim_v8. Some guidelines to follow:
		Create a new case for 'curdir' and change the size/hole inclusion in genMesh within the subdirectory 'Functions\'
	- Note: This code runs for delaunay triangulation within MatLab and hence calculates the deformation gradients for quadratic tetrahedral elements
		+ For hexahedral elements, you'll need to use the files that are in the folders 'Old_hex' to run them
	- If you run into any unforeseen issues, feel free to contact Denislav Nikolov: dnikolov@umich.edu

5. Decoupling parameters and sensitivity metric quantification
	- Open Decomposition_Sensitivity directory
	- k and lambda decoupling of F_t:
		+ Load the MRI-3Ddefs_SimpleShear...mat file from the respective curdir subdirectory into MATLAB workspace
		+ Load og_matprop.mat
		+ Run one of the following inline function (depending on whether you'd like to use parfor or par):
			[k,lam,~,~,~] = param_decoup_nopar(F_t,og_matprop);
			[k,lam,~,~,~] = param_decoup_main(F_t,og_matprop);
			'k' and 'lam' contain information for the k and lambda decoupling for each voxel's deformation gradient
	- Sensitivity plots:
		+ Open 'Plot2DHists_Both.m' (You can change whether you'd want filtered data or none and see the effect of filtering helping significantly with the plots)
		+ Run 'Plot2DHists_Both.m' for the plots in figure 7a and 7b
		+ The figures and data are saved in the Filter or Filter_None subdirectories
		+ Open 'Plot2DHists_Uniaxial_Both_mu_alpha.m' (You can change whether you'd want filtered data or none and see the effect of filtering helping significantly with the plots)
		+ Run 'Plot2DHists_Uniaxial_Both_mu_alpha.m' for the plot in figure 7c (Ensure to run Plot2DHists_Both first to populate the appropriate data)
		+ The figure is save in the directory as 'sens_metric_all_v2.png' and 'sens_metric_all_v2.pdf'

6. Visualization of data (Reproduce processing pipeline, and displacement plots figure)
	- Open Visualization
	- Run 'run_plots_complex.m' for all the processing pipeline images (images saved/located in 5a-Complex folder)
	- Run 'run_plots_displacement.m' for all the displacement field images (images saved/located in 5bcd_V2-Slice_U)
