function [mesh_ref] = optimize_mesh(model_3D,mesh_ref,l,w,h,line_res,edge)
global glb_var
glb_var.mesh_ref = mesh_ref;
glb_var.l = l; glb_var.w = w; glb_var.h = h; glb_var.line_res = line_res; glb_var.edge = edge;
iter2_x0 = model_3D.Mesh.MaxElementSize + (size(model_3D.Mesh.Elements,2)-mesh_ref.num_of_el)*model_3D.Mesh.MaxElementSize/mesh_ref.num_of_el;

options_unc = optimoptions(@fminunc,'Display','iter');
% options_con = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter');
tic
xsol = fminunc(@ (x) generateGeomMesh(x),iter2_x0,options_unc);
% xsol = fmincon(@ (x) generateGeomMesh(x),iter2_x0,[],[],[],[],0,2*iter2_x0,[],options_con);
toc
mesh_ref.maxelsize = xsol;

clear glb_var
close(gcf)
end