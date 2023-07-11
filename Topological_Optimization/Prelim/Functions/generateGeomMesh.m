function f = generateGeomMesh(x)
global glb_var
mesh_ref = glb_var.mesh_ref;
mesh_ref.maxelsize = x(1);
l = glb_var.l; w = glb_var.w; h = glb_var.h; line_res = glb_var.line_res; edge = glb_var.edge;
model_3D = createGeometry(l,w,h,line_res,edge,mesh_ref);

f = abs(mesh_ref.num_of_el-size(model_3D.Mesh.Elements,2));

end