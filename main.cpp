#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/geometry/n_sided_poygon.h>

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void IGM(const uint genus, DrawableQuadmesh<> & igm)
{
    uint c = igm.vert_add(vec3d(0,0,0));
    std::vector<vec3d> verts = n_sided_polygon(4*genus, CIRCLE);
    for(const auto & p : verts) igm.vert_add(p);

    for(uint i=0; i<genus; ++i)
    {
        uint off = 4*i;
        uint v0 = igm.vert_add((verts[off+1] + igm.vert(c) + verts[off+2])/3);
        uint v1 = igm.vert_add(0.70*verts[off+2] + 0.30*igm.vert(c));
        uint v2 = igm.vert_add((verts[off+2] + igm.vert(c) + verts[off+3])/3);

        igm.poly_add({off+1, off+2, v0,  c});
        igm.poly_add({off+2, off+3, v1, v0});
        igm.poly_add({off+3, off+4, v2, v1});
        igm.poly_add({off+4, (off+4)%igm.num_verts()+1,  c, v2});
        igm.poly_add({   v0,    v1, v2,  c});
    }

    igm.updateGL();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main()
{
    DrawableQuadmesh<> m;
    IGM(3,m);
    GLcanvas gui;
    SurfaceMeshControls<DrawableQuadmesh<>> controls(&m,&gui,"IGM");
    gui.push(&m);
    gui.push(&controls);
    return gui.launch();
}
