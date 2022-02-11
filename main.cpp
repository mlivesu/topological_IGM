#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/geometry/n_sided_poygon.h>
#include <cinolib/cut_mesh.h>


using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void IGM(const uint genus, Quadmesh<> & igm)
{
    std::vector<vec3d> verts = n_sided_polygon(4*genus, CIRCLE);
    for(const auto & p : verts) igm.vert_add(p);
    uint origin = igm.vert_add(vec3d(0,0,0)); // center point

    for(uint i=0; i<genus; ++i)
    {
        // 5 consecutive corners of the Canonical Polygonal Schema
        uint cps[5] =
        {
            4*i,
            4*i+1,
            4*i+2,
            4*i+3,
           (4*i+4)%(genus*4)
        };

        // auxiliary inner vertices
        uint aux[3] =
        {
            igm.vert_add((verts[cps[1]] + verts[cps[2]] + igm.vert(origin))/3),
            igm.vert_add(0.7*verts[cps[2]] + 0.3*igm.vert(origin)),
            igm.vert_add((verts[cps[2]] + verts[cps[3]] + igm.vert(origin))/3)
        };

        igm.poly_add({cps[0], cps[1], aux[0], origin});
        igm.poly_add({cps[1], cps[2], aux[1], aux[0]});
        igm.poly_add({cps[2], cps[3], aux[2], aux[1]});
        igm.poly_add({cps[3], cps[4], origin, aux[2]});
        igm.poly_add({aux[0], aux[1], aux[2], origin});
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main()
{
    GLcanvas gui;
    DrawableQuadmesh<> m;
    IGM(3,m);
    m.updateGL();

    // texturing
    for(uint eid=0; eid<m.num_edges(); ++eid)
    {
        m.edge_data(eid).flags[MARKED] = m.edge_valence(eid)>1;
    }
    cut_mesh_along_marked_edges(m);
    for(uint pid=0; pid<m.num_polys(); ++pid)
    {
        m.vert_data(m.poly_vert_id(pid,0)).uvw = vec3d(0,0,0);
        m.vert_data(m.poly_vert_id(pid,1)).uvw = vec3d(1,0,0);
        m.vert_data(m.poly_vert_id(pid,2)).uvw = vec3d(1,1,0);
        m.vert_data(m.poly_vert_id(pid,3)).uvw = vec3d(0,1,0);
    }
    m.show_texture2D(TEXTURE_2D_CHECKERBOARD, 1);
    m.updateGL();

    SurfaceMeshControls<DrawableQuadmesh<>> controls(&m,&gui,"IGM");
    gui.push(&m);
    gui.push(&controls);
    return gui.launch();
}
