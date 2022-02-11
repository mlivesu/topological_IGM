#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/geometry/n_sided_poygon.h>

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void IGM(const uint genus, DrawableQuadmesh<> & igm)
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
        uint inn[3] =
        {
            igm.vert_add((verts[cps[1]] + verts[cps[2]] + igm.vert(origin))/3),
            igm.vert_add(0.7*verts[cps[2]] + 0.3*igm.vert(origin)),
            igm.vert_add((verts[cps[2]] + verts[cps[3]] + igm.vert(origin))/3)
        };

        igm.poly_add({cps[0], cps[1], inn[0], origin});
        igm.poly_add({cps[1], cps[2], inn[1], inn[0]});
        igm.poly_add({cps[2], cps[3], inn[2], inn[1]});
        igm.poly_add({cps[3], cps[4], origin, inn[2]});
        igm.poly_add({inn[0], inn[1], inn[2], origin});
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
