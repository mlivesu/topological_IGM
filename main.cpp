#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/geometry/n_sided_poygon.h>

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main()
{
    uint genus = 2;
    std::vector<vec3d> verts = n_sided_polygon(4*genus, CIRCLE);
    DrawablePolygonmesh<> m(verts);
    GLcanvas gui;
    SurfaceMeshControls<DrawablePolygonmesh<>> controls(&m,&gui,"Polygonal Schema");
    gui.push(&m);
    gui.push(&controls);
    return gui.launch();
}
