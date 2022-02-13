#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/geometry/n_sided_poygon.h>
#include <cinolib/cut_mesh.h>
#include <cinolib/canonical_polygonal_schema.h>
#include <cinolib/homotopy_basis.h>
#include <cinolib/octree.h>
#include <cinolib/geometry/quad_utils.h>
#include <cinolib/geometry/segment_utils.h>

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void map_igm_cps(Quadmesh<> & igm,
                 Trimesh<>  & cps)
{
    Octree igm_edges;
    igm_edges.build_from_mesh_edges(igm);

    // edge IDs are not safe if edge_split is called
    // operate on a queue of pairs of vertex ids
    std::queue<vec2i> q;
    for(uint eid=0; eid<cps.num_edges(); ++eid)
    {
        if(!cps.edge_is_boundary(eid))
        {
            q.push(vec2i(cps.edge_vert_id(eid,0),
                         cps.edge_vert_id(eid,1)));
        }
    }
    while(!q.empty())
    {
        vec2i vids = q.front();
        vec3d s[2] = {cps.vert(vids[0]),
                      cps.vert(vids[1])};
        q.pop();
        std::unordered_set<uint> inters;
        if(igm_edges.intersects_segment(s, true, inters))
        {
            // TODO: handle multiple intersections!
            vec3d p = segment_intersection(cps.vert(vids[0]),
                                           cps.vert(vids[1]),
                                           igm.edge_vert(*inters.begin(),0),
                                           igm.edge_vert(*inters.begin(),1));

            int eid = cps.edge_id(vids[0],vids[1]);
            assert(eid>=0);
            cps.edge_split(eid,p);
        }
    }

    /* TODO:
     * insert IGM inner points
     * map IGM vertices to CPS vertices
     * mark IGM paths on CPS
     * cut mesh along such paths
     * texturing (done below)
     */
}

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
    DrawableTrimesh<> m("/Users/cino/Desktop/tmp_data/eight.off");
    uint genus = m.genus();

    // compute optimal homotopy basis
    HomotopyBasisData hb;
    hb.globally_shortest = false;
    hb.detach_loops      = true;
    hb.split_strategy    = EDGE_SPLIT_STRATEGY;
    homotopy_basis(m,hb);

    DrawableTrimesh<> cps;
    canonical_polygonal_schema(m, hb, cps);

    DrawableQuadmesh<> igm;
    IGM(genus,igm);
    igm.show_marked_edge_width(5);
    igm.show_mesh_points();
    igm.updateGL();

    map_igm_cps(igm,cps);
    cps.updateGL();

    Octree o;
    o.build_from_mesh_polys(igm);
    for(uint vid=0; vid<cps.num_verts(); ++vid)
    {
        uint   pid;
        vec3d  pos;
        double dist;
        o.closest_point(cps.vert(vid), pid, pos, dist);
        assert(dist<1e-10);

        vec4d bary;
        quad_barycentric_coords(igm.poly_vert(pid,0),
                                igm.poly_vert(pid,1),
                                igm.poly_vert(pid,2),
                                igm.poly_vert(pid,3),
                                cps.vert(vid), bary);

        cps.vert_data(vid).uvw = vec3d(0,0,0) * bary[0] +
                                 vec3d(1,0,0) * bary[1] +
                                 vec3d(1,1,0) * bary[2] +
                                 vec3d(0,1,0) * bary[3];
    }

    GLcanvas gui;
    gui.push(&cps);
    gui.push(&igm);
    gui.push(new SurfaceMeshControls<DrawableTrimesh<>> (&cps,&gui,"CPS"));
    gui.push(new SurfaceMeshControls<DrawableQuadmesh<>>(&igm,&gui,"IGM"));
    return gui.launch();
}
