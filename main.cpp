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

    Octree o;
    o.build_from_mesh_edges(cps);
    cps.edge_set_flag(MARKED, false);
    std::vector<vec4i> cps_igm_e_inters; // (i,j) => cps edge verts, (k,l) => igm edge verts
    for(uint eid=0; eid<igm.num_edges(); ++eid)
    {
        if(igm.edge_is_boundary(eid)) continue;
        std::unordered_set<uint> ids;
        o.intersects_segment(igm.edge_verts(eid).data(), true, ids);
        for(uint id: ids)
        {
            if(cps.edge_is_boundary(id)) continue;
            //cps.edge_data(id).flags[MARKED] = true;
            vec4i inters;
            inters[0] = cps.edge_vert_id( id,0);
            inters[1] = cps.edge_vert_id( id,1);
            inters[2] = igm.edge_vert_id(eid,0);
            inters[3] = igm.edge_vert_id(eid,1);
            cps_igm_e_inters.push_back(inters);
        }
    }

    for(const vec4i & inters : cps_igm_e_inters)
    {
        vec3d p = segment_intersection(cps.vert(inters[0]),
                                       cps.vert(inters[1]),
                                       igm.vert(inters[2]),
                                       igm.vert(inters[3]));
        int eid = cps.edge_id(inters[0],inters[1]);
        if(eid>=0) cps.edge_split(eid, p);
    }
    cps.updateGL();


    {
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
    }

    GLcanvas gui;
    SurfaceMeshControls<DrawableTrimesh<>> controls(&cps,&gui,"IGM");
    gui.push(&cps);
//    gui.push(&igm);
    gui.push(&controls);

    return gui.launch();
}
