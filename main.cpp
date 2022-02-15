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
#include <cinolib/dijkstra.h>
#include <cinolib/connected_components.h>

using namespace cinolib;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void overlay_IGM_and_CPS(Quadmesh<> & igm,
                         Trimesh<>  & cps)
{
    // add IGM inner vertices inside the CPS
    Octree o;
    o.build_from_mesh_polys(cps);
    for(uint vid=0; vid<igm.num_verts(); ++vid)
    {
        if(igm.vert_is_boundary(vid)) continue;

        uint pid;
        if(o.contains(igm.vert(vid), true, pid))
        {
            cps.poly_split(pid, igm.vert(vid));
        }
    }

    // add IGM inner edges inside the CPS
    for(uint igm_eid=0; igm_eid<igm.num_edges(); ++igm_eid)
    {
        if(igm.edge_is_boundary(igm_eid)) continue;

        // I hate this, but since the content of the octree is dynamic (due to edge splits in the loop)
        // I have no other choice than making a new octree from scratch every time. Note that adding sub
        // segments in the tree is not a valid choice because  after splits edge ids will become inconsistent...
        Octree o;
        o.build_from_mesh_edges(cps);
        std::unordered_set<uint> eids;
        if(o.intersects_segment(igm.edge_verts(igm_eid).data(), true, eids))
        {
            std::set<uint,std::greater<uint>> ord_eids(eids.begin(),eids.end());
            for(uint cps_eid : ord_eids)
            {
                vec3d p = segment_intersection(cps.edge_vert(cps_eid,0),
                                               cps.edge_vert(cps_eid,1),
                                               igm.edge_vert(igm_eid,0),
                                               igm.edge_vert(igm_eid,1));

                double d0 = p.dist(cps.edge_vert(cps_eid,0));
                double d1 = p.dist(cps.edge_vert(cps_eid,1));
                // avoid creating tiny edges
                if(std::min(d0,d1)>1e-7)
                {
                    cps.edge_split(cps_eid,p);
                }
            }
        }
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void map_igm_cps(Quadmesh<> & igm,
                 Trimesh<>  & cps)
{
    overlay_IGM_and_CPS(igm,cps);

    std::vector<uint> v_igm2cps(igm.num_verts());
    // cut the mesh open in order to have each quad in the IGM
    // to be a separate domain  in the CPS (just for the sake of
    // texturing). Edges in CPS are created connecting corners
    // using Dijkstra's algorithm
    cps.edge_set_flag(MARKED,false);
    for(uint eid=0; eid<igm.num_edges(); ++eid)
    {
        if(igm.edge_is_boundary(eid)) continue;
        uint v0 = v_igm2cps.at(igm.edge_vert_id(eid,0));
        uint v1 = v_igm2cps.at(igm.edge_vert_id(eid,1));
        std::vector<uint> path;
        dijkstra(cps, v0, v1, path);
        for(uint i=0; i<path.size()-1; ++i)
        {
            int eid = cps.edge_id(path[i],path[i+1]);
            assert(eid>=0);
            cps.edge_data(eid).flags[MARKED] = true;
        }
    }
    cut_mesh_along_marked_edges(cps);

//    // label each connected component separately
//    std::vector<std::unordered_set<uint>> ccs;
//    uint n_labels = connected_components(cps, ccs);
//    std::cout << n_labels << " connected_components" << std::endl;
//    int l = 0;
//    for(auto & cc : ccs)
//    {
//        for(uint vid : cc)
//        {
//            cps.vert_data(vid).label = l;
//            for(uint pid : cps.adj_v2p(vid))
//            {
//                cps.poly_data(pid).label = l;
//            }
//        }
//        ++l;
//    }
//    cps.poly_color_wrt_label();

//    // apply bilinear texturing
//    std::vector<uint> lab2igm(n_labels,0);
//    for(uint pid=0; pid<igm.num_polys(); ++pid)
//    {
//        uint id;
//        cps_bvh.contains(igm.poly_centroid(pid), false, id);
//        lab2igm.at(cps.poly_data(id).label) = pid;
//    }
//    for(uint vid=0; vid<cps.num_verts(); ++vid)
//    {
//        uint pid = lab2igm.at(cps.vert_data(vid).label);

//        vec4d bary;
//        quad_barycentric_coords(igm.poly_vert(pid,0),
//                                igm.poly_vert(pid,1),
//                                igm.poly_vert(pid,2),
//                                igm.poly_vert(pid,3),
//                                cps.vert(vid), bary);

//        // quad in uv space
//        cps.vert_data(vid).uvw = vec3d(0,0,0) * bary[0] +
//                                 vec3d(1,0,0) * bary[1] +
//                                 vec3d(1,1,0) * bary[2] +
//                                 vec3d(0,1,0) * bary[3];
//    }
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

    GLcanvas gui;
    gui.push(&cps);
//    gui.push(&igm);
    gui.push(new SurfaceMeshControls<DrawableTrimesh<>> (&cps,&gui,"CPS"));
    gui.push(new SurfaceMeshControls<DrawableQuadmesh<>>(&igm,&gui,"IGM"));
    return gui.launch();
}
