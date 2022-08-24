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
                         Trimesh<>  & cps,
                         Trimesh<>  & obj)
{
    // add IGM inner vertices inside the CPS
    Octree o;
    o.build_from_mesh_polys(cps);
    for(uint vid=0; vid<igm.num_verts(); ++vid)
    {
        if(igm.vert_is_boundary(vid)) continue;

        uint pid;
        if(o.contains(igm.vert(vid), false, pid))
        {
            // NOTE: this is save as long as each triangle in CPS
            // contains at most one vertex of the IGM. If two points
            // are on the same pid, the sub-triangle that inherits
            // the pid may not be the one containing the second point.
            // Here I am assuming this will never happen...
            cps.poly_split(pid, igm.vert(vid));

            // replicate the same move in object space
            double w[3];
            triangle_barycentric_coords(cps.poly_vert(pid,0),
                                        cps.poly_vert(pid,1),
                                        cps.poly_vert(pid,2),
                                        igm.vert(vid), w);

            obj.poly_split(pid, w[0]*obj.poly_vert(pid,0)+
                                w[1]*obj.poly_vert(pid,1)+
                                w[2]*obj.poly_vert(pid,2));
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
                    // replicate move in object space
                    float t = d0/(d0+d1);
                    obj.edge_split(cps_eid, t);
                    cps.edge_split(cps_eid, t);
                }
            }
        }
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void bilinear_texturing(Quadmesh<> & igm,
                        Trimesh<>  & cps,
                        Trimesh<>  & obj)
{
    // create seams along IGM edges so that
    // each quad can be textured separately
    Octree o;
    o.build_from_mesh_points(cps);
    std::vector<uint> vmap(igm.num_verts()); // IGM => CPS vmap
    for(uint vid=0; vid<igm.num_verts(); ++vid)
    {
        uint   cps_vid;
        vec3d  pos;
        double dist;
        o.closest_point(igm.vert(vid), cps_vid, pos, dist);
        vmap.at(vid) = cps_vid;
    }
    cps.edge_set_flag(MARKED,false);
    obj.edge_set_flag(MARKED,false);
    for(uint eid=0; eid<igm.num_edges(); ++eid)
    {
        if(igm.edge_is_boundary(eid)) continue;
        uint v0 = vmap.at(igm.edge_vert_id(eid,0));
        uint v1 = vmap.at(igm.edge_vert_id(eid,1));
        std::vector<uint> path;
        dijkstra(cps, v0, v1, path);
        for(uint i=0; i<path.size()-1; ++i)
        {
            int eid = cps.edge_id(path[i],path[i+1]);
            assert(eid>=0);
            cps.edge_data(eid).flags[MARKED] = true;
            obj.edge_data(eid).flags[MARKED] = true;
        }
    }
    cut_mesh_along_marked_edges(cps);
    cut_mesh_along_marked_edges(obj);

    // use labels to distinguish IGM domains
    std::vector<std::unordered_set<uint>> ccs;
    uint n_labels = connected_components(cps, ccs);
    std::cout << n_labels << " connected_components - " << igm.num_polys() << " IGM domains " << std::endl;
    int l = 0;
    for(auto & cc : ccs)
    {
        for(uint vid : cc)
        {
            cps.vert_data(vid).label = l;
            obj.vert_data(vid).label = l;
            for(uint pid : cps.adj_v2p(vid))
            {
                cps.poly_data(pid).label = l;
                obj.poly_data(pid).label = l;
            }
        }
        ++l;
    }
    cps.poly_color_wrt_label();
    obj.poly_color_wrt_label();

    // apply bilinear texturing
    Octree o_polys;
    o_polys.build_from_mesh_polys(cps);
    std::vector<uint> lab2igm(n_labels,0);
    for(uint pid=0; pid<igm.num_polys(); ++pid)
    {
        uint id;
        o_polys.contains(igm.poly_centroid(pid), false, id);
        lab2igm.at(cps.poly_data(id).label) = pid;
    }
    for(uint vid=0; vid<cps.num_verts(); ++vid)
    {
        uint pid = lab2igm.at(cps.vert_data(vid).label);

        vec4d bary;
        quad_barycentric_coords(igm.poly_vert(pid,0),
                                igm.poly_vert(pid,1),
                                igm.poly_vert(pid,2),
                                igm.poly_vert(pid,3),
                                cps.vert(vid), bary);

        // quad in uv space
        vec3d uvw = vec3d(0,0,0) * bary[0] +
                    vec3d(1,0,0) * bary[1] +
                    vec3d(1,1,0) * bary[2] +
                    vec3d(0,1,0) * bary[3];

        cps.vert_data(vid).uvw = uvw;
        obj.vert_data(vid).uvw = uvw;
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void topological_IGM(const uint genus, Quadmesh<> & igm)
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
    std::string s = "/Users/cino/Desktop/Screen_Bracket.obj";
    //std::string s = "/Users/cino/Desktop/triblader-new-with-holes-3.5mm.stl";
    //std::string s = "/Users/cino/Desktop/tmp_data/eight.off";
    DrawableTrimesh<> obj(s.c_str());
    DrawableTrimesh<> cps;
    uint genus = obj.genus();

    DrawableSegmentSoup obj_loops;
    DrawableSegmentSoup cps_edges;
    obj_loops.draw_joint_spheres = true;
    obj_loops.joint_sphere_subd = 2;
    cps_edges.draw_joint_spheres = true;
    cps_edges.joint_sphere_subd = 2;
    cps_edges.thickness = 7;

    // compute homotopy basis
    HomotopyBasisData data;
    data.root              = 0;
    data.globally_shortest = false;
    data.detach_loops      = true;
    data.split_strategy    = EDGE_SPLIT_STRATEGY;

    auto IGM = [&]()
    {
        homotopy_basis(obj, data);
        canonical_polygonal_schema(obj, data, cps);
        std::cout << data << std::endl;
        // color code basis loops/CPS edges...
        obj_loops.clear();
        cps_edges.clear();
        obj_loops.colors.clear();
        cps_edges.colors.clear();
        for(uint eid=0; eid<obj.num_edges(); ++eid)
        {
            if(obj.edge_is_boundary(eid))
            {
                uint v0 = obj.edge_vert_id(eid,0);
                uint v1 = obj.edge_vert_id(eid,1);
                int loop_id = obj.vert_data(v0).label;
                if(loop_id<0) loop_id = obj.vert_data(v1).label;
                Color c = Color::scatter(data.loops.size(),loop_id,1,1);
                obj_loops.push_seg(obj.vert(v0),obj.vert(v1),c);
                cps_edges.push_seg(cps.vert(v0),cps.vert(v1),c);
            }
        }
        // overlay IGM template and do the IGM mapping
        if(genus>1)
        {
            DrawableQuadmesh<> igm;
            topological_IGM(genus,igm);
            igm.show_mesh_points();
            overlay_IGM_and_CPS(igm,cps,obj);
            bilinear_texturing(igm,cps,obj);
            obj.edge_mark_boundaries();
            cps.edge_mark_boundaries();
        }
    };

    IGM();

    GLcanvas gui_obj;
    GLcanvas gui_cps;
    obj.show_texture2D(TEXTURE_2D_CHECKERBOARD, 1);
    obj.show_wireframe(false);
    obj.updateGL();
    gui_obj.push(&obj);
    gui_obj.push(&obj_loops);
    cps.show_texture2D(TEXTURE_2D_CHECKERBOARD, 1);
    cps.show_wireframe(false);
    cps.updateGL();
    gui_cps.push(&cps);
    gui_cps.push(&cps_edges);

    gui_obj.push(new SurfaceMeshControls<DrawableTrimesh<>> (&obj,&gui_obj,"OBJ"));
    gui_obj.push(new SurfaceMeshControls<DrawableTrimesh<>> (&cps,&gui_obj,"CPS"));

    gui_obj.callback_mouse_left_click = [&](int modifiers) -> bool
    {
        if(modifiers & GLFW_MOD_SHIFT)
        {
            vec3d p;
            vec2d click = gui_obj.cursor_pos();
            if(gui_obj.unproject(click, p)) // transform click in a 3d point
            {
                obj.load(s.c_str());
                data.root = obj.pick_vert(p);
                std::cout << "Center Homotopy Basis at Vertex ID " << data.root << std::endl;
                IGM();
                obj.updateGL();
                cps.updateGL();
            }
        }
        return false;
    };

    // allow users to interactively edit cut thickness and color
    gui_obj.callback_app_controls = [&]()
    {
        ImGui::Text("Cut Thicknesses");
        ImGui::SliderFloat("#thickobj", &obj_loops.thickness,0,30);
        ImGui::SliderFloat("#thickcps", &cps_edges.thickness,0,30);
    };

    return gui_obj.launch({&gui_cps});
}
