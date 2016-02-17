gfx read node Cell_Geom.part0.exnode
gfx read ele Cell_Geom.part0.exelem

gfx define faces
gfx modify g_element "/" general clear;
gfx modify g_element "/" points domain_nodes coordinate Coordinate tessellation default_points LOCAL glyph point size "1*1*1" offset 0,0,0 font default select_on material default selected_material default_selected render_shaded;
gfx modify g_element "/" lines domain_mesh1d coordinate Coordinate face all tessellation default LOCAL line line_base_size 0 select_on material default selected_material default_selected render_shaded;

gfx create window 1 double_buffer;
gfx modify window 1 image scene "/" filter default infinite_viewer_lighting two_sided_lighting;
gfx modify window 1 image add_light default;
gfx modify window 1 image add_light default_ambient;
gfx modify window 1 layout 2d ortho_axes z -y eye_spacing 0.25 width 503 height 512;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 0 0 0 texture none;
gfx modify window 1 view parallel eye_point 0 0 28.8485 interest_point 0 0 0 up_vector -0 1 -0 view_angle 40 near_clipping_plane 1.44243 far_clipping_plane 58.547 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

