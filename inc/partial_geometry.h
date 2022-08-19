// The partial_geometry class provides all the geometry functionality for
// the silicon strip detectors for EMPHATIC. It is partial
// because there are of course more detectors in the experiment.
// This file also contains some utility data classes used in the geometry
// (but not the line_3d class which is defined in common_data_structures and functions),
// as well as some utility functions such as for finding the intersection of
// a straight line and a plane, etc.
//
// Blame: Nikolay Kolev, kolev20n@uregina.ca
//
// 2022
//

#ifndef partial_geometry_h
#define partial_geometry_h

//#include <vector>
//#include <string>

#include <TMath.h>
#include <Math/Vector3D.h>

#include "common_data_structures_and_functions.h"

class plane_3d
{
public:
  // the coefficients in the plane equation ax + by + cz + d = 0
  double a, b, c, d;
};

class ssd_plaque
{
public:
  int plaque_id; // for reference purposes only
  int tracking_region;
  int plate;
  int group;
  int number_of_strips;
  std::string type;
  double efficiency;
  ROOT::Math::XYZVector position; // position in cm
  ROOT::Math::XYZVector rotation; // rotations in degrees
  ROOT::Math::XYZVector size; // size in cm: along strips, perpendicular to strips, thickness
  ROOT::Math::XYZVector strip_direction; // strip direction unit vector
  ROOT::Math::XYZVector normal_direction; // normal direction unit vector
  // ROOT::Math::XYZVector intended_strip_direction; // remove later; this will be the plate_type
  double pitch;
  plane_3d plaque_plane_equation;
};

class ssd_plate
{
public:
  int number_of_member_plaques;
  int member_plaques[max_number_of_plaques_in_a_plate];
  ROOT::Math::XYZVector position; // this is a hack so previous code doesn't break
  int group;
  int tracking_region;
  std::string intended_to_measure; // can be X, Y, D1 or D2
};

class plate_group
{
// Groups are XY or D1D2. Individual D plates do not need to be assigned to a group.
// For the moment groups of single plates are allowed, so there can be X groups and Y groups.
public:
  int number_of_plates_in_group; // if X and Y groups are allowed this is needed
  int plate_index_list[max_number_of_plates_in_group]; // maybe have the actual plates?
  
  int number_of_attached_d_plates; // can be 0 or 1
  int attached_d_plate_index; // this is optional, to keep track; only individual D plates can be attached
  
  std::string group_type;
  // XY, D1D2, X, Y (the last 2 may not be needed)
  
  int tracking_region; // 0 for upstream, 1 for midstream, 2 for downstream
};

class tracking_region
{
// tracking region is a class that provides convenient (for speed) but otherwise redundant
// information about the groups and plates in a tracking region.
public:
  std::string name;
  int number_of_groups;
  int group_index_list[max_number_of_groups_in_tracking_region];
  int number_of_plates;
  int plate_index_list[max_number_of_plates_in_tracking_region];
};

class target
{
// Obvously, the target information.
public:
  ROOT::Math::XYZVector position; // position in cm
  ROOT::Math::XYZVector rotation; // rotations in degrees
  ROOT::Math::XYZVector size; // size in cm
  std::string material; // none, carbon, iron, aluminum
};

class track_point
{
public:
  ROOT::Math::XYZVector position;
  ROOT::Math::XYZVector momentum;
  double momentum_magnitude;
};

class track
{
// A class that bundles together all the information about a reconstructed track:
// clusters that belong to it, track parameters (a stright line in the no-field case)
// start and end group and a sorting function for the cluster lines (because they)
// may not be assigned in any particular order.
public:
  std::vector<line_3d> cluster_lines;
  std::vector<int> plate_index_of_cluster_line;
  
  line_3d track_line;
  
  bool projects_downstream;
  
  // the following are the position and direction at z = 0;
  // the track may not start or end there;
  // this is irrelevant for the calculation though
  // ROOT::Math::XYZVector position_at_0;
  // ROOT::Math::XYZVector direction_at_0;
  
  double fit_goodness;
  
  // the track shouldn't be extrapolated outside
  // the following z coordinates
  double assumed_start_group;
  double assumed_end_group;
  
  //double momentum_magnitude;
  std::vector<track_point> track_points;
  
  // the following concern the case when there is a missing hit on an upstream plate,
  // so do later; need to calculate these in the branches of find_upstream_track,
  // specifically all fit_track_no_field calls
  ROOT::Math::XYZVector position_on_target; // at z = 0
  ROOT::Math::XYZVector error_in_position_on_target; // at z = 0
  
  ROOT::Math::XYZVector track_plate_intersection_points[max_number_of_plates]; // set to 0, 0, 0 if not hit
  
  void sort_cluster_lines_by_plate_index();
  bool find_intersection_points_no_field(); // used for no field and in field for upstream and maybe midstream regions; uses the assumed start and end groups
  bool find_intersection_points_rough(); // using first plaque position on multi-plates, which is approximate
};

class track_options
{
public:
  double option_quality;
  std::vector<track> the_tracks;
};

class multi_track_selections
{
public:
  int number_of_tracks;
  double selection_quality;
  int index_x_track_selections[200]; // use const variable later
  int index_y_track_selections[200]; // ditto
};

class event_reco_output_data_structure
{
public:
  std::vector<line_3d> cluster_lines[max_number_of_plates];
  
  int event_class_by_cluster_multiplicity;
  // 1111 all plates have exactly 1 cluster
  // 1112 almost all have 1 cluster, except allowed by cuts (global and per region)
  // 2100 almost all are 1 upstream, group multiplicity is 1 or 2 midstream, 0 or 1 downstream, always non-increasing by group
  // 2200 ditto, but more than 2 midstream (group) on some groups
  // 5000 multitrack events
  // 9000 rejected events: multiple tracks upstream, or increasing by some criterion (secondary vertex), but make sure these cannot be reconstructed
  // also events with late secondary vertex; still can do partial tracks until then
  
  ROOT::Math::XYZVector guessed_vertex;
  
  int group_cluster_multiplicity[max_number_of_groups];
  
  ROOT::Math::XYZVector vertex;
  std::vector<ROOT::Math::XYZVector> track_directions;
  double global_fit_quality;
  
  std::vector<track> reconstructed_tracks;
  
  track the_upstream_track;
  
  // void print();
};

class partial_geometry
{
// The central and important silicon strip plates geometry class.
public:
  partial_geometry();

// number getters
  int get_number_of_plaques();
  int get_number_of_plates();
  int get_number_of_plates_in_tracking_region(int a_region);
  int get_number_of_tracking_regions();
  int get_number_of_groups();
  int get_number_of_groups_in_tracking_region(int a_region);
  int get_number_of_plates_in_group(int a_group);
  int get_number_of_upstream_plates();
  int get_number_of_midstream_plates();
  int get_number_of_downstream_plates();
  int get_number_of_midstream_plus_downstream_plates();
  int get_number_of_upstream_groups();
  int get_number_of_midstream_groups();
  int get_number_of_downstream_groups();
  int get_number_of_midstream_plus_downstream_groups();
  
  int get_number_of_midstream_diagonal_plates();
  int get_number_of_downstream_diagonal_plates();
  
  target* get_target();
  
// getters for the primary and secondary direction and their names
  ROOT::Math::XYZVector get_primary_direction_1();
  ROOT::Math::XYZVector get_primary_direction_2();
  ROOT::Math::XYZVector get_secondary_direction_1();
  ROOT::Math::XYZVector get_secondary_direction_2();
  std::string get_primary_direction_1_name();
  std::string get_primary_direction_2_name();
  std::string get_secondary_direction_1_name();
  std::string get_secondary_direction_2_name();

// handles to groups and plates
  // tracking_region* get_tracking_region(int a_region);
  plate_group* get_plate_group(int a_group);
  ssd_plate* get_plate(int a_plate_index);
  ssd_plate* get_plate_by_id(int a_plate_id); // for reference purposes
  ssd_plaque* get_plaque(int a_plaque_index);
  
// adding stuff to the geometry
  bool add_plaque(ssd_plaque a_plaque);
  // bool add_plate(ssd_plate a_plate);
  bool add_target(target a_target);
  
// setters for the primary and secondary direction and their names
  bool set_primary_direction_1(ROOT::Math::XYZVector a_direction);
  bool set_primary_direction_2(ROOT::Math::XYZVector a_direction);
  bool set_secondary_direction_1(ROOT::Math::XYZVector a_direction);
  bool set_secondary_direction_2(ROOT::Math::XYZVector a_direction);
  bool set_primary_direction_1_name(std::string a_name);
  bool set_primary_direction_2_name(std::string a_name);
  bool set_secondary_direction_1_name(std::string a_name);
  bool set_secondary_direction_2_name(std::string a_name);
  
// setter for tracking region names
  bool set_tracking_region_names(int a_region, std::string a_name);
  
// plate classification happens during initialization of the geometry
  bool classify_plates();
  
// magentic field getters and setters
  ROOT::Math::XYZVector get_magnetic_field_general_direction();
  bool set_magnetic_field_general_direction(ROOT::Math::XYZVector a_direction);
  bool set_magnetic_field_general_direction_name(std::string a_magnetic_field_general_direction_name);
  
// getters for the extent of the target at x = 0 and y = 0
  double get_target_min_z();
  double get_target_max_z();
  
// getters for the midstream X and sY plates and the X and Y scales
  bool get_midstream_X_plates(int& a_number_of_midstream_X_plates, int a_midstream_X_plates[]);
  bool get_midstream_Y_plates(int& a_number_of_midstream_Y_plates, int a_midstream_Y_plates[]);
  bool get_midstream_X_scales(double a_midstream_X_scales_max[], double a_midstream_X_scales_min[]);
  bool get_midstream_Y_scales(double a_midstream_Y_scales_max[], double a_midstream_Y_scales_min[]);
  bool get_upstream_XY_groups(int& a_number_of_upstream_XY_groups, int a_upstream_XY_groups[]); // this includes triple groups
  bool get_midstream_XY_groups(int& a_number_of_midstream_XY_groups, int a_midstream_XY_groups[]); // this also includes triple groups
  bool get_downstream_XY_groups(int& a_number_of_downstream_XY_groups, int a_downstream_XY_groups[]); // this also includes triple groups
  bool get_types_of_midstream_plates(int& number_of_x_plates, ssd_plate x_plates[], int x_plate_indices[], int& number_of_y_plates, ssd_plate y_plates[], int y_plate_indices[], int& number_of_d1_plates, ssd_plate d1_plates[], int d1_plate_indices[], int& number_of_d2_plates, ssd_plate d2_plates[], int d2_plate_indices[], int& number_of_d_plates, ssd_plate d_plates[], int d_plate_indices[]);
  
private:
  int actual_number_of_plaques;
  ssd_plaque the_plaques[max_number_of_plaques];
  
  int actual_number_of_plates;
  ssd_plate the_plates[max_number_of_plates];
  
  int number_of_tracking_regions;
  tracking_region the_tracking_regions[max_number_of_tracking_regions];
  std::string tracking_region_names[max_number_of_tracking_regions];

  int number_of_groups;
  plate_group the_plate_groups[max_number_of_groups];
  
  int number_of_triple_groups;
  int the_triple_groups[max_number_of_groups];
  
  int number_of_diagonal_plates;
  int diagonal_plates[max_number_of_diagonal_plates];
  
  int number_of_d1_plates_midstream;
  int d1_plates_midstream[max_number_of_diagonal_plates];
  
  int number_of_d2_plates_midstream;
  int d2_plates_midstream[max_number_of_diagonal_plates];
  
  int number_of_diagonal_plates_midstream;
  int diagonal_plates_midstream[max_number_of_diagonal_plates];
  
  int number_of_diagonal_plates_downstream;
  int diagonal_plates_downstream[max_number_of_diagonal_plates];
  
  ROOT::Math::XYZVector primary_direction_1;
  std::string primary_direction_1_name;
  
  ROOT::Math::XYZVector primary_direction_2;
  std::string primary_direction_2_name;
  
  ROOT::Math::XYZVector secondary_direction_1;
  std::string secondary_direction_1_name;
  
  ROOT::Math::XYZVector secondary_direction_2;
  std::string secondary_direction_2_name;

  target the_target;
  
  ROOT::Math::XYZVector magnetic_field_general_direction;
  std::string magnetic_field_general_direction_name;
  
  int number_of_midstream_X_plates;
  int midstream_X_plates[max_number_of_groups_in_tracking_region];
  // for the scales use only the above the diagonal elements of the matrix, e.g. (0, 1), (0, 2) and (1, 2) (if they are three)
  // or just (0, 1) if they are two
  double midstream_X_scales_max[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];
  double midstream_X_scales_min[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];

  int number_of_midstream_Y_plates;
  int midstream_Y_plates[max_number_of_groups_in_tracking_region];
  double midstream_Y_scales_max[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];
  double midstream_Y_scales_min[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];

  int number_of_upstream_XY_groups; // this includes triple groups
  int upstream_XY_groups[max_number_of_groups_in_tracking_region];
  int number_of_midstream_XY_groups; // this also includes triple groups
  int midstream_XY_groups[max_number_of_groups_in_tracking_region];
  int number_of_downstream_XY_groups; // this also includes triple groups
  int downstream_XY_groups[max_number_of_groups_in_tracking_region];
  
  int number_of_upstream_DD_groups;
  int upstream_DD_groups[max_number_of_DD_groups_in_tracking_region];
  int number_of_midstream_DD_groups;
  int midstream_DD_groups[max_number_of_DD_groups_in_tracking_region];
  int number_of_downstream_DD_groups;
  int downstream_DD_groups[max_number_of_DD_groups_in_tracking_region];

  bool check_cosangle(int a_plaque_index);
};

class pattern
{
public:
  int number_of_x_lines;
  int number_of_y_lines;
  
  ROOT::Math::XYZVector guessed_vertex_x;
  ROOT::Math::XYZVector guessed_vertex_y;
  ROOT::Math::XYZVector common_guessed_vertex;
  ROOT::Math::XYZVector vertex;

  double scale_x;
  double scale_y;
  
  std::vector<line_3d> equations_x;
  std::vector<line_3d> equations_y;
  std::vector<ROOT::Math::XYZVector> pattern_x_coordinates;
  std::vector<ROOT::Math::XYZVector> pattern_y_coordinates;
  
  std::vector<track> x_tracks;
  std::vector<track> y_tracks;
  std::vector<track> combined_tracks;
  
  double pattern_quality[5]; 
  
  void reset();
};

// finds the intersection point of a line and a plane
// returns false if plane and line are parallel, even if line is in plane
bool find_line_plate_intersection_point(line_3d a_track_line, plane_3d plate_plane, ROOT::Math::XYZVector& intersection_point);

// this distance always exists
double find_point_to_line_3d_distance(ROOT::Math::XYZVector a_point, line_3d a_line);

// this distance also always exists; needs reworking (including return double for the distance)
bool find_line_to_line_distance(line_3d line_1, line_3d line_2, double& distance);

// find the midpoint of the segment connecting the two closest points in two skew lines;
// if the lines intersect, this is the intersection point; should be reworked to return 
// false if lines are parallel (which shouldn't happen in the cases this is used)
bool find_closest_point(line_3d line_1, line_3d line_2, ROOT::Math::XYZVector& closest_point);

// two parallel lines always determine a plane, unless they coincide; rework for that, although
// it should not happen in the cases used
bool find_plane_of_two_parallel_lines(line_3d line_1, line_3d line_2, plane_3d& a_plane);

void find_3d_line_from_projections(double x[], double z_x[], double y[], double z_y[], ROOT::Math::XYZVector& point, ROOT::Math::XYZVector& direction);

// subspace is xy or yz or zx (and cyclic)
void find_2d_plane_equation(std::string subspace, ROOT::Math::XYZVector point_1, ROOT::Math::XYZVector point_2, line_3d the_equation_line);

ROOT::Math::XYZVector find_point_on_line_by_z(line_3d a_line, double a_z);

void combine_3d_line_from_projections(track x_track, track y_track, line_3d& a_line);

bool compare_multi_track_selections(const multi_track_selections& first, const multi_track_selections& second);

bool is_intersection_point_inside_plaque(ROOT::Math::XYZVector intersection_point, ssd_plaque* p_dummy_plaque);

#endif
