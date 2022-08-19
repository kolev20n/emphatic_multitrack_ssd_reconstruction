// This file contains global constants and data structures used for
// data exchange between classes, as well as input and output data
// structures.
// It also contains the line_3d class, which should belong to the
// partial_geometry class, because it is used as a data member of
// another class.
//
// Blame: Nikolay Kolev, kolev20n@uregina.ca
//
// 2022
//

#ifndef common_data_structures_and_functions_h
#define common_data_structures_and_functions_h

#include <vector>
#include <map>

#include <TMath.h>
#include <Math/Vector3D.h>

const double cm = 1.;     // units used throughout are centimeters
const double mm = 1.e-1;  // millimeter
const double um = 1.e-4; // micrometer

// the actual number of plates will be read from the geometry file
const int max_number_of_plaques = 32;
const int max_number_of_plates = 32;
const int max_number_of_plates_in_tracking_region = 32;
const int max_number_of_plates_in_group = 4;
const int max_number_of_groups_in_tracking_region = 16;
const int max_number_of_tracking_regions = 4;
const int max_number_of_groups = 16;
const int max_number_of_triplet_groups = 4;
const int max_number_of_tracks = 200;
const int max_number_of_plaques_in_a_plate = 4;
const int max_number_of_diagonal_plates = 4;
const int max_number_of_single_orientation_mistream_plates = 4; // may need to change in much more complex geometry
const int max_number_of_DD_groups_in_tracking_region = 1;
const int max_number_of_strips_on_plaque = 640;
const int max_midstream_options = 3;

const double pi = 3.141592653589793;
const double to_radians = pi / 180.;

class global_data_dispatcher;
class pattern;

extern global_data_dispatcher* the_global_data_dispatcher;

class line_3d
{
public:
  ROOT::Math::XYZVector point;
  ROOT::Math::XYZVector direction;
  ROOT::Math::XYZVector error_point;
  ROOT::Math::XYZVector error_direction;
};

enum geometry_type
{
  minimal, minimal_complete, redundant, complete
};

class reconstruction_algorithm
{
// The algorithm depends on three factors: various geometric factors, 
// absence or presence of magnetic field, and cluster multiplicity on
// the individual plates. Only the cluster mutliplicity is event-dependent,
// so the other two items will be set during initialization.
  
public: 
  geometry_type the_geometry_type;
  bool use_magnetic_field;
  int event_class_by_cluster_multiplicity;
  // 1100 single track events:
  //   1111 all plates have exactly 1 cluster
  //   1112 almost all have 1 cluster, except allowed by cuts (global and per region)
  // 2100 almost all are 1 upstream, group multiplicity is 1 or 2 midstream, 0 or 1 downstream, always non-increasing by group
  //   2101 midstream 2 and 2 and downstream 1 and 1
  // 2200 ditto, but more than 2 midstream (group) on some groups
  //   2201 midstream 3+ and 3+ and downstream 1 and 1
  //   2261 midstream 3+ and 2 and downstream 1 and 1
  // 2300 downstream exactly 2 and 2
  //   2301 midstream 2 and 2
  //   2331 midstream 3+ and 3+
  //   2361 midstream 3+ and 2
  // 5000 multitrack events
  // 9000 rejected events: multiple tracks upstream, or increasing by some criterion (secondary vertex), but make sure these cannot be reconstructed
  // also events with late secondary vertex; still can do partial tracks until then
  
  // Need better variables; alignment can be mutual and with respect to the intended directions.
  int alignment[max_number_of_tracking_regions][4];
};

class reconstruction_options
{
public:
  std::string cluster_interpretation_method; // default: weighted
  double maximum_cosangle_between_intended_and_actual_strip_direction;
  int max_upstream_plates_with_multiple_clusters;
  int max_clusters_on_upstream_plates_with_multiple_clusters;
  int max_midstream_plus_downstream_plates_with_multiple_clusters;
  int max_midstream_plus_downstream_plates_with_zero_clusters;
  int max_midstream_plus_downstream_plates_with_zero_or_multiple_clusters;
  int max_midstream_plus_downstream_plates_with_less_than_2_clusters;
  int max_plates_with_multiple_clusters_for_single_track;
  int max_upstream_plates_with_multiple_clusters_for_single_track;
  int max_midstream_plates_with_multiple_clusters_for_single_track;
  int max_downstream_plates_with_multiple_clusters_for_single_track;
  int max_upstream_plates_with_zero_clusters;
  
  double cut_max_distance_to_accept_in_track;
  double guessed_vertex_weight;
  double cut_max_distance_to_count_to_vertex;
  double cut_min_energy_to_be_high_energy_charged_particle;
  
  // options and cuts for event class 1111 and no magnetic field
  
};

class input_root_data_structure
{
public:
  input_root_data_structure();
  ~input_root_data_structure();
  
  // for the ssd energy tree
  // this is the only thing coming from experimental data
  //
  std::vector<int> *plaque_number;
  std::vector<int> *strip_number;
  std::vector<double> *total_energy_values;
  std::vector<double> *non_ionization_energy_values;
  std::vector<std::vector<double>> *contributing_tracks;
  
  // true particle/track information
  // for MC simulation analysis only
  //
  std::vector<int> *track_id;
  std::vector<int> *particle_code;
  std::vector<int> *parent_track_id;
  std::vector<std::string> *creator_process;
  std::vector<double> *x_vertex;
  std::vector<double> *y_vertex;
  std::vector<double> *z_vertex;
  std::vector<double> *px_vertex;
  std::vector<double> *py_vertex;
  std::vector<double> *pz_vertex;
  std::vector<double> *ekin_vertex;
  std::vector<std::vector<double>> *visited_plaques;
  
  // true step information
  // for MC simulation analysis only
  //
  std::vector<int> *step_track_id;
  std::vector<std::string> *step_volume_name;
  std::vector<int> *step_plaque_id;
  std::vector<double> *step_x_i;
  std::vector<double> *step_y_i;
  std::vector<double> *step_z_i;
  std::vector<double> *step_px_i;
  std::vector<double> *step_py_i;
  std::vector<double> *step_pz_i;
  std::vector<double> *step_x_f;
  std::vector<double> *step_y_f;
  std::vector<double> *step_z_f;
  std::vector<double> *step_px_f;
  std::vector<double> *step_py_f;
  std::vector<double> *step_pz_f;
  std::vector<double> *step_initial_energy;
  std::vector<double> *step_final_energy;
  std::vector<std::string> *step_process_name;
};

class hit_pair
{
public:
  int strip_number;
  double energy_value; // adc_value if experimental hit
};

class charged_track
{
public:
  int track_id;
  int particle_code;
  double ekin_vertex;
  std::vector<int> step_plate_index;
  std::vector<ROOT::Math::XYZVector> position;
  std::vector<double> energy_lost;
};

class mc_step
{
public:
  int step_index;
  int step_track_id;
  std::string step_volume_name;
  int step_plate_id;
  ROOT::Math::XYZVector initial_position;
  ROOT::Math::XYZVector final_position;
  double initial_energy;
  double final_energy;
  std::string step_process_name;
};

class run_mc_output_data_structure
{
public:
  int number_of_events;
  
  std::map<int, int> particle_codes_map;
};

class run_reco_output_data_structure
{
public:
  // geometry
  int number_of_plates;
  
  int number_of_events;
};

class event_mc_output_data_structure
{
public:
  ROOT::Math::XYZVector mc_vertex;
  
  std::vector<charged_track> the_charged_tracks; // high energy charged particles
  std::vector<mc_step> the_selected_steps;
  
  // void print();
};



bool compare_hit_pairs(const hit_pair& first, const hit_pair& second);

bool compare_map(const std::pair<std::string, int> &a, const std::pair<std::string, int> &b);

bool compare_cluster_line_plate(const std::pair<line_3d, int> &a, const std::pair<line_3d, int> &b);

bool compare_pattern_quality_0(const pattern& first, const pattern& second);

#endif
