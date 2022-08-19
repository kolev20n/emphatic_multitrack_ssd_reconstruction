// This class assigns clusters to tracks and then fits the tracks.
// This is the core of the ssd reconstruction package.
//
// Blame: Nikolay Kolev, kolev20n@uregina.ca
//
// 2022
//

#ifndef track_finding_manager_h
#define track_finding_manager_h

#include <vector>

#include "common_data_structures_and_functions.h"
#include "partial_geometry.h"

class track_finding_manager
{
public:
  track_finding_manager();
  bool find_tracks(std::vector<line_3d> cluster_lines[], std::vector<track>& the_tracks);
  
private:
  partial_geometry* the_geometry;
  reconstruction_algorithm* the_reconstruction_algorithm;
  
  bool classify_event(std::vector<line_3d> cluster_lines[]);
  bool fit_track_no_field(track& a_track);
  // bool fit_track_in_field(track& a_track);
  bool find_upstream_track(std::vector<line_3d> cluster_lines[], track& upstream_track, int event_class);
  bool find_midstream_plus_downstream_track_no_field(std::vector<line_3d> cluster_lines[], track& middownstream_track, int event_class);
  bool find_total_track_no_field(std::vector<line_3d> cluster_lines[], track& total_track, int event_class);
  bool assign_clusters_midstream_no_field(std::vector<line_3d> cluster_lines[], std::vector<pattern>& pattern_options, int event_class, track upstream_track);
  bool find_midstream_tracks_no_field(std::vector<pattern>& pattern_options, int event_class, track upstream_track);
  bool select_best_multi_track_options(double multi_track_best_distance[][200], int number_of_x_tracks, int number_of_y_tracks, std::vector<multi_track_selections>& the_multi_track_selections);
  bool propagate_tracks_downstream(std::vector<pattern>& pattern_options);
  bool assign_clusters_downstream_no_field(std::vector<pattern>& pattern_options, std::vector<line_3d> cluster_lines[]);
  bool do_global_fit_no_field(std::vector<pattern>& pattern_options, track upstream_track);
  bool assign_clusters_single_track(std::vector<line_3d> cluster_lines[], pattern upstream_single_track_pattern, pattern midstream_plus_downstream_single_track_pattern, pattern total_single_track_pattern);
  bool do_fit_single_track(track a_track);
  bool assign_clusters_downstream(std::vector<pattern>& pattern_options, std::vector<line_3d> cluster_lines[]);
  bool do_global_fit(std::vector<pattern>& pattern_options, track upstream_track);

};

#endif
