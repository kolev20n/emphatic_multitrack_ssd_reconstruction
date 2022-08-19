#include <iostream>
#include <utility>
#include <algorithm>
#include <math.h>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"

#include "track_finding_manager.h"
#include "partial_geometry.h"
#include "global_data_dispatcher.h"
#include "common_data_structures_and_functions.h"
#include "my_straight_line_fcn.h"
#include "my_straight_line_fcn_with_point.h"
#include "my_global_fit_no_field_fcn.h"
#include "matej_minimum_data_structures_and_functions.h"

using namespace std;
using namespace ROOT::Math;
using namespace ROOT::Minuit2;

track_finding_manager::track_finding_manager()
:the_geometry(0), the_reconstruction_algorithm(0)
{
  cout << "Created the_track_finding_manager." << endl;

  the_geometry = the_global_data_dispatcher->get_partial_geometry();
  
  the_reconstruction_algorithm = the_global_data_dispatcher->get_reconstruction_algorithm();
}

bool track_finding_manager::find_tracks(vector<line_3d> cluster_lines[], vector<track>& the_tracks)
{
  int event_class;
  
  bool operation_success;
  
  double fit_goodness[max_number_of_tracks]; // this is a bad size of the array
  
  track upstream_track, downstream_track, total_track;
  
  double dummy_double;
  
  vector<pattern> pattern_options; // allow the possibility here, but assign_clusters_midstream_no_field or assign_clusters_midstream should select a single one
  
  pattern upstream_single_track_pattern, midstream_plus_downstream_single_track_pattern, total_single_track_pattern;
  
  operation_success = classify_event(cluster_lines);

  if (!operation_success) return operation_success;
  
  event_class = the_global_data_dispatcher->get_event_reco_output_data_structure()->event_class_by_cluster_multiplicity;
  
  if (!(event_class < 9000))
  {
    return false;
  }
  
  pattern_options.clear();
  
cout << event_class << endl;
  
  if (the_global_data_dispatcher->get_reconstruction_algorithm()->use_magnetic_field == false)
  {
    if (event_class == 1111 || event_class == 1112)
    {
      // event_class is 1111
      //
      // fit a straight line through all cluster lines; also fit two straight lines upstream and midstream plus downstream
      // compare which one is better; if the first one, then nothing happened; if the second one find the vertex;
      // if both are bad (by the chi square of the fit and the applied cut) either ignore or look for vertex on a plate (determined by option)
      // by fitting separately to and from that plate; the possibility exists that there can be a double vertex, which can be looked for too
      // (it can be target/plate or plate/plate case, the latter pretty rare), but perhaps ignore for now; this does not depend on alignment,
      // but alignment may affect errors and correlations
      //
      //
      // event_class is 1112
      //
      // fit a straight line through all the one-cluster plates, on the multiple cluster-planes add the one that is the closest (if both are really close use average?)
      // and fit again; do the same for midstream plus downstream; for upstream it can have one cluster on all (so it's done then), but it can have multiple clusters
      // on say one plate or zero clusters on say one plate; if multiple choose the one that produces a better vertex (smaller distance between the two fitted skew lines),
      // if zero determine the vertex by intersection of the plane formed by the other two lines and the midstream/downstream fit; this depends on the alignment of these
      // two strip lines, so will need to correct using the position on the third plate (such a procedure will be necessary in general); with multiple clusters upstream can
      // also use the one that will make the direction the closest to the beam direction
    
      // !!!! this needs to be reworked altogether !!!!
      
      operation_success = find_upstream_track(cluster_lines, upstream_track, event_class);
    
      fit_goodness[0] = upstream_track.fit_goodness;
    
      operation_success = find_midstream_plus_downstream_track_no_field(cluster_lines, downstream_track, event_class);
    
      fit_goodness[1] = downstream_track.fit_goodness;
    
      operation_success = find_total_track_no_field(cluster_lines, total_track, event_class);
    
      fit_goodness[2] = total_track.fit_goodness;
    
      XYZVector temp_vertex;
      find_closest_point(upstream_track.track_line, downstream_track.track_line, temp_vertex);
    
      cout << "Fit goodness: " << fit_goodness[0] << " " << fit_goodness[1] << " " << fit_goodness[2] << endl;
      cout << "Temp vertex: " << temp_vertex.X() << " " << temp_vertex.Y() << " " << temp_vertex.Z() << endl;
      
      /*
      if (event_class == 1111 && the_global_data_dispatcher->get_input_root_data_structure()->track_id->size() > 1)
      {
        cout << "from fit: " << upstream_track.track_line.direction << " " << downstream_track.track_line.direction << endl;
        cout << "-----------" << endl;
      }
      */
      
      //if (event_class == 1111 && the_global_data_dispatcher->get_input_root_data_structure()->track_id->size() > 1)
      /*
      if (event_class == 1111 && the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.size() > 1)
      {
        upstream_track.track_line.direction = upstream_track.track_line.direction.Unit();
        downstream_track.track_line.direction = downstream_track.track_line.direction.Unit();
      
        cout << "from fit: " << upstream_track.track_line.direction.Dot(downstream_track.track_line.direction) << endl;
      }
      */
      
      /*
      if (event_class == 1111)
      {
        cout << "--------------------" << endl;
        
        for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates(); i++)
        {
          dummy_double = the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->position.Z();
        
          cout << upstream_track.track_line.point.X() + (dummy_double - upstream_track.track_line.point.Z()) * upstream_track.track_line.direction.X() / upstream_track.track_line.direction.Z() << " "
               << upstream_track.track_line.point.Y() + (dummy_double - upstream_track.track_line.point.Z()) * upstream_track.track_line.direction.Y() / upstream_track.track_line.direction.Z() << " "
               << dummy_double << endl;
        }
        for (int i = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates(); i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); i++)
        {
          dummy_double = the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->position.Z();
          
          cout << downstream_track.track_line.point.X() + (dummy_double - downstream_track.track_line.point.Z()) * downstream_track.track_line.direction.X() / downstream_track.track_line.direction.Z() << " "
               << downstream_track.track_line.point.Y() + (dummy_double - downstream_track.track_line.point.Z()) * downstream_track.track_line.direction.Y() / downstream_track.track_line.direction.Z() << " "
               << dummy_double << endl;
        }
      }
      */
      
    }
    else if (event_class == 1117 || event_class == 2000)
    {
      // event_class is 1117
      //
      // look for one track midstream, but there are more than the allowed plates with multiple clusters (which
      // would qualify it for class 1112)
      //
      // event_class is 2000
      //
      // look for at least two tracks midstream
      
      operation_success = find_upstream_track(cluster_lines, upstream_track, event_class);
      
cout << "Fit goodness: " << upstream_track.fit_goodness << endl;

      if (operation_success)
      {
        operation_success = assign_clusters_midstream_no_field(cluster_lines, pattern_options, event_class, upstream_track);
      }

      if (operation_success)
      {
        operation_success = find_midstream_tracks_no_field(pattern_options, event_class, upstream_track);
      }
  
      /*
      // this is not needed in the no-field case, because there is no magnet
      if (operation_success)
      {
        operation_success = propagate_tracks_downstream(pattern_options);
      }
      */
      
      if (operation_success)
      {
        operation_success = assign_clusters_downstream_no_field(pattern_options, cluster_lines);
      }

      if (operation_success)
      {
        operation_success = do_global_fit_no_field(pattern_options, upstream_track);
      }
    }
  } // end if on magnetic field
  else
  {
    if (event_class == 1111 || event_class == 1112)
    {
      operation_success = find_upstream_track(cluster_lines, upstream_track, event_class);
      
cout << "Fit goodness: " << upstream_track.fit_goodness << endl;
      
      if (event_class == 1111)
      {
        // assign clusters here to the three tracks
        operation_success = do_fit_single_track(upstream_single_track_pattern.combined_tracks.at(0));
        operation_success = do_fit_single_track(midstream_plus_downstream_single_track_pattern.combined_tracks.at(0));
        operation_success = do_fit_single_track(total_single_track_pattern.combined_tracks.at(0));
      }
      else if (event_class == 1112)
      {
        operation_success = assign_clusters_single_track(cluster_lines, upstream_single_track_pattern, midstream_plus_downstream_single_track_pattern, total_single_track_pattern);
        
        if (operation_success)
        {
          operation_success = do_fit_single_track(upstream_single_track_pattern.combined_tracks.at(0));
          operation_success = do_fit_single_track(midstream_plus_downstream_single_track_pattern.combined_tracks.at(0));
          operation_success = do_fit_single_track(total_single_track_pattern.combined_tracks.at(0));

          if (total_single_track_pattern.pattern_quality[0] > 2.0 * (upstream_single_track_pattern.pattern_quality[0] + midstream_plus_downstream_single_track_pattern.pattern_quality[0]) ) // use cut later
          {
            // this is a vertex case; need a variable for the event to state this
          }
          else
          {
            // this is a non vertex case; ditto
          }
        }
      }
    }
    else if (event_class == 1117 || event_class == 2000)
    {
      operation_success = find_upstream_track(cluster_lines, upstream_track, event_class);
      
cout << "Fit goodness: " << upstream_track.fit_goodness << endl;
      
      if (operation_success)
      {
        // this may need to change if results are not good, particularly for the low momentum runs
        //
        operation_success = assign_clusters_midstream_no_field(cluster_lines, pattern_options, event_class, upstream_track);
      }
      
      if (operation_success)
      {
        // may not be necessary or may need to be done differently with field
        //
        operation_success = find_midstream_tracks_no_field(pattern_options, event_class, upstream_track);
      }
      
      if (operation_success)
      {
        operation_success = propagate_tracks_downstream(pattern_options);
      }
      
      if (operation_success)
      {
        operation_success = assign_clusters_downstream(pattern_options, cluster_lines);
      }

      if (operation_success)
      {
        operation_success = do_global_fit(pattern_options, upstream_track);
      }
    }
  }
    
  return operation_success;
} 

bool track_finding_manager::classify_event(vector<line_3d> cluster_lines[])
{
  int temp_event_class = 0; // all event classes are positive integers
  
  bool is_9001, is_9002, is_9003, is_9004, is_9005;
  // not non-increasing (by new criteria), upstream too few plates hit, upstream too many plates with > 1 hit,
  // upstream too many hits on last plate (implying an early vertex)
  
  int number_of_plates, number_of_upstream_plates, number_of_midstream_plates, number_of_downstream_plates;
  int number_of_groups, number_of_upstream_groups, number_of_midstream_groups, number_of_downstream_groups;
  int number_of_upstream_XY_groups, number_of_midstream_XY_groups, number_of_downstream_XY_groups;
  int midstream_XY_groups[max_number_of_groups_in_tracking_region];
  int dummy_array[max_number_of_groups_in_tracking_region];
  int my_counter;
  
  int zero_cluster_counter_all_regions = 0;
  int zero_cluster_counter_for_region[max_number_of_tracking_regions];
  int single_cluster_counter_all_regions = 0;
  int single_cluster_counter_for_region[max_number_of_tracking_regions];
  int double_cluster_counter_all_regions = 0;
  int double_cluster_counter_for_region[max_number_of_tracking_regions];
  int multiple_cluster_counter_all_regions = 0;
  int multiple_cluster_counter_for_region[max_number_of_tracking_regions];
  
  int min_clusters;
  int dummy_int_1, dummy_int_2;
  
  int group_cluster_multiplicity[max_number_of_groups];
  
  // this may be a better indicator of where the vertex is
  int total_upstream_cluster_lines = 0;
  int total_midstream_cluster_lines = 0;
  int total_downstream_cluster_lines = 0;

  /*
  int number_of_xy_groups_midstream; // must be 2 or 3
  int xy_groups_midstream[max_number_of_groups_in_tracking_region];
  
  int number_of_dd_groups_midstream; // must be 0 or 1
  int dd_groups_midstream[max_number_of_groups_in_tracking_region];

  int number_of_diagonal_plates_midstream; // must be 1 or 2
  int diagonal_plates_midstream[2];
  */
  
  number_of_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates();
  number_of_upstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates();
  number_of_midstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_midstream_plates();
  number_of_downstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_downstream_plates();
  
  number_of_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups();
  number_of_upstream_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_groups();
  number_of_midstream_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_midstream_groups();
  number_of_downstream_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_downstream_groups();
  
  the_global_data_dispatcher->get_partial_geometry()->get_upstream_XY_groups(number_of_upstream_XY_groups, dummy_array);
  the_global_data_dispatcher->get_partial_geometry()->get_midstream_XY_groups(number_of_midstream_XY_groups, midstream_XY_groups);
  the_global_data_dispatcher->get_partial_geometry()->get_downstream_XY_groups(number_of_downstream_XY_groups, dummy_array);
  
  for (int i = 0; i < max_number_of_tracking_regions; i++)
  {
    zero_cluster_counter_for_region[i] = 0;
    single_cluster_counter_for_region[i] = 0;
    double_cluster_counter_for_region[i] = 0;
    multiple_cluster_counter_for_region[i] = 0;
  }
  
  for (int i = 0; i < max_number_of_groups; i++)
  {
    group_cluster_multiplicity[i] = 0;
  }
  
  for (int i = 0; i < number_of_plates; i++)
  {
    if (i < number_of_upstream_plates)
    {
      total_upstream_cluster_lines += cluster_lines[i].size();
    }
    else if (i < number_of_upstream_plates + number_of_midstream_plates)
    {
      total_midstream_cluster_lines += cluster_lines[i].size();
    }
    else
    {
      total_downstream_cluster_lines += cluster_lines[i].size();
    }
    
    if (cluster_lines[i].size() == 0)
    {
      zero_cluster_counter_all_regions++;
      
      if (i < number_of_upstream_plates)
      {
        zero_cluster_counter_for_region[0]++;
      }
      else if (i < number_of_upstream_plates + number_of_midstream_plates)
      {
        zero_cluster_counter_for_region[1]++;
      }
      else
      {
        zero_cluster_counter_for_region[2]++;
      }
    }
    else if (cluster_lines[i].size() == 1)
    {
      single_cluster_counter_all_regions++;
      
      if (i < number_of_upstream_plates)
      {
        single_cluster_counter_for_region[0]++;
      }
      else if (i < number_of_upstream_plates + number_of_midstream_plates)
      {
        single_cluster_counter_for_region[1]++;
      }
      else
      {
        single_cluster_counter_for_region[2]++;
      }
    }
    else if (cluster_lines[i].size() == 2)
    {
      double_cluster_counter_all_regions++;
      
      if (i < number_of_upstream_plates)
      {
        double_cluster_counter_for_region[0]++;
      }
      else if (i < number_of_upstream_plates + number_of_midstream_plates)
      {
        double_cluster_counter_for_region[1]++;
      }
      else
      {
        double_cluster_counter_for_region[2]++;
      }
    }
    else if (cluster_lines[i].size() > 2)
    {
      multiple_cluster_counter_all_regions++;
      
      if (i < number_of_upstream_plates)
      {
        multiple_cluster_counter_for_region[0]++;
      }
      else if (i < number_of_upstream_plates + number_of_midstream_plates)
      {
        multiple_cluster_counter_for_region[1]++;
      }
      else
      {
        multiple_cluster_counter_for_region[2]++;
      }
    }
  }
  
  // ?? this may need to change
  for (int i = 0; i < number_of_groups; i++)
  {
    min_clusters = 1000;
    
    dummy_int_1 = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(i)->number_of_plates_in_group;
    
    for (int j = 0; j < dummy_int_1; j++)
    {
      dummy_int_2 = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(i)->plate_index_list[j];
      
      if (the_global_data_dispatcher->get_partial_geometry()->get_plate(dummy_int_2)->intended_to_measure == "X" ||
          the_global_data_dispatcher->get_partial_geometry()->get_plate(dummy_int_2)->intended_to_measure == "Y")
      {
        if (cluster_lines[dummy_int_2].size() < min_clusters)
        {
          min_clusters = cluster_lines[dummy_int_2].size();
        }
      }
    }
    
    group_cluster_multiplicity[i] = min_clusters;
  }
  
  if (single_cluster_counter_all_regions == number_of_plates)
  {
    temp_event_class = 1111;
  }
  else if ((single_cluster_counter_for_region[0] + single_cluster_counter_for_region[1] + single_cluster_counter_for_region[2]) >=
           number_of_plates - the_global_data_dispatcher->get_reconstruction_options()->max_plates_with_multiple_clusters_for_single_track &&
           single_cluster_counter_for_region[1] >= number_of_midstream_plates - the_global_data_dispatcher->get_reconstruction_options()->max_midstream_plates_with_multiple_clusters_for_single_track &&
           single_cluster_counter_for_region[2] >= number_of_downstream_plates - the_global_data_dispatcher->get_reconstruction_options()->max_downstream_plates_with_multiple_clusters_for_single_track)
  {
    temp_event_class = 1112;
  }
  
  if (zero_cluster_counter_for_region[0] > the_global_data_dispatcher->get_reconstruction_options()->max_upstream_plates_with_zero_clusters)
  {
    is_9002 = true;
  }
  else
  {
    is_9002 = false;
  }
  
  if (double_cluster_counter_for_region[0] + multiple_cluster_counter_for_region[0] > the_global_data_dispatcher->get_reconstruction_options()->max_upstream_plates_with_multiple_clusters)
  {
    is_9003 = true;
  }
  else
  {
    is_9003 = false;
  }
  
  if (cluster_lines[number_of_upstream_plates - 1].size() > 2) // make it a cut; maybe 3 is better
  {
    is_9004 = true;
  }
  else
  {
    is_9004 = false;
  }
  
  if (total_upstream_cluster_lines > 6) // make it a cut
  {
    is_9005 = true;
  }
  else
  {
    is_9005 = false;
  }
  
  if ( (((double)total_midstream_cluster_lines) / ((double)(number_of_midstream_plates))) /
       (((double)total_upstream_cluster_lines) / ((double)(number_of_upstream_plates))) >= 1.6 &&
         cluster_lines[number_of_upstream_plates].size() > 1)
    // the number 1.6 needs to be a cut
    // 1.6 includes 1 1 1 2  2 2 2 2 2
  {
    is_9001 = false;
  }
  else
  {
    if (temp_event_class < 1000)
    {
      is_9001 = true;
    }
    else
    {
      is_9001 = false;
    }
  }
  
  my_counter = 0;
  
  for (int i = 0; i < number_of_midstream_XY_groups; i++)
  {
    if (group_cluster_multiplicity[midstream_XY_groups[i]] >= 2)
    {
      my_counter++;
    }
  }
  
  if (my_counter < 2 && temp_event_class < 1000) // make 2 a cut
  {
    temp_event_class = 1117;
  }
  
  /*
  cout << "^^^^^^^^^^^" << endl;
  for (int i = 0; i < number_of_plates; i++)
  {
    cout << cluster_lines[i].size() << " ";
  }
  cout << endl;
  
  if (is_9001) cout << 9001 << " ";
  if (is_9002) cout << 9002 << " ";
  if (is_9003) cout << 9003 << " ";
  if (is_9004) cout << 9004 << " ";
  if (is_9005) cout << 9005 << " ";
  if (temp_event_class == 1111) cout << 1111 << " " ;
  if (temp_event_class == 1112) cout << 1112 << " " ;
  if (temp_event_class == 1117) cout << 1117 << " " ;
  if (temp_event_class == 2000) cout << 2000 << " " ;
  cout << endl;
  */
  
  if (temp_event_class != 1111 && temp_event_class != 1112 && temp_event_class != 1117) temp_event_class = 2000;
    
  if (is_9001 || is_9002 || is_9003 || is_9004 || is_9005) temp_event_class = 9000;
  
  // cout << "Event class: " << temp_event_class << endl;
  
  
  
  
  
  the_global_data_dispatcher->get_event_reco_output_data_structure()->event_class_by_cluster_multiplicity = temp_event_class;
  
  
  
  return true;
}

bool track_finding_manager::find_upstream_track(vector<line_3d> cluster_lines[], track& upstream_track, int event_class)
{
  bool operation_success;
  
  int plate_index_multiple_clusters;
  int plate_index_zero_clusters;
  
  track temp_track, temp_best_track;
  double best_direction, temp_direction;
  XYZVector z_direction;
  
  string zero_cluster_plate_type;
  line_3d temp_cluster_line;
  
  // general enough ??
  z_direction.SetCoordinates(0., 0., 1.);
  
  int count_zero_cluster_plates = 0;
  int count_multiple_cluster_plates = 0;
  
  int number_of_upstream_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_groups();
  int number_of_upstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates();
  
  upstream_track.cluster_lines.clear();
  upstream_track.plate_index_of_cluster_line.clear();
  
  for (int i = 0; i < number_of_upstream_plates; i++)
  {
    if (cluster_lines[i].size() == 1)
    {
      upstream_track.cluster_lines.push_back(cluster_lines[i].at(0));
      upstream_track.plate_index_of_cluster_line.push_back(i);
    }
    else if (cluster_lines[i].size() == 0)
    {
      count_zero_cluster_plates++;
      plate_index_zero_clusters = i;
    }
    else
    {
      count_multiple_cluster_plates++;
      plate_index_multiple_clusters = i;
    }
  }
  
  // check cuts; !! do properly later
  if (count_zero_cluster_plates + count_multiple_cluster_plates > 1)
  {
    cout << "Problem: too many upstream plates with cluster number different than 1." << endl;
    return false;
  }
  else if (count_zero_cluster_plates + count_multiple_cluster_plates == 1 && event_class == 1111)
  {
    cout << "Problem: event class 1111 with clusters != 1 on an upstream plate." << endl;
    return false;
  }
  
  upstream_track.assumed_start_group = 0;
  upstream_track.assumed_end_group = number_of_upstream_groups - 1;
  
  if (count_zero_cluster_plates + count_multiple_cluster_plates == 0)
  {
    operation_success = fit_track_no_field(upstream_track);
    // calculate the position on target and error (both at z = 0)
//cout << "yes" << endl;
//the_global_data_dispatcher->get_event_reco_output_data_structure()->reconstructed_tracks.push_back(upstream_track);

/*
cout << "uuuuuuuu" << endl;
for (int i = 0; i < 4; i++)
{
  cout << upstream_track.track_plate_intersection_points[i] << endl;
}
*/
    
    the_global_data_dispatcher->get_event_reco_output_data_structure()->the_upstream_track = upstream_track;
    
    return operation_success;
  }
  else if (count_multiple_cluster_plates == 1)
  {
    best_direction = 10.;
    
    for (int i = 0; i < cluster_lines[plate_index_multiple_clusters].size(); i++)
    {
      temp_track.cluster_lines.clear();
      temp_track.plate_index_of_cluster_line.clear();
      
      temp_track = upstream_track;
      temp_track.cluster_lines.push_back(cluster_lines[plate_index_multiple_clusters].at(i));
      temp_track.plate_index_of_cluster_line.push_back(plate_index_multiple_clusters);
      
      temp_track.sort_cluster_lines_by_plate_index();
      operation_success = fit_track_no_field(temp_track);
      // calculate the position on target and error (both at z = 0)
      temp_direction = fabs(sqrt(temp_track.track_line.direction.Unit().Cross(z_direction).Mag2()));
      
      if (temp_direction < best_direction)
      {
        best_direction = temp_direction;
        temp_best_track = temp_track;
      }
    }
    
    // do with cut later !!
    if (best_direction < 0.01)
    {
      upstream_track = temp_best_track;
//cout << "yes" << endl;
//the_global_data_dispatcher->get_event_reco_output_data_structure()->reconstructed_tracks.push_back(upstream_track);
      return true;
/*
cout << "uuuuuuuu" << endl;
for (int i = 0; i < 4; i++)
{
  cout << upstream_track.track_plate_intersection_points[i] << endl;
}
*/
      the_global_data_dispatcher->get_event_reco_output_data_structure()->the_upstream_track = upstream_track;
    }
    else
    {
      return false;
    }
  }
  else if (count_zero_cluster_plates == 1)
  {
    zero_cluster_plate_type = the_global_data_dispatcher->get_partial_geometry()->get_plate(plate_index_zero_clusters)->intended_to_measure;
    
    for (int i = 0; i < number_of_upstream_plates; i++)
    {
      if (the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->intended_to_measure == zero_cluster_plate_type)
      {
        for (int j = 0; j < upstream_track.cluster_lines.size(); j++)
        {
          if (upstream_track.plate_index_of_cluster_line.at(j) == i)
          {
            temp_cluster_line.point = upstream_track.cluster_lines.at(j).point;
            temp_cluster_line.point.SetZ(the_global_data_dispatcher->get_partial_geometry()->get_plate(plate_index_zero_clusters)->position.Z());
            temp_cluster_line.direction = upstream_track.cluster_lines.at(j).direction;
          }
        }
      }
    }
    
    upstream_track.cluster_lines.push_back(temp_cluster_line);
    upstream_track.plate_index_of_cluster_line.push_back(plate_index_zero_clusters);
    upstream_track.sort_cluster_lines_by_plate_index();
    
    operation_success = fit_track_no_field(upstream_track);
    // calculate the position on target and error (both at z = 0)
//cout << "yes" << endl;
//the_global_data_dispatcher->get_event_reco_output_data_structure()->reconstructed_tracks.push_back(upstream_track);

/*
cout << "uuuuuuuu" << endl;
for (int i = 0; i < 4; i++)
{
  cout << upstream_track.track_plate_intersection_points[i] << endl;
}
*/

    the_global_data_dispatcher->get_event_reco_output_data_structure()->the_upstream_track = upstream_track;
    
    return operation_success;
  }
  else
  {
    cout << "A track_finding_manager::find_upstream_track exception ??. This shouldn't happen..." << endl;
    return false;
  }
}

bool track_finding_manager::find_midstream_plus_downstream_track_no_field(std::vector<line_3d> cluster_lines[], track& middownstream_track, int event_class)
{
  bool operation_success;
  
  int number_of_upstream_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_groups();
  int number_of_upstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates();
  int number_of_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups();
  int number_of_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates();
  
  vector<int> unused_plates_list;
  int closest_cluster;
  double closest_cluster_distance, dummy_distance;
  
  middownstream_track.cluster_lines.clear();
  middownstream_track.plate_index_of_cluster_line.clear();
  
  middownstream_track.assumed_start_group = number_of_upstream_groups;
  middownstream_track.assumed_end_group = number_of_groups - 1;
  
  if (event_class == 1111)
  {
    for (int i = number_of_upstream_plates; i < number_of_plates; i++)
    {
      if (cluster_lines[i].size() == 1)
      {
        middownstream_track.cluster_lines.push_back(cluster_lines[i].at(0));
        middownstream_track.plate_index_of_cluster_line.push_back(i);
      }
      else
      {
        cout << "Problem: event class 1111 with clusters != 1 on a midstream/downstream plate." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    operation_success = fit_track_no_field(middownstream_track);
    
    return operation_success;
  }
  else if (event_class == 1112)
  {
    middownstream_track.cluster_lines.clear();
    middownstream_track.plate_index_of_cluster_line.clear();
    unused_plates_list.clear();
    
    int primary_direction_1_plate_counter = 0;
    int primary_direction_2_plate_counter = 0;
    
    for (int i = number_of_upstream_plates; i < number_of_plates; i++)
    {
      if (cluster_lines[i].size() == 1)
      {
        middownstream_track.cluster_lines.push_back(cluster_lines[i].at(0));
        middownstream_track.plate_index_of_cluster_line.push_back(i);
        
        if (the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->intended_to_measure == the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_1_name())
        {
          primary_direction_1_plate_counter++;
        }
        else if (the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->intended_to_measure == the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_2_name())
        {
          primary_direction_2_plate_counter++;
        }
      }
      else if (cluster_lines[i].size() > 1)
      {
        unused_plates_list.push_back(i);
      }
    }
    
    if (primary_direction_1_plate_counter < 2 || primary_direction_2_plate_counter < 2)
    {
      cout << "Track Finding Manager Error ???? See code. Quitting..." << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      operation_success = fit_track_no_field(middownstream_track);
    }
    
    for (int i = 0; i < unused_plates_list.size(); i++)
    {
      closest_cluster = -1;
      closest_cluster_distance = 1000.;
      
      if (cluster_lines[unused_plates_list.at(i)].size() > 1)
      {
        for (int j = 0; j < cluster_lines[unused_plates_list.at(i)].size(); j++)
        {
          if (find_line_to_line_distance(cluster_lines[unused_plates_list.at(i)].at(j), middownstream_track.track_line, dummy_distance))
          {
            if (dummy_distance < closest_cluster_distance)
            {
              closest_cluster = j;
              closest_cluster_distance = dummy_distance;
            }
          }
        }
      }
      
      if (closest_cluster_distance < 0.01) // !! use cut later
      {
        middownstream_track.cluster_lines.push_back(cluster_lines[unused_plates_list.at(i)].at(closest_cluster));
        middownstream_track.plate_index_of_cluster_line.push_back(unused_plates_list.at(i));
      }
    }
    
    middownstream_track.sort_cluster_lines_by_plate_index();
    
    operation_success = fit_track_no_field(middownstream_track);
    
    return operation_success;
  }
}

bool track_finding_manager::find_total_track_no_field(std::vector<line_3d> cluster_lines[], track& total_track, int event_class)
{
  bool operation_success;
  
  vector<int> unused_plates_list;
  int closest_cluster;
  double closest_cluster_distance, dummy_distance;
  
  int number_of_groups = the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups();
  int number_of_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates();
  
  total_track.cluster_lines.clear();
  total_track.plate_index_of_cluster_line.clear();
  
  total_track.assumed_start_group = 0;
  total_track.assumed_end_group = the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups() - 1;
  
  if (event_class == 1111)
  {
    for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); i++)
    {
      if (cluster_lines[i].size() == 1)
      {
        total_track.cluster_lines.push_back(cluster_lines[i].at(0));
        total_track.plate_index_of_cluster_line.push_back(i);
      }
      else
      {
        cout << "Problem: event class 1111 with clusters != 1 on an upstream/midstream/downstream plate." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    operation_success = fit_track_no_field(total_track);
    return operation_success;
  }
  else if (event_class == 1112)
  {
    total_track.cluster_lines.clear();
    total_track.plate_index_of_cluster_line.clear();
    unused_plates_list.clear();
    
    // used to make sure there are at least 2 plates in each primary direction for the initial fit
    int primary_direction_1_plate_counter = 0;
    int primary_direction_2_plate_counter = 0;
    
    for (int i = 0; i < number_of_plates; i++)
    {
      if (cluster_lines[i].size() == 1)
      {
        total_track.cluster_lines.push_back(cluster_lines[i].at(0));
        total_track.plate_index_of_cluster_line.push_back(i);
        
        if (the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->intended_to_measure == the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_1_name())
        {
          primary_direction_1_plate_counter++;
        }
        else if (the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->intended_to_measure == the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_2_name())
        {
          primary_direction_2_plate_counter++;
        }
      }
      else if (cluster_lines[i].size() > 1)
      {
        unused_plates_list.push_back(i);
      }
    }
    
    if (primary_direction_1_plate_counter < 2 || primary_direction_2_plate_counter < 2)
    {
      cout << "Track Finding Manager Error ???? See code. Quitting..." << endl;
      exit(EXIT_FAILURE);
    }
    else
    {
      operation_success = fit_track_no_field(total_track);
    }
    
    for (int i = 0; i < unused_plates_list.size(); i++)
    {
      closest_cluster = -1;
      closest_cluster_distance = 1000.;
      
      if (cluster_lines[unused_plates_list.at(i)].size() > 1)
      {
        for (int j = 0; j < cluster_lines[unused_plates_list.at(i)].size(); j++)
        {
          if (find_line_to_line_distance(cluster_lines[unused_plates_list.at(i)].at(j), total_track.track_line, dummy_distance))
          {
            if (dummy_distance < closest_cluster_distance)
            {
              closest_cluster = j;
              closest_cluster_distance = dummy_distance;
            }
          }
        }
      }
      
      if (closest_cluster_distance < 0.01) // !! use cut later
      {
        total_track.cluster_lines.push_back(cluster_lines[unused_plates_list.at(i)].at(closest_cluster));
        total_track.plate_index_of_cluster_line.push_back(unused_plates_list.at(i));
      }
    }
    
    total_track.sort_cluster_lines_by_plate_index();
    operation_success = fit_track_no_field(total_track);
    
    return operation_success;
  }
}

bool track_finding_manager::fit_track_no_field(track& a_track)
{
  bool fitting_success = true;
  
  double average_x = 0., average_y = 0., average_z = 0.;
  double dx_dz, dy_dz;
  
  string primary_direciton_names[2];
  int primary_direction_plates[2][max_number_of_plates];
  int number_of_primary_direction_plates[2] = {0, 0};
  
  vector<double> initial_values_of_parameters;
  vector<double> errors_in_initial_values_of_parameters;
  
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < max_number_of_plates; j++)
    {
      primary_direction_plates[i][j] = -1;
    }
  }
  
  primary_direciton_names[0] = the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_1_name();
  primary_direciton_names[1] = the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_2_name();
  
  for (int i = 0; i < a_track.cluster_lines.size(); i++)
  {
    if (the_global_data_dispatcher->get_partial_geometry()->get_plate(a_track.plate_index_of_cluster_line.at(i))->intended_to_measure == primary_direciton_names[0])
    {
      primary_direction_plates[0][number_of_primary_direction_plates[0]] = i;
      number_of_primary_direction_plates[0]++;
    }
    else if (the_global_data_dispatcher->get_partial_geometry()->get_plate(a_track.plate_index_of_cluster_line.at(i))->intended_to_measure == primary_direciton_names[1])
    {
      primary_direction_plates[1][number_of_primary_direction_plates[1]] = i;
      number_of_primary_direction_plates[1]++;
    }
  }
  
  if (number_of_primary_direction_plates[0] < 2 || number_of_primary_direction_plates[1] < 2)
  {
    cout << "Too few plates in a given direction to fit a track..." << endl;
    return false;
  }
  
  // also only implemented for X and Y as primary
  for (int i = 0; i < number_of_primary_direction_plates[0]; i++)
  {
    average_z += a_track.cluster_lines.at(primary_direction_plates[0][i]).point.z();
    
    if (primary_direciton_names[0] == "X")
    {
      average_x += a_track.cluster_lines.at(primary_direction_plates[0][i]).point.x();
    }
    else if (primary_direciton_names[0] == "Y")
    {
      average_y += a_track.cluster_lines.at(primary_direction_plates[0][i]).point.y();
    }
    else
    {
      cout << "These primary directions have not been implemented yet" << endl;
      return false;
    }
  }
  
  if (primary_direciton_names[0] == "X")
  {
    average_x /= (double) number_of_primary_direction_plates[0];
  }
  else if (primary_direciton_names[0] == "Y")
  {
    average_y /= (double) number_of_primary_direction_plates[0];
  }
  
  for (int i = 0; i < number_of_primary_direction_plates[1]; i++)
  {
    average_z += a_track.cluster_lines.at(primary_direction_plates[1][i]).point.z();
    
    if (primary_direciton_names[1] == "X")
    {
      average_x += a_track.cluster_lines.at(primary_direction_plates[1][i]).point.x();
    }
    else if (primary_direciton_names[1] == "Y")
    {
      average_y += a_track.cluster_lines.at(primary_direction_plates[1][i]).point.y();
    }
    else
    {
      cout << "These primary directions have not been implemented yet" << endl;
      return false;
    }
  }
  
  if (primary_direciton_names[1] == "X")
  {
    average_x /= (double) number_of_primary_direction_plates[1];
  }
  else if (primary_direciton_names[1] == "Y")
  {
    average_y /= (double) number_of_primary_direction_plates[1];
  }
  
  average_z /= (double)(number_of_primary_direction_plates[0] + number_of_primary_direction_plates[1]);
  
  // this is a particular case; work out the general case later
  //
  if (primary_direciton_names[0] == "X" && primary_direciton_names[1] == "Y")
  {
    dy_dz = (a_track.cluster_lines.at(primary_direction_plates[1][number_of_primary_direction_plates[1] - 1]).point.y() - a_track.cluster_lines.at(primary_direction_plates[1][0]).point.y()) / (a_track.cluster_lines.at(primary_direction_plates[1][number_of_primary_direction_plates[1] - 1]).point.z() - a_track.cluster_lines.at(primary_direction_plates[1][0]).point.z());
    dx_dz = (a_track.cluster_lines.at(primary_direction_plates[0][number_of_primary_direction_plates[0] - 1]).point.x() - a_track.cluster_lines.at(primary_direction_plates[0][0]).point.x()) / (a_track.cluster_lines.at(primary_direction_plates[0][number_of_primary_direction_plates[0] - 1]).point.z() - a_track.cluster_lines.at(primary_direction_plates[0][0]).point.z());
  }
  else if (primary_direciton_names[0] == "Y" && primary_direciton_names[1] == "X")
  {
    dy_dz = (a_track.cluster_lines.at(primary_direction_plates[0][number_of_primary_direction_plates[0] - 1]).point.y() - a_track.cluster_lines.at(primary_direction_plates[0][0]).point.y()) / (a_track.cluster_lines.at(primary_direction_plates[0][number_of_primary_direction_plates[0] - 1]).point.z() - a_track.cluster_lines.at(primary_direction_plates[0][0]).point.z());
    dx_dz = (a_track.cluster_lines.at(primary_direction_plates[1][number_of_primary_direction_plates[1] - 1]).point.x() - a_track.cluster_lines.at(primary_direction_plates[1][0]).point.x()) / (a_track.cluster_lines.at(primary_direction_plates[1][number_of_primary_direction_plates[1] - 1]).point.z() - a_track.cluster_lines.at(primary_direction_plates[1][0]).point.z());
  }
  else
  {
    cout << "This configuration of primary directions has not been implemented yet..." << endl;
    return false;
  }
  
  initial_values_of_parameters.clear();
  errors_in_initial_values_of_parameters.clear();
  
  initial_values_of_parameters.push_back(average_x);
  initial_values_of_parameters.push_back(average_y);
  initial_values_of_parameters.push_back(dx_dz);
  initial_values_of_parameters.push_back(dy_dz);
  
  // need better estimates for these
  errors_in_initial_values_of_parameters.push_back(0.1);
  errors_in_initial_values_of_parameters.push_back(0.1);
  errors_in_initial_values_of_parameters.push_back(0.01);
  errors_in_initial_values_of_parameters.push_back(0.01);
  
  my_straight_line_fcn the_fcn(a_track.cluster_lines, average_z);
  
  VariableMetricMinimizer the_minimizer;
  
  FunctionMinimum min = the_minimizer.Minimize(the_fcn, initial_values_of_parameters, errors_in_initial_values_of_parameters);
  
  a_track.track_line.error_point.SetCoordinates(min.UserState().Error(0), min.UserState().Error(1), 0.);
  // dangerous, may lead to division by 0 !!
  a_track.track_line.error_direction.SetCoordinates(min.UserState().Error(2) / min.UserState().Value(2), min.UserState().Error(3) / min.UserState().Value(3), 0.);
  a_track.track_line.point.SetCoordinates(min.UserState().Value(0), min.UserState().Value(1), average_z);
  a_track.track_line.direction.SetCoordinates(min.UserState().Value(2), min.UserState().Value(3), 1.);
  a_track.track_line.direction = a_track.track_line.direction.Unit();
  
  a_track.fit_goodness = min.Fval();
  
  a_track.find_intersection_points_no_field();
  
  return fitting_success;
}

bool track_finding_manager::assign_clusters_midstream_no_field(vector<line_3d> cluster_lines[], vector<pattern>& pattern_options, int event_class, track upstream_track)
{
  bool operation_success = true;

  // for the 2022 runs there are 2 X, 2 Y and 2 D (1 D1 and 1 D2) plates; the sturcture is general here,
  // allowing more plates, but the implementation only works with 2 plates, so some of the code will need
  // to be modified in runs with higher numbers of plates of each type
  //
  
  // get geometry information
  //
  // ssd plates first
  //
  int number_of_x_plates, number_of_y_plates, number_of_d1_plates, number_of_d2_plates, number_of_d_plates;
  ssd_plate x_plates[max_number_of_plates_in_tracking_region], y_plates[max_number_of_plates_in_tracking_region], d1_plates[max_number_of_plates_in_tracking_region], d2_plates[max_number_of_plates_in_tracking_region], d_plates[max_number_of_plates_in_tracking_region];
  int x_plate_indices[max_number_of_plates_in_tracking_region],  y_plate_indices[max_number_of_plates_in_tracking_region], d1_plate_indices[max_number_of_plates_in_tracking_region], d2_plate_indices[max_number_of_plates_in_tracking_region], d_plate_indices[max_number_of_plates_in_tracking_region];
  
  // now the scales between plates in X and Y so that the vertex has a Z inside the target
  //
  double midstream_X_scales_max[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];
  double midstream_X_scales_min[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];
  double midstream_Y_scales_max[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];
  double midstream_Y_scales_min[max_number_of_single_orientation_mistream_plates * max_number_of_single_orientation_mistream_plates];
  
  // these will hold the indices to the upper and downer X or Y midstream plates
  //
  int x_u_plate, x_d_plate, y_u_plate, y_d_plate;
  // these will hold the number of clusters on each of these plates
  //
  int n_x_u_clusters, n_x_d_clusters, n_y_u_clusters, n_y_d_clusters;
  
  // number of points on the upstream track inside the target - candidate vertices; do with cut later
  //
  const int number_of_divisions = 100;

  // Z coordinates in the target and step for the loop
  //
  double start_z, stop_z, step_z, current_z;
  
  // placeholder for potential vertices
  //
  XYZVector target_point, best_target_point;
  
  // distance accumulators in order to find the best match and temporary such
  //
  double current_total_distance, best_total_distance;
  double temp_best_distance;
  
  // temporary track and line holders
  //
  line_3d temp_line_3d;
  track temp_track_x, temp_track_y, temp_track, temp_combined_track;
  
  // temporary indices to diagonal clusters belonging to a track
  //
  int temp_closest_d_strip[2]; // 2 diagonal plates maximum
  
  // indices of clusters that have a match in the pattern for the current iterated vertex
  //
  vector<int> useful_x_clusters_indices, useful_y_clusters_indices;
  
  // temporary pattern variable
  //
  pattern temp_pattern;
  
  // pattern projection variables for multitrack
  //
  int corresponding_u_points_x[200], corresponding_u_points_y[200]; // this should be max_clusters_per_plate
  double distances[200]; // ditto
  double projection_coordinates[200]; // ditto
  double total_distance_x, total_distance_y;
  double total_track_to_diagonal_strip_distance[max_number_of_tracks][max_number_of_tracks];
  int closest_d_strip[max_number_of_tracks][max_number_of_tracks][max_number_of_diagonal_plates];
  vector<multi_track_selections> the_multi_track_selections;

  
  
  // miscellaneous variables
  //
  double dummy_double, dummy_distance;
  XYZVector temp_point;
  
  
  
  
  
  // begin actual processing
  //
  
  the_global_data_dispatcher->get_partial_geometry()->get_types_of_midstream_plates(number_of_x_plates, x_plates, x_plate_indices, number_of_y_plates, y_plates, y_plate_indices, number_of_d1_plates, d1_plates, d1_plate_indices, number_of_d2_plates, d2_plates, d2_plate_indices, number_of_d_plates, d_plates, d_plate_indices);
  
  the_global_data_dispatcher->get_partial_geometry()->get_midstream_X_scales(midstream_X_scales_max, midstream_X_scales_min);
  the_global_data_dispatcher->get_partial_geometry()->get_midstream_Y_scales(midstream_Y_scales_max, midstream_Y_scales_min);
  
  if (number_of_y_plates != 2 || number_of_x_plates != 2)
  {
    // more than 2 midstream plates of the same direction has to be done later
    //
    cout << "Midstream Y plates or midstream X plates different than 2. Can't do cluster assignment. Quitting..." << endl;
    exit(EXIT_FAILURE);
  }

  y_u_plate = y_plate_indices[0];
  y_d_plate = y_plate_indices[1];
  x_u_plate = x_plate_indices[0];
  x_d_plate = x_plate_indices[1];
  
  n_y_u_clusters = cluster_lines[y_u_plate].size();
  n_y_d_clusters = cluster_lines[y_d_plate].size();
  n_x_u_clusters = cluster_lines[x_u_plate].size();
  n_x_d_clusters = cluster_lines[x_d_plate].size();
  
  start_z = the_global_data_dispatcher->get_partial_geometry()->get_target_min_z();
  stop_z = the_global_data_dispatcher->get_partial_geometry()->get_target_max_z();
  step_z = (stop_z - start_z) / ((double)(number_of_divisions - 1));
  
  best_total_distance = 1000.;
  
  pattern_options.clear();
  
  if (n_y_d_clusters == 1 || n_x_d_clusters == 1)
  {
    // currently keeps only one option; cuts and weights may be done in various ways
    // and so maybe more than one option should be kept
    
    // in all that follows x means 2 or more and n means more than 2
    //
    // cases in this if:
    // (11)(11) (11)(1x) (11)(x1) (1x)(11) (x1)(11) (x1)(1x) (1x)(x1) (x1)(x1)
    
    // in all these cases only a single track can be looked for, so that's what's gonna be done
    //
    
    for (int i = 0; i < n_y_u_clusters; i++)
    {
      for (int j = 0; j < n_y_d_clusters; j++)
      {
    
        temp_track_y.cluster_lines.clear();
        temp_track_y.plate_index_of_cluster_line.clear();
        temp_track_y.cluster_lines.push_back(cluster_lines[y_u_plate].at(i));
        temp_track_y.plate_index_of_cluster_line.push_back(y_u_plate);
        temp_track_y.cluster_lines.push_back(cluster_lines[y_d_plate].at(j));
        temp_track_y.plate_index_of_cluster_line.push_back(y_d_plate);
        
        for (int m = 0; m < n_x_u_clusters; m++)
        {
          for (int n = 0; n < n_x_d_clusters; n++)
          {
            temp_track_x.cluster_lines.clear();
            temp_track_x.plate_index_of_cluster_line.clear();
            temp_track_x.cluster_lines.push_back(cluster_lines[x_u_plate].at(m));
            temp_track_x.plate_index_of_cluster_line.push_back(x_u_plate);
            temp_track_x.cluster_lines.push_back(cluster_lines[x_d_plate].at(n));
            temp_track_x.plate_index_of_cluster_line.push_back(x_d_plate);
            
            combine_3d_line_from_projections(temp_track_x, temp_track_y, temp_line_3d);
            current_total_distance = 0.;
            
            for (int k = 0; k < number_of_d_plates; k++)
            {
              temp_best_distance = 1000.;
              
              for (int l = 0; l < cluster_lines[d_plate_indices[k]].size(); l++)
              {
                find_line_to_line_distance(cluster_lines[d_plate_indices[k]].at(l), temp_line_3d, dummy_distance);
                
                if (dummy_distance < temp_best_distance)
                {
                  temp_best_distance = dummy_distance;
                  temp_closest_d_strip[k] = l;
                }
              }
              
              current_total_distance += temp_best_distance; // use weighting later
            }

            // !! make sure upstream_track.track_line has been calculated
            find_line_to_line_distance(upstream_track.track_line, temp_line_3d, dummy_distance);
            current_total_distance += dummy_distance; // !! add weight here
            
            find_closest_point(temp_line_3d, upstream_track.track_line, temp_point);
            
            if (temp_point.Z() < the_global_data_dispatcher->get_partial_geometry()->get_target_min_z())
            {
              // !! use weights to add
              current_total_distance += fabs(temp_point.Z() - the_global_data_dispatcher->get_partial_geometry()->get_target_min_z());
            }
            else if (temp_point.Z() > the_global_data_dispatcher->get_partial_geometry()->get_target_max_z())
            {
              // !! use weights to add
              current_total_distance += fabs(temp_point.Z() - the_global_data_dispatcher->get_partial_geometry()->get_target_max_z());
            }
            else
            {
              // do nothing here
            }

            if (current_total_distance < best_total_distance)
            {
              best_total_distance = current_total_distance;
              best_target_point = temp_point;
              
              temp_pattern.reset();

              temp_pattern.common_guessed_vertex = best_target_point;
              temp_pattern.number_of_y_lines = 1;
              temp_pattern.number_of_x_lines = 1;
              temp_pattern.y_tracks.push_back(temp_track_y);
              temp_pattern.x_tracks.push_back(temp_track_x);

              find_2d_plane_equation("yz", cluster_lines[y_d_plate].at(j).point, best_target_point, temp_line_3d);
              temp_pattern.equations_y.push_back(temp_line_3d);
              // record pattern; for the moment just the downstream plate
              temp_point.SetCoordinates(0., cluster_lines[y_d_plate].at(j).point.Y(), cluster_lines[y_d_plate].at(j).point.Z());
              temp_pattern.pattern_y_coordinates.push_back(temp_point);

              find_2d_plane_equation("xz", cluster_lines[x_d_plate].at(n).point, best_target_point, temp_line_3d);
              temp_pattern.equations_x.push_back(temp_line_3d);
              // record pattern
              temp_point.SetCoordinates(cluster_lines[x_d_plate].at(n).point.X(), 0., cluster_lines[x_d_plate].at(n).point.Z());
              temp_pattern.pattern_x_coordinates.push_back(temp_point);

              temp_combined_track.cluster_lines.clear();
              temp_combined_track.plate_index_of_cluster_line.clear();
              temp_combined_track.cluster_lines.push_back(cluster_lines[y_u_plate].at(i));
              temp_combined_track.plate_index_of_cluster_line.push_back(y_u_plate);
              temp_combined_track.cluster_lines.push_back(cluster_lines[y_d_plate].at(j));
              temp_combined_track.plate_index_of_cluster_line.push_back(y_d_plate);
              temp_combined_track.cluster_lines.push_back(cluster_lines[x_u_plate].at(m));
              temp_combined_track.plate_index_of_cluster_line.push_back(x_u_plate);
              temp_combined_track.cluster_lines.push_back(cluster_lines[x_d_plate].at(n));
              temp_combined_track.plate_index_of_cluster_line.push_back(x_d_plate);

              for (int k = 0; k < number_of_d_plates; k++)
              {
                temp_combined_track.cluster_lines.push_back(cluster_lines[d_plate_indices[k]].at(temp_closest_d_strip[k]));
                temp_combined_track.plate_index_of_cluster_line.push_back(d_plate_indices[k]);
              }
                                         
              temp_combined_track.sort_cluster_lines_by_plate_index();
              
              temp_pattern.combined_tracks.push_back(temp_combined_track);
              temp_pattern.pattern_quality[0] = best_total_distance;
            }
          }
        }
      }
    }

    // also check if overall goodness cut is obeyed (temp_best_distance < cut), otherwise return false
    if (true) // !! use cuts etc. later
    {
      pattern_options.push_back(temp_pattern);
      operation_success = true;
    }
    else
    {
      operation_success = false;
    }
  }
  else
  {
    // this is the general case; will look for multiple tracks (at least two)
    //
    
    best_total_distance = 1000.;
    
    for (int i = 0; i < number_of_divisions; i++)
    {
      current_z = start_z + step_z * ((double) i);
      temp_point = find_point_on_line_by_z(upstream_track.track_line, current_z);
      
      useful_x_clusters_indices.clear();
      useful_y_clusters_indices.clear();
 
      for (int j = 0; j < n_y_u_clusters; j++)
      {
        dummy_double = temp_point.Y() + (cluster_lines[y_u_plate].at(j).point.Y() - temp_point.Y()) * (y_plates[1].position.Z() - temp_point.Z()) / (y_plates[0].position.Z() - temp_point.Z());
        
        for (int k = 0; k < n_y_d_clusters; k++)
        {
          if (fabs(cluster_lines[y_d_plate].at(k).point.Y() - dummy_double) < 0.02) // do with cut later
          {
            useful_y_clusters_indices.push_back(k);
          }
        }
      }
      
      for (int j = 0; j < n_x_u_clusters; j++)
      {
        dummy_double = temp_point.X() + (cluster_lines[x_u_plate].at(j).point.X() - temp_point.X()) * (x_plates[1].position.Z() - temp_point.Z()) / (x_plates[0].position.Z() - temp_point.Z());
        
        for (int k = 0; k < n_x_d_clusters; k++)
        {
          if (fabs(cluster_lines[x_d_plate].at(k).point.x() - dummy_double) < 0.02) // do with cut later
          {
            useful_x_clusters_indices.push_back(k);
          }
        }
      }

      sort(useful_y_clusters_indices.begin(), useful_y_clusters_indices.end());
      sort(useful_x_clusters_indices.begin(), useful_x_clusters_indices.end());
      
      for (int j = 0; j < useful_y_clusters_indices.size(); j++)
      {
        corresponding_u_points_y[j] = -1;
        distances[j] = 100.;
        
        projection_coordinates[j] = temp_point.Y() + (cluster_lines[y_d_plate].at(useful_y_clusters_indices.at(j)).point.Y() - temp_point.Y()) * (y_plates[0].position.Z() - temp_point.Z()) / (y_plates[1].position.Z() - temp_point.Z());
        
        for (int k = 0; k < n_y_u_clusters; k++)
        {
          //cout << "-->" << fabs(projection_coordinates[j] - cluster_lines[u_plate].at(k).point.Y()) << endl;
          if (fabs(projection_coordinates[j] - cluster_lines[y_u_plate].at(k).point.Y()) < distances[j])
          {
            distances[j] = fabs(projection_coordinates[j] - cluster_lines[y_u_plate].at(k).point.Y());
            corresponding_u_points_y[j] = k;
          }
        }
      }
      
      total_distance_y = 0.;
      
      for (int j = 0; j < useful_y_clusters_indices.size(); j++)
      {
        total_distance_y += distances[j];
      }

      // !! need to weigh by useful_y_clusters_indices.size()
      for (int j = 0; j < useful_x_clusters_indices.size(); j++)
      {
        corresponding_u_points_x[j] = -1;
        distances[j] = 100.;
        
        projection_coordinates[j] = temp_point.X() + (cluster_lines[x_d_plate].at(useful_x_clusters_indices.at(j)).point.X() - temp_point.X()) * (y_plates[0].position.Z() - temp_point.Z()) / (y_plates[1].position.Z() - temp_point.Z());
      
        for (int k = 0; k < n_x_u_clusters; k++)
        {
          //cout << "-->" << fabs(projection_coordinates[j] - cluster_lines[u_plate].at(k).point.Y()) << endl;
          if (fabs(projection_coordinates[j] - cluster_lines[x_u_plate].at(k).point.Y()) < distances[j])
          {
            distances[j] = fabs(projection_coordinates[j] - cluster_lines[x_u_plate].at(k).point.Y());
            corresponding_u_points_x[j] = k;
          }
        }
      }

      total_distance_x = 0.;
      
      for (int j = 0; j < useful_x_clusters_indices.size(); j++)
      {
        total_distance_x += distances[j];
      }
      
      // !! need to weigh by useful_x_clusters_indices.size()

      if (total_distance_x + total_distance_y < best_total_distance)
      {
        best_total_distance = total_distance_x + total_distance_y;
        best_target_point = temp_point;
        
        temp_pattern.reset();
        temp_pattern.common_guessed_vertex = best_target_point;
        
        for (int j = 0; j < useful_y_clusters_indices.size(); j++)
        {
          temp_track.cluster_lines.clear();
          temp_track.plate_index_of_cluster_line.clear();
          temp_track.cluster_lines.push_back(cluster_lines[y_u_plate].at(corresponding_u_points_y[j]));
          temp_track.plate_index_of_cluster_line.push_back(y_u_plate);
          temp_track.cluster_lines.push_back(cluster_lines[y_d_plate].at(useful_y_clusters_indices.at(j)));
          temp_track.plate_index_of_cluster_line.push_back(y_d_plate);
          
          if (n_y_u_clusters < useful_y_clusters_indices.size())
          {
            temp_pattern.number_of_y_lines = n_y_u_clusters;
          }
          else
          {
            temp_pattern.number_of_y_lines = useful_y_clusters_indices.size();
          }
          
          temp_pattern.y_tracks.push_back(temp_track);
          
          // taking the vertex and the downstream plate for the track 2D equation
          find_2d_plane_equation("yz", cluster_lines[y_d_plate].at(useful_y_clusters_indices.at(j)).point, target_point, temp_line_3d);
          temp_pattern.equations_y.push_back(temp_line_3d);
          // record pattern; for the moment just the downstream plate
          temp_point.SetCoordinates(0., cluster_lines[y_d_plate].at(useful_y_clusters_indices.at(j)).point.Y(), cluster_lines[y_d_plate].at(useful_y_clusters_indices.at(j)).point.Z());
          temp_pattern.pattern_y_coordinates.push_back(temp_point);
        }
           
        for (int j = 0; j < useful_x_clusters_indices.size(); j++)
        {
          temp_track.cluster_lines.clear();
          temp_track.plate_index_of_cluster_line.clear();
          temp_track.cluster_lines.push_back(cluster_lines[x_u_plate].at(corresponding_u_points_x[j]));
          temp_track.plate_index_of_cluster_line.push_back(x_u_plate);
          temp_track.cluster_lines.push_back(cluster_lines[x_d_plate].at(useful_x_clusters_indices.at(j)));
          temp_track.plate_index_of_cluster_line.push_back(x_d_plate);
         
          if (n_x_u_clusters < useful_x_clusters_indices.size())
          {
            temp_pattern.number_of_x_lines = n_x_u_clusters;
          }
          else
          {
            temp_pattern.number_of_x_lines = useful_x_clusters_indices.size();
          }
         
          temp_pattern.x_tracks.push_back(temp_track);
         
          // taking the vertex and the downstream plate for the track 2D equation
          find_2d_plane_equation("xz", cluster_lines[x_d_plate].at(useful_x_clusters_indices.at(j)).point, target_point, temp_line_3d);
          temp_pattern.equations_x.push_back(temp_line_3d);
          // record pattern
          temp_point.SetCoordinates(cluster_lines[x_d_plate].at(useful_x_clusters_indices.at(j)).point.X(), 0., cluster_lines[x_d_plate].at(useful_x_clusters_indices.at(j)).point.Z());
          temp_pattern.pattern_x_coordinates.push_back(temp_point);
        }
      }

      for (int i = 0; i < max_number_of_tracks; i++)
      {
        for (int j = 0; j < max_number_of_tracks; j++)
        {
          total_track_to_diagonal_strip_distance[i][j] = 0.;
        }
      }
                                     
      for (int i = 0; i < temp_pattern.x_tracks.size(); i++)
      {
        for (int j = 0; j < temp_pattern.y_tracks.size(); j++)
        {
          combine_3d_line_from_projections(temp_pattern.x_tracks.at(i), temp_pattern.y_tracks.at(j), temp_line_3d);
          
          for (int k = 0; k < number_of_d_plates; k++)
          {
            temp_best_distance = 100.;
            
            for (int kk = 0; kk < cluster_lines[d_plate_indices[k]].size(); kk++)
            {
              find_line_to_line_distance(cluster_lines[d_plate_indices[k]].at(kk), temp_line_3d, dummy_distance);
              
              if (dummy_distance < temp_best_distance)
              {
                temp_best_distance = dummy_distance;
                closest_d_strip[i][j][k] = kk;
              }
            }
            
            total_track_to_diagonal_strip_distance[i][j] += temp_best_distance;
          }
        }
      }

      the_multi_track_selections.clear();
      select_best_multi_track_options(total_track_to_diagonal_strip_distance, temp_pattern.x_tracks.size(), temp_pattern.y_tracks.size(), the_multi_track_selections);
      
      if (the_multi_track_selections.size() == 0)
      {
        cout << "No multitrack selection in a multitrack event." << endl;
        return false;
      }
      
      // ?? how many selections will be returned?; for now use 1
      if (the_multi_track_selections.at(0).selection_quality < 0.1 && best_total_distance < 0.01) // use cuts later; these numbers are totally meaningless now
      {
        for (int i = 0; i < the_multi_track_selections.at(0).number_of_tracks; i++)
        {
          temp_combined_track.cluster_lines.clear();
          temp_combined_track.plate_index_of_cluster_line.clear();

          for (int j = 0; j < temp_pattern.y_tracks.at(the_multi_track_selections.at(0).index_y_track_selections[i]).cluster_lines.size(); j++)
          {
            temp_combined_track.cluster_lines.push_back(temp_pattern.y_tracks.at(the_multi_track_selections.at(0).index_y_track_selections[i]).cluster_lines.at(j));
            temp_combined_track.plate_index_of_cluster_line.push_back(temp_pattern.y_tracks.at(the_multi_track_selections.at(0).index_y_track_selections[i]).plate_index_of_cluster_line.at(j));
          }

          for (int j = 0; j < temp_pattern.x_tracks.at(the_multi_track_selections.at(0).index_x_track_selections[i]).cluster_lines.size(); j++)
          {
            temp_combined_track.cluster_lines.push_back(temp_pattern.x_tracks.at(the_multi_track_selections.at(0).index_x_track_selections[i]).cluster_lines.at(j));
            temp_combined_track.plate_index_of_cluster_line.push_back(temp_pattern.x_tracks.at(the_multi_track_selections.at(0).index_x_track_selections[i]).plate_index_of_cluster_line.at(j));
          }
          
          for (int k = 0; k < number_of_d_plates; k++)
          {
            temp_combined_track.cluster_lines.push_back(cluster_lines[d_plate_indices[k]].at(closest_d_strip[the_multi_track_selections.at(0).index_x_track_selections[i]][the_multi_track_selections.at(0).index_y_track_selections[i]][k]));
            temp_combined_track.plate_index_of_cluster_line.push_back(d_plate_indices[k]);
          }

          temp_combined_track.sort_cluster_lines_by_plate_index();
          
          temp_pattern.combined_tracks.push_back(temp_combined_track);
        }
      }
    } // end loop on target point
    
    // now decide if this is a better option
    //
    if (the_multi_track_selections.at(0).selection_quality < pattern_options.at(0).pattern_quality[1] && best_total_distance < pattern_options.at(0).pattern_quality[0])
    {
      temp_pattern.pattern_quality[0] = best_total_distance;
      temp_pattern.pattern_quality[1] = the_multi_track_selections.at(0).selection_quality;
      pattern_options.clear();
      pattern_options.push_back(temp_pattern);
      // !! or can go for saving more than one option, particularly if they have different number of tracks or some other major difference
    }
                                     
  } // end if on combinations of number of clusters
  
  if (pattern_options.at(0).pattern_quality[0] < 0.001 && pattern_options.at(0).pattern_quality[1] < 0.01) // !! if only one option is selected; these numbers here are completely arbitrary, so change with cuts later
  {
    operation_success = true;
  }
  else
  {
    operation_success = false;
  }
  
  return operation_success;
}

bool track_finding_manager::select_best_multi_track_options(double multi_track_best_distance[][200], int number_of_x_tracks, int number_of_y_tracks, vector<multi_track_selections>& the_multi_track_selections)
{
  bool operation_success;
  
  double total_distance;
  
  multi_track_selections dummy_selection;
  
  if (number_of_x_tracks == 3 && number_of_y_tracks == 3)
  {
    total_distance = multi_track_best_distance[0][0] + multi_track_best_distance[1][1] + multi_track_best_distance[2][2];
// cout << total_distance << endl;
    if (total_distance < 0.3) // rough number; use cut later
    {
      dummy_selection.number_of_tracks = 3;
      dummy_selection.selection_quality = total_distance;
      dummy_selection.index_x_track_selections[0] = 0;
      dummy_selection.index_y_track_selections[0] = 0;
      dummy_selection.index_x_track_selections[1] = 1;
      dummy_selection.index_y_track_selections[1] = 1;
      dummy_selection.index_x_track_selections[2] = 2;
      dummy_selection.index_y_track_selections[2] = 2;
      
      the_multi_track_selections.push_back(dummy_selection);
    }
    
    total_distance = multi_track_best_distance[0][1] + multi_track_best_distance[1][2] + multi_track_best_distance[2][0];
// cout << total_distance << endl;

    if (total_distance < 0.3) // rough number; use cut later
    {
      dummy_selection.number_of_tracks = 3;
      dummy_selection.selection_quality = total_distance;
      dummy_selection.index_x_track_selections[0] = 0;
      dummy_selection.index_y_track_selections[0] = 1;
      dummy_selection.index_x_track_selections[1] = 1;
      dummy_selection.index_y_track_selections[1] = 2;
      dummy_selection.index_x_track_selections[2] = 2;
      dummy_selection.index_y_track_selections[2] = 0;
      
      the_multi_track_selections.push_back(dummy_selection);
    }
    
    total_distance = multi_track_best_distance[0][2] + multi_track_best_distance[1][0] + multi_track_best_distance[2][1];
// cout << total_distance << endl;

    if (total_distance < 0.3) // rough number; use cut later
    {
      dummy_selection.number_of_tracks = 3;
      dummy_selection.selection_quality = total_distance;
      dummy_selection.index_x_track_selections[0] = 0;
      dummy_selection.index_y_track_selections[0] = 2;
      dummy_selection.index_x_track_selections[1] = 1;
      dummy_selection.index_y_track_selections[1] = 0;
      dummy_selection.index_x_track_selections[2] = 2;
      dummy_selection.index_y_track_selections[2] = 1;
      
      the_multi_track_selections.push_back(dummy_selection);
    }
    
    total_distance = multi_track_best_distance[0][2] + multi_track_best_distance[1][1] + multi_track_best_distance[2][0];
// cout << total_distance << endl;

    if (total_distance < 0.3) // rough number; use cut later
    {
      dummy_selection.number_of_tracks = 3;
      dummy_selection.selection_quality = total_distance;
      dummy_selection.index_x_track_selections[0] = 0;
      dummy_selection.index_y_track_selections[0] = 2;
      dummy_selection.index_x_track_selections[1] = 1;
      dummy_selection.index_y_track_selections[1] = 1;
      dummy_selection.index_x_track_selections[2] = 2;
      dummy_selection.index_y_track_selections[2] = 0;
      
      the_multi_track_selections.push_back(dummy_selection);
    }
    
    total_distance = multi_track_best_distance[0][1] + multi_track_best_distance[1][0] + multi_track_best_distance[2][2];
// cout << total_distance << endl;

    if (total_distance < 0.3) // rough number; use cut later
    {
      dummy_selection.number_of_tracks = 3;
      dummy_selection.selection_quality = total_distance;
      dummy_selection.index_x_track_selections[0] = 0;
      dummy_selection.index_y_track_selections[0] = 1;
      dummy_selection.index_x_track_selections[1] = 1;
      dummy_selection.index_y_track_selections[1] = 0;
      dummy_selection.index_x_track_selections[2] = 2;
      dummy_selection.index_y_track_selections[2] = 2;
      
      the_multi_track_selections.push_back(dummy_selection);
    }
    
    total_distance = multi_track_best_distance[0][0] + multi_track_best_distance[1][2] + multi_track_best_distance[2][1];
// cout << total_distance << endl;

    if (total_distance < 0.3) // rough number; use cut later
    {
      dummy_selection.number_of_tracks = 3;
      dummy_selection.selection_quality = total_distance;
      dummy_selection.index_x_track_selections[0] = 0;
      dummy_selection.index_y_track_selections[0] = 0;
      dummy_selection.index_x_track_selections[1] = 1;
      dummy_selection.index_y_track_selections[1] = 2;
      dummy_selection.index_x_track_selections[2] = 2;
      dummy_selection.index_y_track_selections[2] = 1;
      
      the_multi_track_selections.push_back(dummy_selection);
    }
    
    sort(the_multi_track_selections.begin(), the_multi_track_selections.begin(), compare_multi_track_selections);
    
    operation_success = true;
  }
  else
  {
    operation_success = false;
  }
  
// cout << "++++++ " << the_multi_track_selections.size() << endl;
  
  return operation_success;
}

bool track_finding_manager::find_midstream_tracks_no_field(vector<pattern>& pattern_options, int event_class, track upstream_track)
{
  // It makes sense here to start all tracks at the vertex, so the point is common (as it really is),
  // which reduces the number of parameters. Therefore there are 3 (vertex position) + 3 x number_of_tracks
  // parameters in the fit. The vertex needs to have more weight, but doesn't it already?
  //
  
  
  
  
  return true;
}


bool track_finding_manager::propagate_tracks_downstream(vector<pattern>& pattern_options)
{
  // this function will eliminate the most divergent tracks that for sure will not reach downstream;
  // the rule here is if in doubt, leave it, it will not have clusters belonging to it anyway if
  // it doesn't actually propagate
  
  bool operation_success;
  
  plane_3d magnet_front_face;
  XYZVector intersection_point;
  
  double distance_to_axis; // distance between hit point on magnet front face plane and beam axis
  
  magnet_front_face.a = 0.;
  magnet_front_face.b = 0.;
  magnet_front_face.c = 1.;
  magnet_front_face.d = -8.; // magnet's front face is at 8 cm; do with file input later

  
  for (int i = 0; i < pattern_options.size(); i++)
  {
    cout << "tracks:" << pattern_options.at(i).combined_tracks.size() << endl;
    
    for (int j = 0; j < pattern_options.at(i).combined_tracks.size(); j++)
    {
      // cout << pattern_options.at(i).combined_tracks.at(j).track_line.point << " " << pattern_options.at(i).combined_tracks.at(j).track_line.direction << endl;
      
      operation_success = find_line_plate_intersection_point(pattern_options.at(i).combined_tracks.at(j).track_line, magnet_front_face, intersection_point);
      
      if (operation_success)
      {
        distance_to_axis = sqrt(intersection_point.X() * intersection_point.X() + intersection_point.Y() * intersection_point.Y());
        cout << distance_to_axis << " ";
        
        if (distance_to_axis > 3.) // !! fix later
        {
          pattern_options.at(i).combined_tracks.at(j).projects_downstream = false;
        }
        else
        {
          pattern_options.at(i).combined_tracks.at(j).projects_downstream = true;
        }
      }
      else
      {
        cout << "This shouldn't happen! Check the track_finding_manager::propagate_tracks_downstream code." << endl;
        cout << "Quitting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    
    cout << endl;
  }
  
  return operation_success;
}

bool track_finding_manager::assign_clusters_downstream_no_field(vector<pattern>& pattern_options, vector<line_3d> cluster_lines[])
{
  bool operation_success;
  bool geometric_calculation_success;
  
  int number_of_downstream_plates;
  int number_of_upstream_plus_midstream_plates;
  int number_of_total_plates;
  
  double dummy_distance;
  XYZVector intersection_point;
  ssd_plaque* p_dummy_plaque;
  
  number_of_downstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_downstream_plates();
  number_of_upstream_plus_midstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates() + the_global_data_dispatcher->get_partial_geometry()->get_number_of_midstream_plates();;
  number_of_total_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates();
  
  // for now work with pattern option 0
  for (int i = 0; i < pattern_options.at(0).combined_tracks.size(); i++)
  {
    for (int j = number_of_upstream_plus_midstream_plates; j < number_of_total_plates; j++)
    {
      for (int k = 0; k < cluster_lines[j].size(); k++)
      {
        geometric_calculation_success = find_line_to_line_distance(pattern_options.at(0).combined_tracks.at(i).track_line, cluster_lines[j].at(k), dummy_distance);
        
        if (geometric_calculation_success && dummy_distance < 0.02) // !! do with cut; this is 200 microns
        {
          geometric_calculation_success = false;
          
          for (int l = 0; l < the_global_data_dispatcher->get_partial_geometry()->get_plate(j)->number_of_member_plaques; l++)
          {
            p_dummy_plaque = the_global_data_dispatcher->get_partial_geometry()->get_plaque(the_global_data_dispatcher->get_partial_geometry()->get_plate(j)->member_plaques[l]);
            
            if (find_line_plate_intersection_point(pattern_options.at(0).combined_tracks.at(i).track_line, p_dummy_plaque->plaque_plane_equation, intersection_point))
            {
              if (is_intersection_point_inside_plaque(intersection_point, p_dummy_plaque))
              {
                geometric_calculation_success = true;
              }
            }
          }
          
          if (geometric_calculation_success)
          {
            pattern_options.at(0).combined_tracks.at(i).cluster_lines.push_back(cluster_lines[j].at(k));
            pattern_options.at(0).combined_tracks.at(i).plate_index_of_cluster_line.push_back(j);
            
            // remove or don't remove cluster?
          }
        }
      }
    }
    
    pattern_options.at(0).combined_tracks.at(i).sort_cluster_lines_by_plate_index();
  }
  
  operation_success = true;
  
  return operation_success;
}

bool track_finding_manager::do_global_fit_no_field(vector<pattern>& pattern_options, track upstream_track)
{
  bool operation_success = true;
  
  vector<double> initial_values_of_parameters;
  vector<double> errors_in_initial_values_of_parameters;
  
  double x_vertex, y_vertex, z_vertex;
  double dx_dz, dy_dz;
  
  // parameters:
  // 0: vertex.Z()
  // 1, 2 and later groups of two: direction of each track (px/pz and py/pz with dz = 1)
  // the number of parameters depends on the number of tracks
  
  initial_values_of_parameters.clear();
  errors_in_initial_values_of_parameters.clear();

  initial_values_of_parameters.push_back(pattern_options[0].common_guessed_vertex.Z());

  for (int i = 0; i < pattern_options[0].combined_tracks.size(); i++)
  {
    initial_values_of_parameters.push_back(pattern_options[0].combined_tracks.at(i).track_line.direction.X() / pattern_options[0].combined_tracks.at(i).track_line.direction.Z());
    initial_values_of_parameters.push_back(pattern_options[0].combined_tracks.at(i).track_line.direction.Y() / pattern_options[0].combined_tracks.at(i).track_line.direction.Z());
  }
  
  errors_in_initial_values_of_parameters.push_back(0.2); // use input later
  
  for (int i = 0; i < pattern_options[0].combined_tracks.size(); i++)
  {
    errors_in_initial_values_of_parameters.push_back(0.02); // ditto
    errors_in_initial_values_of_parameters.push_back(0.02);
  }
  
  // the zeroth option should be the best; this must have been
  // sorted out by now, but there is still the possibility to use
  // any option - will be figured out with experimentation
  //
  my_global_fit_no_field_fcn the_fcn(pattern_options[0], upstream_track);
  
  VariableMetricMinimizer the_minimizer;
  
  FunctionMinimum min = the_minimizer.Minimize(the_fcn, initial_values_of_parameters, errors_in_initial_values_of_parameters);
  
  cout << min.Fval() << endl;
  
  z_vertex = min.UserState().Value(0);
  
  x_vertex = upstream_track.track_line.point.X() + (z_vertex - upstream_track.track_line.point.Z()) * upstream_track.track_line.direction.X() / upstream_track.track_line.direction.Z();
  y_vertex = upstream_track.track_line.point.Y() + (z_vertex - upstream_track.track_line.point.Z()) * upstream_track.track_line.direction.Y() / upstream_track.track_line.direction.Z();
  
  pattern_options[0].vertex.SetCoordinates(x_vertex, y_vertex, z_vertex);
  
  for (int i = 0; i < pattern_options[0].combined_tracks.size(); i++)
  {
    dx_dz = min.UserState().Value(i * 2 + 1);
    dy_dz = min.UserState().Value(i * 2 + 2);

    pattern_options[0].combined_tracks.at(i).track_line.direction.SetCoordinates(dx_dz, dy_dz, 1.);
    pattern_options[0].combined_tracks.at(i).track_line.direction = pattern_options[0].combined_tracks.at(i).track_line.direction.Unit();
  }
  
  pattern_options[0].pattern_quality[1] = min.Fval();
  
  the_global_data_dispatcher->get_event_reco_output_data_structure()->vertex = pattern_options[0].vertex;
  the_global_data_dispatcher->get_event_reco_output_data_structure()->track_directions.clear();
  for (int i = 0; i < pattern_options[0].combined_tracks.size(); i++)
  {
    the_global_data_dispatcher->get_event_reco_output_data_structure()->track_directions.push_back(pattern_options[0].combined_tracks.at(i).track_line.direction);
  }
  the_global_data_dispatcher->get_event_reco_output_data_structure()->global_fit_quality = pattern_options[0].pattern_quality[1];
  
  for (int i = 0; i < pattern_options.size(); i++)
  {
    for (int j = 0; j < pattern_options.at(i).combined_tracks.size(); j++)
    {
      pattern_options.at(i).combined_tracks.at(j).find_intersection_points_rough();
    }
  }
  
  the_global_data_dispatcher->get_event_reco_output_data_structure()->reconstructed_tracks.clear();
  for (int i = 0; i < pattern_options[0].combined_tracks.size(); i++)
  {
    the_global_data_dispatcher->get_event_reco_output_data_structure()->reconstructed_tracks.push_back(pattern_options[0].combined_tracks.at(i));
  }
  
  return operation_success;
}

bool track_finding_manager::assign_clusters_single_track(vector<line_3d> cluster_lines[], pattern upstream_single_track_pattern, pattern midstream_plus_downstream_single_track_pattern, pattern total_single_track_pattern)
{
  
  
  return true;
}

bool track_finding_manager::do_fit_single_track(track the_track)
{
  vector<double> xRec;
  vector<double> yRec;
  vector<double> txRec;
  vector<double> tyRec;
  vector<double> qdpRec;
  vector<vector<double>> covMat;
  
  double par[5] = {0, 0, 0.0, 0.0, 1/7.};
  double cov[15] = {10,  0,  0,  0,  0,
    10,  0,  0,  0,
    1,  0,  0,
    1,  0,
    1};
  
  double p;
  
  int current_plaque;
  double a, b, c, d;
  plane_3d a_plane_equation;
  
  // need to clear all the local vectors !!
  
  TrackExtrapolation *fExtrap = new TrackExtrapolation(true, 0.05, the_global_data_dispatcher->get_magnetic_field());
  
  // TrackPar trpar(eKisel, fPlanes.at(0).GetZPosition()-fPlanes.at(0).GetSizeZ(), par, cov);
  
  current_plaque = the_global_data_dispatcher->get_partial_geometry()->get_plate(the_track.plate_index_of_cluster_line.at(0))->member_plaques[0]; // !! fix for multiplaque plates
  
  TrackPar trpar(eKisel, the_global_data_dispatcher->get_partial_geometry()->get_plaque(current_plaque)->position.Z() - the_global_data_dispatcher->get_partial_geometry()->get_plaque(current_plaque)->size.Z() / 2., par, cov);
  
  fExtrap->SetTrackPar(trpar);
  
  for( int i = 0; i < the_track.cluster_lines.size(); i++)
  {
    fExtrap->SetRadLength(30390);
    
    current_plaque = the_global_data_dispatcher->get_partial_geometry()->get_plate(the_track.plate_index_of_cluster_line.at(i))->member_plaques[0]; // !! fix for multiplaque plates
    a_plane_equation = the_global_data_dispatcher->get_partial_geometry()->get_plaque(current_plaque)->plaque_plane_equation;
    
    // double a = fPlanes.at(i).GetPlaneNormVector()(0);
    // double b = fPlanes.at(i).GetPlaneNormVector()(1);
    // double c = fPlanes.at(i).GetPlaneNormVector()(2);
    // double d = fPlanes.at(i).GetDPars().at(0);
    
    fExtrap->ExtrapolateToPlane(a_plane_equation.a, a_plane_equation.b, a_plane_equation.c, a_plane_equation.d);
  }
  
  TrackPar &trackPar = fExtrap->GetStopTrackParam();
  
  p = TMath::Abs(1/trackPar.GetPar(4)); // ??
  
  xRec.push_back(trackPar.GetPar(0));
  yRec.push_back(trackPar.GetPar(1));
  txRec.push_back(trackPar.GetPar(2));
  tyRec.push_back(trackPar.GetPar(3));
  qdpRec.push_back(trackPar.GetPar(4));
  
  double *a_cov = trackPar.GetCov();
  
  vector<double> vCov;
  
  for (int i = 0; i < 15; i++)
  {
    vCov.push_back(a_cov[i]);
  }
  
  covMat.push_back(vCov);
  
  return true;
}


bool track_finding_manager::assign_clusters_downstream(vector<pattern>& pattern_options, vector<line_3d> cluster_lines[])
{
  // use a single option for the moment; fix later if necessary !!
  
  
  
  return true;
}


bool track_finding_manager::do_global_fit(vector<pattern>& pattern_options, track upstream_track)
{
  
  
  return true;
}
