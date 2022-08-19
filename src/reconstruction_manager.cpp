#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
//#include <stdlib.h>
//#include <cmath>

#include <TFile.h>
#include <TTree.h>

//#include <TMath.h>
//#include <Math/Vector3D.h>
//#include <Math/GenVector/Rotation3D.h>
//#include <Math/GenVector/RotationX.h>
//#include <Math/GenVector/RotationY.h>
//#include <Math/GenVector/RotationZ.h>

#include "reconstruction_manager.h"
#include "input_output_manager.h"
#include "cluster_manager.h"
#include "transformation_manager.h"
#include "track_finding_manager.h"
//#include "histogram_manager.h"
#include "global_data_dispatcher.h"
#include "partial_geometry.h"

using namespace std;
using namespace ROOT::Math;

reconstruction_manager::reconstruction_manager(input_output_manager* an_input_output_manager)
:the_input_output_manager(an_input_output_manager), the_geometry(0), the_input_root_file(0), the_input_ssd_energy_tree(0), the_input_particle_truth_tree(0), the_input_steps_truth_tree(0), the_input_root_data_structure(0),
    the_run_mc_output_data_structure(0), the_run_reco_output_data_structure(0), the_event_mc_output_data_structure(0), the_event_reco_output_data_structure(0), the_cluster_manager(0), the_transformation_manager(0), the_track_finding_manager(0)
{
  cout << "Created the_reconstruction_manager." << endl;
  
  the_geometry = the_global_data_dispatcher->get_partial_geometry();
  
  the_input_root_data_structure = new input_root_data_structure();
  the_global_data_dispatcher->register_input_root_data_structure(the_input_root_data_structure);
  
  the_run_mc_output_data_structure = new run_mc_output_data_structure();
  the_global_data_dispatcher->register_run_mc_output_data_structure(the_run_mc_output_data_structure);
  
  the_run_reco_output_data_structure = new run_reco_output_data_structure();
  the_global_data_dispatcher->register_run_reco_output_data_structure(the_run_reco_output_data_structure);
  
  the_event_mc_output_data_structure = new event_mc_output_data_structure();
  the_global_data_dispatcher->register_event_mc_output_data_structure(the_event_mc_output_data_structure);
  
  the_event_reco_output_data_structure = new event_reco_output_data_structure();
  the_global_data_dispatcher->register_event_reco_output_data_structure(the_event_reco_output_data_structure);
  
  the_input_output_manager->initialize_root_trees(the_input_root_file, the_input_ssd_energy_tree, the_input_particle_truth_tree, the_input_steps_truth_tree, the_input_root_data_structure);
  
  the_cluster_manager = new cluster_manager();
  
  the_transformation_manager = new transformation_manager(the_geometry);
  
  the_track_finding_manager = new track_finding_manager();
  
  
  
}

void reconstruction_manager::process_events()
{
  // there can be additional clustering after looking at the patterns !!!
  bool operation_success;
  // keep track if it's worth continuing with the analysis of this event
  bool good_event;
  
  int number_of_events_to_process = the_input_output_manager->get_number_of_events_to_process();
  int actual_number_of_plates = the_geometry->get_number_of_plates();
  int actual_number_of_plaques = the_geometry->get_number_of_plaques();

  hit_pair dummy_hit_pair;
  vector<hit_pair> hit_pair_list[max_number_of_plaques];
  vector<double> clusters[max_number_of_plaques];
  
  vector<line_3d> cluster_lines[max_number_of_plates];
  line_3d dummy_line;
  vector<line_3d> plaque_cluster_lines[max_number_of_plaques];
  
  vector<track> the_tracks;

  // ?? temporary
  string dummy_string;
  stringstream dummy_stream;
  
  map<vector<int>, int> cluster_distribution_map;
  vector<int> cluster_distribution_vector;
  vector<pair<string, int>> sorted_cluster_distributions;
  
  map<vector<int>, int> upstream_cluster_distribution_map;
  vector<int> upstream_cluster_distribution_vector;
  vector<pair<string, int>> sorted_upstream_cluster_distributions;
  
  map<vector<int>, int> midstream_cluster_distribution_map;
  vector<int> midstream_cluster_distribution_vector;
  vector<pair<string, int>> sorted_midstream_cluster_distributions;
  
  map<vector<int>, int> downstream_cluster_distribution_map;
  vector<int> downstream_cluster_distribution_vector;
  vector<pair<string, int>> sorted_downstream_cluster_distributions;
  
  int dummy_event_class;
  int counters_event_class[100];

  for (int i = 0; i < 100; i++)
  {
    counters_event_class[i] = 0;
  }
  
  cluster_distribution_map.clear();
  sorted_cluster_distributions.clear();
  
  upstream_cluster_distribution_map.clear();
  sorted_upstream_cluster_distributions.clear();
  
  midstream_cluster_distribution_map.clear();
  sorted_midstream_cluster_distributions.clear();
  
  downstream_cluster_distribution_map.clear();
  sorted_downstream_cluster_distributions.clear();
  // ?? end temporary
  
  // some MC and other outputs
  the_run_mc_output_data_structure->number_of_events = number_of_events_to_process;
  the_run_reco_output_data_structure->number_of_events = number_of_events_to_process;
  the_run_reco_output_data_structure->number_of_plates = actual_number_of_plates;
  the_run_mc_output_data_structure->particle_codes_map.clear();
  
  ofstream out_file("output_vertices.txt");
  
  cout << "Processing events..." << endl;
  
  for (int i = 0; i < number_of_events_to_process; i++)
  //for (int i = 0; i < 1; i++)
  {
    good_event = true;
    
    // ?? temporary
    cluster_distribution_vector.clear();
    upstream_cluster_distribution_vector.clear();
    midstream_cluster_distribution_vector.clear();
    downstream_cluster_distribution_vector.clear();
    // ?? end temporary
    
    for (int j = 0; j < max_number_of_plates; j++)
    {
      cluster_lines[j].clear();
    }
    
    for (int j = 0; j < max_number_of_plaques; j++)
    {
      hit_pair_list[j].clear();
      clusters[j].clear();
      plaque_cluster_lines[j].clear();
    }
    
    //if ((i + 1) % 1000 == 0)
    {
      // cout << i + 1 << " events processed." << endl;
    }
    
    the_input_ssd_energy_tree->GetEntry(i);
    the_input_particle_truth_tree->GetEntry(i);
    the_input_steps_truth_tree->GetEntry(i);
    
    // remove this if geant gives cm
    for (int j = 0; j < the_input_root_data_structure->x_vertex->size(); j++)
    {
      the_input_root_data_structure->x_vertex->at(j) *= mm;
      the_input_root_data_structure->y_vertex->at(j) *= mm;
      the_input_root_data_structure->z_vertex->at(j) *= mm;
    }
    
    // ditto
    for (int j = 0; j < the_input_root_data_structure->step_track_id->size(); j++)
    {
      the_input_root_data_structure->step_x_i->at(j) *= mm;
      the_input_root_data_structure->step_y_i->at(j) *= mm;
      the_input_root_data_structure->step_z_i->at(j) *= mm;
      
      the_input_root_data_structure->step_x_f->at(j) *= mm;
      the_input_root_data_structure->step_y_f->at(j) *= mm;
      the_input_root_data_structure->step_z_f->at(j) *= mm;
    }
    
    calculate_particle_data();
    
    
    for (int j = 0; j < the_input_root_data_structure->plaque_number->size(); j++)
    {
      dummy_hit_pair.strip_number = the_input_root_data_structure->strip_number->at(j);
      dummy_hit_pair.energy_value = the_input_root_data_structure->total_energy_values->at(j);
      
      //if (dummy_hit_pair.energy_value > 0.05)
      {
        hit_pair_list[the_input_root_data_structure->plaque_number->at(j)].push_back(dummy_hit_pair);
      }
    }
    
    operation_success = the_cluster_manager->cluster(hit_pair_list, clusters, actual_number_of_plaques, the_global_data_dispatcher->get_reconstruction_options()->cluster_interpretation_method);
    
    if (!operation_success)
    {
      good_event = false;
    }

    if (good_event)
    {
      operation_success = the_transformation_manager->transform_to_global_tracking_frame(clusters, plaque_cluster_lines);
      
      if (!operation_success)
      {
        good_event = false;
      }
      else
      {
        operation_success = the_transformation_manager->merge_plaques_in_same_plate(plaque_cluster_lines, cluster_lines);
        
        if(!operation_success)
        {
          good_event = false;
        }
        else
        {
          for (int j = 0; j < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); j++)
          {
            the_global_data_dispatcher->get_event_reco_output_data_structure()->cluster_lines[j].clear();
            the_global_data_dispatcher->get_event_reco_output_data_structure()->cluster_lines[j] = cluster_lines[j];
          }
        }
      }
    }
    
    
    
    
    
    

    
    
    the_tracks.clear();

    if (good_event)
    {
      operation_success = the_track_finding_manager->find_tracks(cluster_lines, the_tracks);
      
      if (!operation_success)
      {
        good_event = false;
      }
    }

    // ?? temporary
    //cout << "--------------------" << endl;
    for (int j = 0; j < actual_number_of_plates; j++)
    {
cout << cluster_lines[j].size() << " ";
      
      cluster_distribution_vector.push_back(cluster_lines[j].size());
    }
cout << endl;
    
    for (int j = 0; j < the_geometry->get_number_of_upstream_plates(); j++)
    {
      upstream_cluster_distribution_vector.push_back(cluster_lines[j].size());
    }
    
    for (int j = the_geometry->get_number_of_upstream_plates(); j < (the_geometry->get_number_of_upstream_plates() + the_geometry->get_number_of_midstream_plates()); j++)
    {
      midstream_cluster_distribution_vector.push_back(cluster_lines[j].size());
    }
    
    for (int j = (the_geometry->get_number_of_upstream_plates() + the_geometry->get_number_of_midstream_plates()); j < actual_number_of_plates; j++)
    {
      downstream_cluster_distribution_vector.push_back(cluster_lines[j].size());
    }
    
    //if (the_global_data_dispatcher->get_event_reco_output_data_structure()->group_cluster_multiplicity[2] == 2 && the_global_data_dispatcher->get_event_reco_output_data_structure()->group_cluster_multiplicity[3] == 2)
    {
      cluster_distribution_map[cluster_distribution_vector]++;
      
      upstream_cluster_distribution_map[upstream_cluster_distribution_vector]++;
      
      midstream_cluster_distribution_map[midstream_cluster_distribution_vector]++;
      
      downstream_cluster_distribution_map[downstream_cluster_distribution_vector]++;
    }
    // ?? end temporary
    
    the_input_output_manager->record_event();
    
    dummy_event_class = the_global_data_dispatcher->get_event_reco_output_data_structure()->event_class_by_cluster_multiplicity;
    
    if (dummy_event_class == 1111)
    {
      counters_event_class[0]++;
    }
    else if (dummy_event_class == 1112)
    {
      counters_event_class[1]++;
    }
    else if (dummy_event_class == 1117)
    {
      counters_event_class[2]++;
    }
    else if (dummy_event_class == 2000)
    {
      counters_event_class[3]++;
    }
    else if (dummy_event_class == 9000)
    {
      counters_event_class[4]++;
    }
    
// cout << dummy_event_class << endl;
    XYZVector simple_vertex;

    /*
    if (dummy_event_class == 1111 && the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.size() > 1)
    {
      cout << "--------------------" << endl;
      cout << "HAS VERTEX" << endl;
    }
    */
    
// cout << "mc_vertex" << the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex << endl;

    //if (dummy_event_class == 1111 && the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.size() > 1)
    if (dummy_event_class == 1111)
    {
      //
      /*
      cout << "begin 1111 event" << endl;
      
      cout << "mc vertex: " << endl;
      cout << the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex << endl;
      
      cout << "tracks:" << endl;
      for (int k = 0; k < the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.size(); k++)
      {
        cout << "track id: " << the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(k).track_id << endl;
        cout << "particle code: " << the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(k).particle_code << endl;
        cout << "mc hits:" << endl;
        for (int l = 0; l < the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(k).position.size(); l++)
        {
          cout << the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(k).position.at(l) << " ";
        }
        cout << endl;
      }
      
      cout << "reco silicon strip cluster lines: " << endl;
      for (int k = 0; k < the_geometry->get_number_of_plates(); k++)
      {
        for (int l = 0; l < the_global_data_dispatcher->get_event_reco_output_data_structure()->cluster_lines[k].size(); l++)
        {
          cout << the_global_data_dispatcher->get_event_reco_output_data_structure()->cluster_lines[k].at(l).point << " ";
        }
      }
      cout << endl;
      */
      
      double x[2], z_x[2], y[2], z_y[2];
      XYZVector point;
      XYZVector direction;
      line_3d line_upstream, line_midstream;

      x[0] = cluster_lines[1].at(0).point.X();
      x[1] = cluster_lines[3].at(0).point.X();
      y[0] = cluster_lines[0].at(0).point.Y();
      y[1] = cluster_lines[2].at(0).point.Y();
      z_x[0] = cluster_lines[1].at(0).point.Z();
      z_x[1] = cluster_lines[3].at(0).point.Z();
      z_y[0] = cluster_lines[0].at(0).point.Z();
      z_y[1] = cluster_lines[2].at(0).point.Z();

      find_3d_line_from_projections(x, z_x, y, z_y, point, direction);
      
      line_upstream.point = point;
      line_upstream.direction = direction;

      x[0] = cluster_lines[5].at(0).point.X();
      x[1] = cluster_lines[8].at(0).point.X();
      y[0] = cluster_lines[4].at(0).point.Y();
      y[1] = cluster_lines[7].at(0).point.Y();
      z_x[0] = cluster_lines[5].at(0).point.Z();
      z_x[1] = cluster_lines[8].at(0).point.Z();
      z_y[0] = cluster_lines[4].at(0).point.Z();
      z_y[1] = cluster_lines[7].at(0).point.Z();
      
      find_3d_line_from_projections(x, z_x, y, z_y, point, direction);
      
      line_midstream.point = point;
      line_midstream.direction = direction;
      
      find_closest_point(line_upstream, line_midstream, simple_vertex);
      
      //cout << "Simple vertex: " << simple_vertex << endl;
      // cout << "Line directions: " << line_upstream.direction << " " << line_midstream.direction << endl;
      
      line_upstream.direction = line_upstream.direction.Unit();
      line_midstream.direction = line_midstream.direction.Unit();
      
      //cout << "from lines: " << line_upstream.direction.Dot(line_midstream.direction) << endl;
    }
    
    
    /*
    if (the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex.X() != -10. & !(simple_vertex.X() == 0. && simple_vertex.Y() == 0. && simple_vertex.Z() == 0.))
    {
      out_file << the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex.X() << " "
               << the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex.Y() << " "
               << the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex.Z() << " "
               << simple_vertex.X() << " "
               << simple_vertex.Y() << " "
               << simple_vertex.Z() << endl;
    }
    */
    
    out_file << the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex.X() << " "
             << the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex.Y() << " "
             << the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex.Z() << " "
             << the_global_data_dispatcher->get_event_reco_output_data_structure()->guessed_vertex.X() << " "
             << the_global_data_dispatcher->get_event_reco_output_data_structure()->guessed_vertex.Y() << " "
             << the_global_data_dispatcher->get_event_reco_output_data_structure()->guessed_vertex.Z() << endl;

    
    
  } // end of event loop
  
  
  
  out_file.close();
  
  the_input_output_manager->record_run();
  
  // ?? temporary
  for (map<vector<int>, int>::iterator it = cluster_distribution_map.begin(); it != cluster_distribution_map.end(); ++it)
  {
    dummy_string = "";
    dummy_stream.str("");
    
    //if (it->second >= 5)
    {
      for (int j = 0; j < actual_number_of_plates; j++)
      {
        dummy_stream << it->first.at(j) << " ";
      }
      dummy_string = dummy_stream.str();
      
      sorted_cluster_distributions.push_back(make_pair(dummy_string, it->second));
    }
  }
  
  sort(sorted_cluster_distributions.begin(), sorted_cluster_distributions.end(), compare_map);
  
  
  cout << "All plates: " << endl;
  for (int i = 0; i < sorted_cluster_distributions.size(); i++)
  {
    cout << sorted_cluster_distributions.at(i).first << "=> " << sorted_cluster_distributions.at(i).second << endl;
  }

  /*
  for (map<vector<int>, int>::iterator it = upstream_cluster_distribution_map.begin(); it != upstream_cluster_distribution_map.end(); ++it)
  {
    dummy_string = "";
    dummy_stream.str("");
    
    //if (it->second >= 5)
    {
      for (int j = 0; j < the_geometry->get_number_of_upstream_plates(); j++)
      {
        dummy_stream << it->first.at(j) << " ";
      }
      dummy_string = dummy_stream.str();
      
      sorted_upstream_cluster_distributions.push_back(make_pair(dummy_string, it->second));
    }
  }
  
  sort(sorted_upstream_cluster_distributions.begin(), sorted_upstream_cluster_distributions.end(), compare_map);
  
  
 
  cout << "Upstream plates: " << endl;
  
  for (int i = 0; i < sorted_upstream_cluster_distributions.size(); i++)
  {
    cout << sorted_upstream_cluster_distributions.at(i).first << "=> " << sorted_upstream_cluster_distributions.at(i).second << endl;
  }
  
  cout << the_geometry->get_number_of_upstream_plates() << " "
  << the_geometry->get_number_of_midstream_plates() << " "
  << the_geometry->get_number_of_downstream_plates() << endl;
  
  for (map<vector<int>, int>::iterator it = midstream_cluster_distribution_map.begin(); it != midstream_cluster_distribution_map.end(); ++it)
  {
    dummy_string = "";
    dummy_stream.str("");
    
    //if (it->second >= 5)
    {
      for (int j = 0; j < the_geometry->get_number_of_midstream_plates(); j++)
      {
        dummy_stream << it->first.at(j) << " ";
      }
      dummy_string = dummy_stream.str();
      
      sorted_midstream_cluster_distributions.push_back(make_pair(dummy_string, it->second));
    }
  }
  
  sort(sorted_midstream_cluster_distributions.begin(), sorted_midstream_cluster_distributions.end(), compare_map);
  
  cout << "Midstream plates: " << endl;
  
  for (int i = 0; i < sorted_midstream_cluster_distributions.size(); i++)
  {
    cout << sorted_midstream_cluster_distributions.at(i).first << "=> " << sorted_midstream_cluster_distributions.at(i).second << endl;
  }
  
  for (map<vector<int>, int>::iterator it = downstream_cluster_distribution_map.begin(); it != downstream_cluster_distribution_map.end(); ++it)
  {
    dummy_string = "";
    dummy_stream.str("");
    
    //if (it->second >= 5)
    {
      for (int j = 0; j < the_geometry->get_number_of_downstream_plates(); j++)
      {
        dummy_stream << it->first.at(j) << " ";
      }
      dummy_string = dummy_stream.str();
      
      sorted_downstream_cluster_distributions.push_back(make_pair(dummy_string, it->second));
    }
  }
  
  sort(sorted_downstream_cluster_distributions.begin(), sorted_downstream_cluster_distributions.end(), compare_map);
  
  cout << "Downstream plates: " << endl;
  
  for (int i = 0; i < sorted_downstream_cluster_distributions.size(); i++)
  {
    cout << sorted_downstream_cluster_distributions.at(i).first << "=> " << sorted_downstream_cluster_distributions.at(i).second << endl;
  }
  // ?? end temporary
  */
  
  cout << "Event class -> Number of events" << endl;
  cout << "1111 -> " << counters_event_class[0] << endl;
  cout << "1112 -> " << counters_event_class[1] << endl;
  cout << "1117 -> " << counters_event_class[2] << endl;
  cout << "2000 -> " << counters_event_class[3] << endl;
  cout << "9000 -> " << counters_event_class[4] << endl;


  int temp_sum = 0;
  
  for (int i = 0; i < 100; i++)
  {
    temp_sum += counters_event_class[i];
  }
  cout << "Total cases: " << temp_sum << endl;
  
  return;
}

void reconstruction_manager::calculate_particle_data()
{
  bool same_vertex = true;
  double dummy_double;
  
  int dummy_id;
  
  charged_track dummy_track;
  map<int, int> track_id_to_vector_element_map;
  
  double x_sum, y_sum, z_sum, de_sum;
  int counter_steps_in_strip;
  XYZVector dummy_vector;
  bool first_counted = true;
  int previous_step_index;
  
  mc_step dummy_step;
  
  // calculate MC vertex
  
  int dummy_counter = 0;
  
  for (int i = 0; i < the_global_data_dispatcher->get_input_root_data_structure()->track_id->size(); i++)
  {
    dummy_double = -1000.;
    
    if (the_global_data_dispatcher->get_input_root_data_structure()->parent_track_id->at(i) == 1 && the_global_data_dispatcher->get_input_root_data_structure()->ekin_vertex->at(i) > 500.)
    {
      dummy_counter++;
      /*
       cout << the_input_root_data_structure->track_id->at(i) << " "
       << the_input_root_data_structure->particle_code->at(i) << " "
       << the_input_root_data_structure->ekin_vertex->at(i) << " "
       << the_input_root_data_structure->x_vertex->at(i) << " "
       << the_input_root_data_structure->y_vertex->at(i) << " "
       << the_input_root_data_structure->z_vertex->at(i) << endl;
       */
      
      if (dummy_double < -100.)
      {
        dummy_double = the_global_data_dispatcher->get_input_root_data_structure()->z_vertex->at(i);
        the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex.SetCoordinates(the_global_data_dispatcher->get_input_root_data_structure()->x_vertex->at(i), the_global_data_dispatcher->get_input_root_data_structure()->y_vertex->at(i), the_global_data_dispatcher->get_input_root_data_structure()->z_vertex->at(i));
      }
      else
      {
        if (the_global_data_dispatcher->get_input_root_data_structure()->z_vertex->at(i) != dummy_double)
        {
          // cout << dummy_double << " " << the_input_root_data_structure->z_vertex->at(j) << endl;
          // cout << "aaaaa" << endl;
          same_vertex = false;
          
        }
      }
    }
  }
  if (!same_vertex)
  {
    the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex.SetCoordinates(-10., -10., -10.);
  }
  
  if (dummy_counter == 0) the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex.SetCoordinates(-10., -10., -10.);
  
  // now calculate MC track things
  
  the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.clear();
  the_global_data_dispatcher->get_event_mc_output_data_structure()->the_selected_steps.clear();
  track_id_to_vector_element_map.clear();
    
  for (int i = 0; i < the_global_data_dispatcher->get_input_root_data_structure()->track_id->size(); i++)
  {
    the_global_data_dispatcher->get_run_mc_output_data_structure()->particle_codes_map[the_global_data_dispatcher->get_input_root_data_structure()->particle_code->at(i)]++;
    
    if (the_global_data_dispatcher->get_input_root_data_structure()->ekin_vertex->at(i) > the_global_data_dispatcher->get_reconstruction_options()->cut_min_energy_to_be_high_energy_charged_particle)
    {
      dummy_id = the_global_data_dispatcher->get_input_root_data_structure()->particle_code->at(i);
      
      if (abs(dummy_id) == 11 || abs(dummy_id) == 13 || abs(dummy_id) == 211 || abs(dummy_id) == 321 || abs(dummy_id) == 2212 || abs(dummy_id) == 3112 || abs(dummy_id) == 3222 || abs(dummy_id) == 3312 || abs(dummy_id) == 3334)
      {
        dummy_track.track_id = the_global_data_dispatcher->get_input_root_data_structure()->track_id->at(i);
        dummy_track.particle_code = dummy_id;
        dummy_track.ekin_vertex = the_global_data_dispatcher->get_input_root_data_structure()->ekin_vertex->at(i);
        the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.push_back(dummy_track);
        track_id_to_vector_element_map[the_global_data_dispatcher->get_input_root_data_structure()->track_id->at(i)] = the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.size() - 1;
        
        /*
         cout << the_input_root_data_structure->track_id->at(i) << " " << dummy_id << " " << the_input_root_data_structure->ekin_vertex->at(i) << endl;
         */
        
        /*
         if (the_input_root_data_structure->visited_plates->at(i).size() > 0)
         {
         cout << "Visited plates: ";
         for (int j = 0; j < the_input_root_data_structure->visited_plates->at(i).size(); j++)
         {
         cout << the_input_root_data_structure->visited_plates->at(i).at(j) << " ";
         }
         cout << endl;
         }
         */
      }
    }
  }
  
  /*
   cout << "Track id to vector element map size:" << track_id_to_vector_element_map.size() << endl;
   for (map<int, int>::iterator it = track_id_to_vector_element_map.begin(); it != track_id_to_vector_element_map.end(); it++)
   {
   cout << it->first << " " << it->second << endl;
   }
   */
  
  counter_steps_in_strip = 0;
  x_sum = 0.;
  y_sum = 0.;
  z_sum = 0.;
  de_sum = 0.;
  
  for (int i = 0; i < the_global_data_dispatcher->get_input_root_data_structure()->step_track_id->size(); i++)
  {
    if (the_global_data_dispatcher->get_input_root_data_structure()->step_volume_name->at(i) == "physical_ssd_strip" && track_id_to_vector_element_map.find(the_global_data_dispatcher->get_input_root_data_structure()->step_track_id->at(i)) != track_id_to_vector_element_map.end())
    {
      if (first_counted)
      {
        x_sum  = (the_global_data_dispatcher->get_input_root_data_structure()->step_x_f->at(i) + the_global_data_dispatcher->get_input_root_data_structure()->step_x_i->at(i)) / 2.;
        y_sum  = (the_global_data_dispatcher->get_input_root_data_structure()->step_y_f->at(i) + the_global_data_dispatcher->get_input_root_data_structure()->step_y_i->at(i)) / 2.;
        z_sum  = (the_global_data_dispatcher->get_input_root_data_structure()->step_z_f->at(i) + the_global_data_dispatcher->get_input_root_data_structure()->step_z_i->at(i)) / 2.;
        de_sum =  the_global_data_dispatcher->get_input_root_data_structure()->step_final_energy->at(i) - the_global_data_dispatcher->get_input_root_data_structure()->step_initial_energy->at(i);
        counter_steps_in_strip++;
        first_counted = false;
        previous_step_index = i;
      }
      else if (the_global_data_dispatcher->get_input_root_data_structure()->step_track_id->at(i) == the_global_data_dispatcher->get_input_root_data_structure()->step_track_id->at(i - 1) && the_global_data_dispatcher->get_input_root_data_structure()->step_plaque_id->at(i) == the_global_data_dispatcher->get_input_root_data_structure()->step_plaque_id->at(i - 1))
      {
        x_sum  += (the_global_data_dispatcher->get_input_root_data_structure()->step_x_f->at(i) + the_global_data_dispatcher->get_input_root_data_structure()->step_x_i->at(i)) / 2.;
        y_sum  += (the_global_data_dispatcher->get_input_root_data_structure()->step_y_f->at(i) + the_global_data_dispatcher->get_input_root_data_structure()->step_y_i->at(i)) / 2.;
        z_sum  += (the_global_data_dispatcher->get_input_root_data_structure()->step_z_f->at(i) + the_global_data_dispatcher->get_input_root_data_structure()->step_z_i->at(i)) / 2.;
        de_sum +=  the_global_data_dispatcher->get_input_root_data_structure()->step_initial_energy->at(i) - the_global_data_dispatcher->get_input_root_data_structure()->step_final_energy->at(i);
        counter_steps_in_strip++;
        previous_step_index = i;
      }
      else
      {
        dummy_id = track_id_to_vector_element_map[the_global_data_dispatcher->get_input_root_data_structure()->step_track_id->at(previous_step_index)];
        
        the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(dummy_id).step_plate_index.push_back(the_global_data_dispatcher->get_input_root_data_structure()->step_plaque_id->at(previous_step_index));
        dummy_vector.SetCoordinates(x_sum / ((double) counter_steps_in_strip), y_sum / ((double) counter_steps_in_strip), z_sum / ((double) counter_steps_in_strip));
        the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(dummy_id).position.push_back(dummy_vector);
        
        the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(dummy_id).energy_lost.push_back(de_sum);
        
        counter_steps_in_strip = 0;
        
        x_sum  = (the_global_data_dispatcher->get_input_root_data_structure()->step_x_f->at(i) + the_global_data_dispatcher->get_input_root_data_structure()->step_x_i->at(i)) / 2.;
        y_sum  = (the_global_data_dispatcher->get_input_root_data_structure()->step_y_f->at(i) + the_global_data_dispatcher->get_input_root_data_structure()->step_y_i->at(i)) / 2.;
        z_sum  = (the_global_data_dispatcher->get_input_root_data_structure()->step_z_f->at(i) + the_global_data_dispatcher->get_input_root_data_structure()->step_z_i->at(i)) / 2.;
        de_sum =  the_global_data_dispatcher->get_input_root_data_structure()->step_initial_energy->at(i) - the_global_data_dispatcher->get_input_root_data_structure()->step_final_energy->at(i);
        counter_steps_in_strip++;
        previous_step_index = i;
      }
      
      dummy_step.step_index = i;
      dummy_step.step_track_id = the_global_data_dispatcher->get_input_root_data_structure()->step_track_id->at(i);
      dummy_step.step_volume_name = the_global_data_dispatcher->get_input_root_data_structure()->step_volume_name->at(i);
      dummy_step.step_plate_id = the_global_data_dispatcher->get_input_root_data_structure()->step_plaque_id->at(i);
      dummy_step.initial_position.SetCoordinates(the_global_data_dispatcher->get_input_root_data_structure()->step_x_i->at(i), the_global_data_dispatcher->get_input_root_data_structure()->step_y_i->at(i), the_global_data_dispatcher->get_input_root_data_structure()->step_z_i->at(i));
      dummy_step.final_position.SetCoordinates(the_global_data_dispatcher->get_input_root_data_structure()->step_x_f->at(i), the_global_data_dispatcher->get_input_root_data_structure()->step_y_f->at(i), the_global_data_dispatcher->get_input_root_data_structure()->step_z_f->at(i));
      dummy_step.initial_energy = the_global_data_dispatcher->get_input_root_data_structure()->step_initial_energy->at(i);
      dummy_step.final_energy = the_global_data_dispatcher->get_input_root_data_structure()->step_final_energy->at(i);
      dummy_step.step_process_name = the_global_data_dispatcher->get_input_root_data_structure()->step_process_name->at(i);
      the_global_data_dispatcher->get_event_mc_output_data_structure()->the_selected_steps.push_back(dummy_step);
    }
  }
  
  dummy_id = track_id_to_vector_element_map[the_global_data_dispatcher->get_input_root_data_structure()->step_track_id->at(previous_step_index)];
  
  /*
   cout << dummy_id << " " << the_input_root_data_structure->step_track_id->at(previous_step_index) << endl;
   cout << counter_steps_in_strip << " " << the_input_root_data_structure->step_track_id->at(previous_step_index) << " " << the_input_root_data_structure->step_plate_id->at(previous_step_index) << " " << x_sum << endl;
   */
  
  the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(dummy_id).step_plate_index.push_back(the_global_data_dispatcher->get_input_root_data_structure()->step_plaque_id->at(previous_step_index));
  dummy_vector.SetCoordinates(x_sum / ((double) counter_steps_in_strip), y_sum / ((double) counter_steps_in_strip), z_sum / ((double) counter_steps_in_strip));
  the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(dummy_id).position.push_back(dummy_vector);
  
  the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(dummy_id).energy_lost.push_back(de_sum);
}

