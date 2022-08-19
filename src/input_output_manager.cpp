#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sstream>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <Math/Vector3D.h>

#include "input_output_manager.h"
#include "jsmn.h"
#include "global_data_dispatcher.h"
#include "partial_geometry.h"
#include "matej_minimum_data_structures_and_functions.h"

using namespace std;
using namespace ROOT::Math;

input_output_manager::input_output_manager()
:output_message(""), the_input_root_file(0), the_input_ssd_energy_tree(0),
  the_input_particle_truth_tree(0), the_input_steps_truth_tree(0),
  number_of_events_to_process(0), the_reconstruction_options(0),
  the_algorithm(0), the_magnetic_field(0),
  the_output_root_file(0), the_run_output_tree(0),
  the_event_mc_output_tree(0), the_event_reco_output_tree(0)
{
  if (DEBUG_PRINT_LEVEL > 0)
  {
    cout << "Created the_input_output_manager." << endl;
  }
  
  the_reconstruction_options = new reconstruction_options();
  the_global_data_dispatcher->register_reconstruction_options(the_reconstruction_options);
  
  the_algorithm = new reconstruction_algorithm();
  the_global_data_dispatcher->register_reconstruction_algorithm(the_algorithm);

  // Why are these being initialized here????
  the_reconstruction_options->max_upstream_plates_with_multiple_clusters = 0;
  the_reconstruction_options->max_clusters_on_upstream_plates_with_multiple_clusters = 100; // arb large number
  the_reconstruction_options->max_midstream_plus_downstream_plates_with_multiple_clusters = 0;
  the_reconstruction_options->max_midstream_plus_downstream_plates_with_zero_clusters = 0;
  the_reconstruction_options->max_midstream_plus_downstream_plates_with_zero_or_multiple_clusters = 0;
  the_reconstruction_options->max_midstream_plus_downstream_plates_with_less_than_2_clusters = 0;
  
  // temporary output file
  out_file.open("output_debug_22_4.txt");
}

input_output_manager::~input_output_manager()
{
  the_input_root_file->Close();
  //the_output_root_file->Close();
  
  if (the_input_root_file) delete the_input_root_file;
  if (the_input_ssd_energy_tree) delete the_input_ssd_energy_tree;
  if (the_input_particle_truth_tree) delete the_input_particle_truth_tree;
  if (the_input_steps_truth_tree) delete the_input_steps_truth_tree;
  if (the_output_root_file) delete the_output_root_file;
  if (the_run_output_tree) delete the_run_output_tree;
  if (the_event_mc_output_tree) delete the_event_mc_output_tree;
  if (the_event_reco_output_tree) delete the_event_reco_output_tree;
  
  out_file.close();
  
  if (DEBUG_PRINT_LEVEL > 0)
  {
    cout << "Deleted the_input_output_manager." << endl;
  }
}

bool input_output_manager::initialize_files(int argc, char** argv)
{
  int number_of_tokens;
  string current_token;
  
  bool initialization_succeeded = true;
  
  string command_input_file_name = "";
  magnetic_field_file_name = "";
  geometry_file_name = "";
  input_root_file_name = "";
  
  string dummy_string;
  
  struct stat buffer;
  
  ifstream command_input_file_stream;
  
  string line;
  string card;
  string string_card_value;
  int int_card_value;
  
  if (DEBUG_PRINT_LEVEL > 0)
  {
    cout << "File initialization started in the_input_output_manager." << endl;
  }
  
  if (argc == 1)
  {
    output_message += "Too few command line arguments. Usage:\n";
    output_message += argv[0];
    output_message += " command_input_file.txt\n";
    initialization_succeeded = false;
  }
  else if (argc > 2)
  {
    output_message += "Too many command line arguments. Usage:\n";
    output_message += argv[0];
    output_message += " command_input_file.txt\n";
    
    initialization_succeeded = false;
  }
  else
  {
    command_input_file_name = argv[1];
    
    if (stat(command_input_file_name.c_str(), &buffer) != 0)
    {
      output_message += "The command_input_file does not exist: ";
      output_message += command_input_file_name;
      output_message += "\n";
      
      initialization_succeeded = false;
    }
    else
    {
      output_message += "Using command input file name: ";
      output_message += command_input_file_name;
      output_message += "\n";
      
      command_input_file_stream.open(command_input_file_name.c_str(), ifstream::in);
      
      while (getline(command_input_file_stream, line))
      {
        bool empty_line = false;
        
        if (line.empty())
        {
          empty_line = true;
        }
        else
        {
          bool line_of_spaces = true;
          
          for (int i = 0; i < line.length(); i++)
          {
            if (line[i] != ' ')
            {
              line_of_spaces = false;
            }
          }
          
          if (!line_of_spaces) empty_line = false;
        }
        
        if (!empty_line)
        {
          if (line[0] != '#')
          {
            istringstream temp_string_stream(line);
            temp_string_stream >> card;
            
            if (card == "use_magnetic_field")
            {
              temp_string_stream >> dummy_string;
              
              if (dummy_string == "true")
              {
                the_algorithm->use_magnetic_field = true;
              }
              else if (dummy_string == "false")
              {
                the_algorithm->use_magnetic_field = false;
              }
              else
              {
                output_message += "Uknown option for use_magnetic_field: ";
                output_message += dummy_string;
                output_message += "\n";
                output_message += "Will use default value: false";
                the_algorithm->use_magnetic_field = false;
              }
            }
            else if (card == "magnetic_field_file_name")
            {
              temp_string_stream >> magnetic_field_file_name;
            }
            else if (card == "geometry_file_name")
            {
              temp_string_stream >> geometry_file_name;
            }
            else if (card == "input_root_file_name")
            {
              temp_string_stream >> input_root_file_name;
            }
            else if (card == "input_root_ssd_energy_tree_name")
            {
              temp_string_stream >> input_root_ssd_energy_tree_name;
            }
            else if (card == "input_root_particle_truth_tree_name")
            {
              temp_string_stream >> input_root_particle_truth_tree_name;
            }
            else if (card == "input_root_steps_truth_tree_name")
            {
              temp_string_stream >> input_root_steps_truth_tree_name;
            }
            else if (card == "output_root_file_name")
            {
              temp_string_stream >> output_root_file_name;
            }
            else if (card == "number_of_events_to_process")
            {
              temp_string_stream >> dummy_string;
              
              if (dummy_string == "all")
              {
                number_of_events_to_process = 0;
              }
              else
              {
                number_of_events_to_process = stoi(dummy_string);
              }
            }
            else if (card == "cluster_interpretation_method")
            {
              temp_string_stream >> dummy_string;
              
              if (dummy_string == "weighted")
              {
                the_reconstruction_options->cluster_interpretation_method = dummy_string;
              }
              else
              {
                output_message += "Unknown cluster interpretation method: ";
                output_message += dummy_string;
                output_message += "\n Will use the default weighted method...\n";
                the_reconstruction_options->cluster_interpretation_method = "weighted";
              }
            }
            else if (card == "maximum_angle_between_intended_and_actual_strip_direction")
            {
              temp_string_stream >> dummy_string;
              
              the_reconstruction_options->maximum_cosangle_between_intended_and_actual_strip_direction = cos(stod(dummy_string) * to_radians);
            }
            else if (card == "max_upstream_plates_with_multiple_clusters")
            {
              temp_string_stream >> dummy_string;

              the_reconstruction_options->max_upstream_plates_with_multiple_clusters = stoi(dummy_string);
            }
            else if (card == "max_clusters_on_upstream_plates_with_multiple_clusters")
            {
              temp_string_stream >> dummy_string;
              the_reconstruction_options->max_clusters_on_upstream_plates_with_multiple_clusters = stoi(dummy_string);
            }
            else if (card == "max_midstream_plus_downstream_plates_with_multiple_clusters")
            {
              temp_string_stream >> dummy_string;
              the_reconstruction_options->max_midstream_plus_downstream_plates_with_multiple_clusters = stoi(dummy_string);
            }
            else if (card == "max_midstream_plus_downstream_plates_with_zero_clusters")
            {
              temp_string_stream >> dummy_string;
              the_reconstruction_options->max_midstream_plus_downstream_plates_with_zero_clusters = stoi(dummy_string);
            }
            else if (card == "max_midstream_plus_downstream_plates_with_zero_or_multiple_clusters")
            {
              temp_string_stream >> dummy_string;
              the_reconstruction_options->max_midstream_plus_downstream_plates_with_zero_or_multiple_clusters = stoi(dummy_string);
            }
            else if (card == "max_midstream_plus_downstream_plates_with_less_than_2_clusters")
            {
              temp_string_stream >> dummy_string;
              the_reconstruction_options->max_midstream_plus_downstream_plates_with_less_than_2_clusters = stoi(dummy_string);
            }
            else if (card == "max_plates_with_multiple_clusters_for_single_track")
            {
              temp_string_stream >> dummy_string;
              the_reconstruction_options->max_plates_with_multiple_clusters_for_single_track = stoi(dummy_string);
            }
            else if (card == "max_upstream_plates_with_multiple_clusters_for_single_track")
            {
              temp_string_stream >> dummy_string;
              the_reconstruction_options->max_upstream_plates_with_multiple_clusters_for_single_track = stoi(dummy_string);
            }
            else if (card == "max_midstream_plates_with_multiple_clusters_for_single_track")
            {
              temp_string_stream >> dummy_string;
              the_reconstruction_options->max_midstream_plates_with_multiple_clusters_for_single_track = stoi(dummy_string);
            }
            else if (card == "max_downstream_plates_with_multiple_clusters_for_single_track")
            {
              temp_string_stream >> dummy_string;
              the_reconstruction_options->max_downstream_plates_with_multiple_clusters_for_single_track = stoi(dummy_string);
            }
            else if (card == "max_upstream_plates_with_zero_clusters")
            {
              temp_string_stream >> dummy_string;
              the_reconstruction_options->max_upstream_plates_with_zero_clusters = stoi(dummy_string);
            }
            else if (card == "cut_max_distance_to_accept_in_track")
            {
              temp_string_stream >> dummy_string;
              the_reconstruction_options->cut_max_distance_to_accept_in_track = stod(dummy_string);
            }
            else if (card == "guessed_vertex_weight")
            {
              temp_string_stream >> dummy_string;
              the_reconstruction_options->guessed_vertex_weight = stod(dummy_string);
            }
            else if (card == "cut_max_distance_to_count_to_vertex")
            {
              temp_string_stream >> dummy_string;
              the_reconstruction_options->cut_max_distance_to_count_to_vertex = stod(dummy_string);
            }
            else if (card == "cut_min_energy_to_be_high_energy_charged_particle")
            {
              temp_string_stream >> dummy_string;
              the_reconstruction_options->cut_min_energy_to_be_high_energy_charged_particle = stod(dummy_string);
            }
            else
            {
              output_message += "Unknown input item: ";
              output_message += card;
              output_message += "\n";
              output_message += "Ignored...\n";
            }
          }
        }
      }
    }
    
    output_message += "Use magnetic field: ";
    
    if (the_algorithm->use_magnetic_field)
    {
      output_message += "yes.\n";
    }
    else
    {
      output_message += "no.\n";
    }
    
    if (the_algorithm->use_magnetic_field)
    {
      dummy_string = "config/";
      dummy_string += magnetic_field_file_name;
      magnetic_field_file_name = dummy_string;
      output_message += "Magnetic field file: ";
      output_message += magnetic_field_file_name;
      output_message += "\n";
    }
    
    dummy_string = "config/";
    dummy_string += geometry_file_name;
    geometry_file_name = dummy_string;
    output_message += "Geometry file: ";
    output_message += geometry_file_name;
    output_message += "\n";
    
    dummy_string = "../data/";
    dummy_string += input_root_file_name;
    input_root_file_name = dummy_string;
    output_message += "Input root file: ";
    output_message += input_root_file_name;
    output_message += "\n";
    
    std::ostringstream oss;
    oss << number_of_events_to_process;
 
    output_message += "Requested to process ";
    if (number_of_events_to_process == 0)
    {
      output_message += "all";
    }
    else
    {
      output_message += oss.str();
    }
    output_message += " events.\n";
  }
  
  if (initialization_succeeded)
  {
    if (the_algorithm->use_magnetic_field)
    {
      if (stat(magnetic_field_file_name.c_str(), &buffer) != 0)
      {
        output_message += "Magnetic field file does not exist: ";
        output_message += magnetic_field_file_name;
        output_message += "\n";
        initialization_succeeded = false;
        
        cout << "Error 2001: Magnetic field file does not exist: " << magnetic_field_file_name << endl;
      }
    }
    
    if (stat(geometry_file_name.c_str(), &buffer) != 0)
    {
      output_message += "Geometry file does not exist: ";
      output_message += geometry_file_name;
      output_message += "\n";
      initialization_succeeded = false;
      
      cout << "Error 2002: Geometry file does not exist: " << geometry_file_name << endl;
    }
    
    if (stat(input_root_file_name.c_str(), &buffer) != 0)
    {
      output_message += "Input root file does not exist: ";
      output_message += input_root_file_name;
      output_message += "\n";
      initialization_succeeded = false;
      
      cout << "Error 2003: Input root file does not exist: " << input_root_file_name << endl;
    }
  }
  
  output_message += "Output histograms and ntuples will be written to file: ";
  output_message += output_root_file_name;
  output_message +="\n";
  
  if (DEBUG_PRINT_LEVEL > 1)
  {
    cout << output_message;
  }
  
  return initialization_succeeded;
}

bool input_output_manager::initialize_geometry(partial_geometry* the_geometry)
{
  bool initialization_succeeded = true;

  int number_of_tokens;
  int p1_token;
  int p1_size;
  string current_token;
  
  string dummy_string;
  
  ssd_plaque temp_plaque;
  target temp_target;
  double temp_x, temp_y, temp_z;
  XYZVector temp_vector;
  
  output_message = "";
  
  ifstream geometry_file(geometry_file_name.c_str(), ifstream::in);
  
  string geometry_file_string((std::istreambuf_iterator<char>(geometry_file)), (std::istreambuf_iterator<char>()));
  
  geometry_file.close();
  
  const char* json_string_1 = geometry_file_string.c_str();
  
  jsmn_parser p;
  jsmntok_t t[2048];
  
  jsmn_init(&p);
  number_of_tokens = jsmn_parse(&p, json_string_1, strlen(json_string_1), t, sizeof(t) / sizeof(t[0]));
  
  if (number_of_tokens < 0)
  {
    output_message += "Failed to parse geometry json file: ";
    output_message += geometry_file_name;
    output_message += "\n";
    
    initialization_succeeded = false;
    
    cout << "Error 2011: Failed to parse geometry json file: " << geometry_file_name << endl;
  }
  else
  {
    output_message += "Number of tokens in geometry json file: ";
    std::ostringstream oss;
    oss << number_of_tokens;
    output_message += oss.str();
    output_message += "\n";
  }
  
  if (number_of_tokens < 1 || t[0].type != JSMN_OBJECT)
  {
    output_message += "Object expected as the root of the geometry json file.\n";
    
    initialization_succeeded = false;
    
    cout << "Error 2012: Object expected as the root of the geometry json file." << endl;
  }
  
  int i = 1;
  
  while (i < number_of_tokens - 1)
  {
    // i++;
    
    current_token = geometry_file_string.substr(t[i].start, t[i].end - t[i].start);
    
    cout << "->->->-> " << i << " " << t[i].type << " " << t[i].size << " " << current_token << endl;
    
    if (current_token == "primary_direction_1")
    {
      current_token = geometry_file_string.substr(t[i + 2].start, t[i + 2].end - t[i + 2].start);
      temp_x = stod(current_token);
      current_token = geometry_file_string.substr(t[i + 3].start, t[i + 3].end - t[i + 3].start);
      temp_y = stod(current_token);
      current_token = geometry_file_string.substr(t[i + 4].start, t[i + 4].end - t[i + 4].start);
      temp_z = stod(current_token);
      temp_vector.SetCoordinates(temp_x, temp_y, temp_z);
      temp_vector = temp_vector.Unit();
      the_geometry->set_primary_direction_1(temp_vector);
      i += 5;
    }
    else if (current_token == "primary_direction_2")
    {
      current_token = geometry_file_string.substr(t[i + 2].start, t[i + 2].end - t[i + 2].start);
      temp_x = stod(current_token);
      current_token = geometry_file_string.substr(t[i + 3].start, t[i + 3].end - t[i + 3].start);
      temp_y = stod(current_token);
      current_token = geometry_file_string.substr(t[i + 4].start, t[i + 4].end - t[i + 4].start);
      temp_z = stod(current_token);
      temp_vector.SetCoordinates(temp_x, temp_y, temp_z);
      temp_vector = temp_vector.Unit();
      the_geometry->set_primary_direction_2(temp_vector);
      i += 5;
    }
    else if (current_token == "secondary_direction_1")
    {
      current_token = geometry_file_string.substr(t[i + 2].start, t[i + 2].end - t[i + 2].start);
      temp_x = stod(current_token);
      current_token = geometry_file_string.substr(t[i + 3].start, t[i + 3].end - t[i + 3].start);
      temp_y = stod(current_token);
      current_token = geometry_file_string.substr(t[i + 4].start, t[i + 4].end - t[i + 4].start);
      temp_z = stod(current_token);
      temp_vector.SetCoordinates(temp_x, temp_y, temp_z);
      temp_vector = temp_vector.Unit();
      the_geometry->set_secondary_direction_1(temp_vector);
      i += 5;
    }
    else if (current_token == "secondary_direction_2")
    {
      current_token = geometry_file_string.substr(t[i + 2].start, t[i + 2].end - t[i + 2].start);
      temp_x = stod(current_token);
      current_token = geometry_file_string.substr(t[i + 3].start, t[i + 3].end - t[i + 3].start);
      temp_y = stod(current_token);
      current_token = geometry_file_string.substr(t[i + 4].start, t[i + 4].end - t[i + 4].start);
      temp_z = stod(current_token);
      temp_vector.SetCoordinates(temp_x, temp_y, temp_z);
      temp_vector = temp_vector.Unit();
      the_geometry->set_secondary_direction_2(temp_vector);
      i += 5;
    }
    else if (current_token == "primary_direction_1_name")
    {
      dummy_string = geometry_file_string.substr(t[i + 1].start, t[i + 1].end - t[i + 1].start);
      i += 2;
      the_geometry->set_primary_direction_1_name(dummy_string);
    }
    else if (current_token == "primary_direction_2_name")
    {
      dummy_string = geometry_file_string.substr(t[i + 1].start, t[i + 1].end - t[i + 1].start);
      i += 2;
      the_geometry->set_primary_direction_2_name(dummy_string);
    }
    else if (current_token == "secondary_direction_1_name")
    {
      dummy_string = geometry_file_string.substr(t[i + 1].start, t[i + 1].end - t[i + 1].start);
      i += 2;
      the_geometry->set_secondary_direction_1_name(dummy_string);
    }
    else if (current_token == "secondary_direction_2_name")
    {
      dummy_string = geometry_file_string.substr(t[i + 1].start, t[i + 1].end - t[i + 1].start);
      i += 2;
      the_geometry->set_secondary_direction_2_name(dummy_string);
    }
    else if (current_token == "magnetic_field_general_direction")
    {
      current_token = geometry_file_string.substr(t[i + 2].start, t[i + 2].end - t[i + 2].start);
      temp_x = stod(current_token);
      current_token = geometry_file_string.substr(t[i + 3].start, t[i + 3].end - t[i + 3].start);
      temp_y = stod(current_token);
      current_token = geometry_file_string.substr(t[i + 4].start, t[i + 4].end - t[i + 4].start);
      temp_z = stod(current_token);
      temp_vector.SetCoordinates(temp_x, temp_y, temp_z);
      temp_vector = temp_vector.Unit();
      the_geometry->set_magnetic_field_general_direction(temp_vector);
      i += 5;
    }
    else if (current_token == "magnetic_field_general_direction_name")
    {
      dummy_string = geometry_file_string.substr(t[i + 1].start, t[i + 1].end - t[i + 1].start);
      i += 2;
      the_geometry->set_magnetic_field_general_direction_name(dummy_string);
    }
    else if (current_token == "tracking_region_0_name")
    {
      dummy_string = geometry_file_string.substr(t[i + 1].start, t[i + 1].end - t[i + 1].start);
      i += 2;
      the_geometry->set_tracking_region_names(0, dummy_string);
    }
    else if (current_token == "tracking_region_1_name")
    {
      dummy_string = geometry_file_string.substr(t[i + 1].start, t[i + 1].end - t[i + 1].start);
      i += 2;
      the_geometry->set_tracking_region_names(1, dummy_string);
    }
    else if (current_token == "tracking_region_2_name")
    {
      dummy_string = geometry_file_string.substr(t[i + 1].start, t[i + 1].end - t[i + 1].start);
      i += 2;
      the_geometry->set_tracking_region_names(2, dummy_string);
    }
    else if (current_token == "tracking_region_3_name")
    {
      dummy_string = geometry_file_string.substr(t[i + 1].start, t[i + 1].end - t[i + 1].start);
      i += 2;
      the_geometry->set_tracking_region_names(3, dummy_string);
    }
    else if (current_token == "plaque")
    {
      p1_token = i + 1;
      p1_size = t[p1_token].size;
      p1_token++;
    
      for (int j = 0; j < p1_size; j++)
      {
        current_token = geometry_file_string.substr(t[p1_token].start, t[p1_token].end - t[p1_token].start);
     
        if (current_token == "id")
        {
          current_token = geometry_file_string.substr(t[p1_token + 1].start, t[p1_token + 1].end - t[p1_token + 1].start);
          temp_plaque.plaque_id = stoi(current_token);
          p1_token += 2;
        }
        else if (current_token == "tracking_region")
        {
          current_token = geometry_file_string.substr(t[p1_token + 1].start, t[p1_token + 1].end - t[p1_token + 1].start);
          temp_plaque.tracking_region = stoi(current_token);
          p1_token += 2;
        }
        else if (current_token == "plate")
        {
          current_token = geometry_file_string.substr(t[p1_token + 1].start, t[p1_token + 1].end - t[p1_token + 1].start);
          temp_plaque.plate = stoi(current_token);
          p1_token += 2;
        }
        else if (current_token == "group")
        {
          current_token = geometry_file_string.substr(t[p1_token + 1].start, t[p1_token + 1].end - t[p1_token + 1].start);
          temp_plaque.group = stoi(current_token);
          p1_token += 2;
        }
        else if (current_token == "intended_to_measure")
        {
          current_token = geometry_file_string.substr(t[p1_token + 1].start, t[p1_token + 1].end - t[p1_token + 1].start);
          temp_plaque.type = current_token;
          p1_token += 2;
        }
        else if (current_token == "efficiency")
        {
          current_token = geometry_file_string.substr(t[p1_token + 1].start, t[p1_token + 1].end - t[p1_token + 1].start);
          temp_plaque.efficiency = stod(current_token);
          p1_token += 2;
        }
        else if (current_token == "position")
        {
          current_token = geometry_file_string.substr(t[p1_token + 2].start, t[p1_token + 2].end - t[p1_token + 2].start);
          temp_x = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 3].start, t[p1_token + 3].end - t[p1_token + 3].start);
          temp_y = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 4].start, t[p1_token + 4].end - t[p1_token + 4].start);
          temp_z = stod(current_token);
          temp_plaque.position.SetXYZ(temp_x, temp_y, temp_z);
          p1_token += 5;
        }
        else if (current_token == "rotation")
        {
          current_token = geometry_file_string.substr(t[p1_token + 2].start, t[p1_token + 2].end - t[p1_token + 2].start);
          temp_x = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 3].start, t[p1_token + 3].end - t[p1_token + 3].start);
          temp_y = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 4].start, t[p1_token + 4].end - t[p1_token + 4].start);
          temp_z = stod(current_token);
          temp_plaque.rotation.SetXYZ(temp_x, temp_y, temp_z);
          p1_token += 5;
        }
        else if (current_token == "size")
        {
          current_token = geometry_file_string.substr(t[p1_token + 2].start, t[p1_token + 2].end - t[p1_token + 2].start);
          temp_x = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 3].start, t[p1_token + 3].end - t[p1_token + 3].start);
          temp_y = stod(current_token);
          temp_plaque.size.SetXYZ(temp_x, temp_y, 0.);
          p1_token += 4;
        }
        else if (current_token == "number_of_strips")
        {
          current_token = geometry_file_string.substr(t[p1_token + 1].start, t[p1_token + 1].end - t[p1_token + 1].start);
          temp_plaque.number_of_strips = stoi(current_token);
          p1_token += 2;
        }
      }
      
      the_geometry->add_plaque(temp_plaque);
      i = p1_token;
      
//cout << ">>>> " << geometry_file_string.substr(t[i].start, t[i].end - t[i].start) << endl;
    }
    else if (current_token == "target")
    {
      p1_size = t[i + 1].size;
      i += 2;
      
      for (int j = 0; j < p1_size; j++)
      {
        current_token = geometry_file_string.substr(t[i].start, t[i].end - t[i].start);
        
        if (current_token == "position")
        {
          current_token = geometry_file_string.substr(t[i + 2].start, t[i + 2].end - t[i + 2].start);
          temp_x = stod(current_token);
          current_token = geometry_file_string.substr(t[i + 3].start, t[i + 3].end - t[i + 3].start);
          temp_y = stod(current_token);
          current_token = geometry_file_string.substr(t[i + 4].start, t[i + 4].end - t[i + 4].start);
          temp_z = stod(current_token);
          temp_target.position.SetXYZ(temp_x, temp_y, temp_z);
          i += 5;
        }
        if (current_token == "rotation")
        {
          current_token = geometry_file_string.substr(t[i + 2].start, t[i + 2].end - t[i + 2].start);
          temp_x = stod(current_token);
          current_token = geometry_file_string.substr(t[i + 3].start, t[i + 3].end - t[i + 3].start);
          temp_y = stod(current_token);
          current_token = geometry_file_string.substr(t[i + 4].start, t[i + 4].end - t[i + 4].start);
          temp_z = stod(current_token);
          temp_target.rotation.SetXYZ(temp_x, temp_y, temp_z);
          i += 5;
        }
        else if (current_token == "size")
        {
          current_token = geometry_file_string.substr(t[i + 2].start, t[i + 2].end - t[i + 2].start);
          temp_x = stod(current_token);
          current_token = geometry_file_string.substr(t[i + 3].start, t[i + 3].end - t[i + 3].start);
          temp_y = stod(current_token);
          current_token = geometry_file_string.substr(t[i + 4].start, t[i + 4].end - t[i + 4].start);
          temp_z = stod(current_token);
          temp_target.size.SetXYZ(temp_x, temp_y, temp_z);
          i += 5;
        }
        else if (current_token == "material")
        {
          temp_target.material = geometry_file_string.substr(t[p1_token + 1].start, t[p1_token + 1].end - t[p1_token + 1].start);
          i += 2;
        }

        the_geometry->add_target(temp_target);
      }
    }
  }
  
  // create lists of x, y and d plates, and calculate number of upstream, midstream and downstream plates
  the_geometry->classify_plates();
  
  if (the_algorithm->use_magnetic_field)
  {
    the_magnetic_field = new MagneticField(magnetic_field_file_name);
    the_global_data_dispatcher->register_magnetic_field(the_magnetic_field);
  }
  
  if (DEBUG_PRINT_LEVEL > 1)
  {
    cout << output_message;
  }
  
  return initialization_succeeded;
}

bool input_output_manager::initialize_root_trees(TFile*& an_input_root_file, TTree*& an_input_ssd_energy_tree, TTree*& an_input_particle_truth_tree, TTree*& an_input_steps_truth_tree, input_root_data_structure*& the_input_root_data_structure)
{
  bool initialization_succeeded = true;
  
  output_message = "";
  
  the_input_root_file = new TFile(input_root_file_name.c_str(), "READ");
  the_input_ssd_energy_tree = (TTree*) the_input_root_file->Get(input_root_ssd_energy_tree_name.c_str());
  the_input_particle_truth_tree = (TTree*) the_input_root_file->Get(input_root_particle_truth_tree_name.c_str());
  the_input_steps_truth_tree = (TTree*) the_input_root_file->Get(input_root_steps_truth_tree_name.c_str());

  the_output_root_file = new TFile(output_root_file_name.c_str(),"RECREATE");
  the_run_output_tree = new TTree(run_output_tree_name.c_str(), run_output_tree_name.c_str());
  the_event_mc_output_tree = new TTree(event_mc_output_tree_name.c_str(), event_mc_output_tree_name.c_str());
  the_event_reco_output_tree  = new TTree(event_reco_output_tree_name.c_str(), event_reco_output_tree_name.c_str());
  
  if (!the_input_ssd_energy_tree)
  {
    output_message += "Root tree not found: ";
    output_message += input_root_ssd_energy_tree_name;
    output_message += "\n";
    initialization_succeeded = false;
    
    cout << "Error 2051: Root tree not found: " << input_root_ssd_energy_tree_name << endl;
  }
  
  if (!the_input_particle_truth_tree)
  {
    output_message += "Root tree not found: ";
    output_message += input_root_particle_truth_tree_name;
    output_message += "\n";
    initialization_succeeded = false;
    
    cout << "Error 2052: Root tree not found: " << input_root_particle_truth_tree_name << endl;
  }

  if (!the_input_steps_truth_tree)
  {
    output_message += "Root tree not found: ";
    output_message += input_root_steps_truth_tree_name;
    output_message += "\n";
    initialization_succeeded = false;
    
    cout << "Error 2053: Root tree not found: " << input_root_steps_truth_tree_name << endl;
  }
  
  if (!(the_input_ssd_energy_tree->GetEntries() == the_input_particle_truth_tree->GetEntries() && the_input_ssd_energy_tree->GetEntries() == the_input_steps_truth_tree->GetEntries()))
  {
    output_message += "The root trees have different number of events: ";
    output_message += the_input_ssd_energy_tree->GetEntries();
    output_message += " ";
    output_message += the_input_particle_truth_tree->GetEntries();
    output_message += " ";
    output_message += the_input_steps_truth_tree->GetEntries();
    output_message += "\n Quitting...\n";
    initialization_succeeded = false;
    
    cout << "Error 2053: The root trees have different number of events: " << the_input_ssd_energy_tree->GetEntries() << " " << the_input_particle_truth_tree->GetEntries() << " " << the_input_steps_truth_tree->GetEntries() << endl;
  }
  
  if (initialization_succeeded)
  {
    if (number_of_events_to_process == 0)
    {
      number_of_events_to_process = the_input_ssd_energy_tree->GetEntries();
      output_message += "Will process: ";
      std::ostringstream oss;
      oss << number_of_events_to_process;
      output_message += oss.str();
      output_message += " events.\n";
    }
    else
    {
      if (number_of_events_to_process > the_input_ssd_energy_tree->GetEntries())
      {
        output_message += "Requested number of events to process higher than existing events in the root file. Will process ";
        number_of_events_to_process = the_input_ssd_energy_tree->GetEntries();
        std::ostringstream oss;
        oss << number_of_events_to_process;
        output_message += oss.str();
        output_message += " events.\n";
      }
    }
  }
  
  an_input_root_file = the_input_root_file;
  an_input_ssd_energy_tree = the_input_ssd_energy_tree;
  an_input_particle_truth_tree = the_input_particle_truth_tree;
  an_input_steps_truth_tree = the_input_steps_truth_tree;
  
  the_input_ssd_energy_tree->SetBranchAddress("plaque_id", &the_input_root_data_structure->plaque_number);
  the_input_ssd_energy_tree->SetBranchAddress("strip_id", &the_input_root_data_structure->strip_number);
  the_input_ssd_energy_tree->SetBranchAddress("ssd_total_energy_deposited", &the_input_root_data_structure->total_energy_values);
  the_input_ssd_energy_tree->SetBranchAddress("ssd_non_ionization_energy_deposited", &the_input_root_data_structure->non_ionization_energy_values);
  the_input_ssd_energy_tree->SetBranchAddress("contributing_tracks", &the_input_root_data_structure->contributing_tracks);
  
  the_input_particle_truth_tree->SetBranchAddress("track_id", &the_input_root_data_structure->track_id);
  the_input_particle_truth_tree->SetBranchAddress("particle_code", &the_input_root_data_structure->particle_code);
  the_input_particle_truth_tree->SetBranchAddress("parent_track_id", &the_input_root_data_structure->parent_track_id);
  the_input_particle_truth_tree->SetBranchAddress("creator_process", &the_input_root_data_structure->creator_process);
  the_input_particle_truth_tree->SetBranchAddress("x_vertex", &the_input_root_data_structure->x_vertex);
  the_input_particle_truth_tree->SetBranchAddress("y_vertex", &the_input_root_data_structure->y_vertex);
  the_input_particle_truth_tree->SetBranchAddress("z_vertex", &the_input_root_data_structure->z_vertex);
  the_input_particle_truth_tree->SetBranchAddress("px_vertex", &the_input_root_data_structure->px_vertex);
  the_input_particle_truth_tree->SetBranchAddress("py_vertex", &the_input_root_data_structure->py_vertex);
  the_input_particle_truth_tree->SetBranchAddress("pz_vertex", &the_input_root_data_structure->pz_vertex);
  the_input_particle_truth_tree->SetBranchAddress("ekin_vertex", &the_input_root_data_structure->ekin_vertex);
  the_input_particle_truth_tree->SetBranchAddress("visited_plaques", &the_input_root_data_structure->visited_plaques);

  the_input_steps_truth_tree->SetBranchAddress("step_track_id", &the_input_root_data_structure->step_track_id);
  the_input_steps_truth_tree->SetBranchAddress("step_volume_name", &the_input_root_data_structure->step_volume_name);
  the_input_steps_truth_tree->SetBranchAddress("step_plaque_id", &the_input_root_data_structure->step_plaque_id);
  the_input_steps_truth_tree->SetBranchAddress("step_x_i", &the_input_root_data_structure->step_x_i);
  the_input_steps_truth_tree->SetBranchAddress("step_y_i", &the_input_root_data_structure->step_y_i);
  the_input_steps_truth_tree->SetBranchAddress("step_z_i", &the_input_root_data_structure->step_z_i);
  the_input_steps_truth_tree->SetBranchAddress("step_px_i", &the_input_root_data_structure->step_px_i);
  the_input_steps_truth_tree->SetBranchAddress("step_py_i", &the_input_root_data_structure->step_py_i);
  the_input_steps_truth_tree->SetBranchAddress("step_pz_i", &the_input_root_data_structure->step_pz_i);
  the_input_steps_truth_tree->SetBranchAddress("step_x_f", &the_input_root_data_structure->step_x_f);
  the_input_steps_truth_tree->SetBranchAddress("step_y_f", &the_input_root_data_structure->step_y_f);
  the_input_steps_truth_tree->SetBranchAddress("step_z_f", &the_input_root_data_structure->step_z_f);
  the_input_steps_truth_tree->SetBranchAddress("step_px_f", &the_input_root_data_structure->step_px_f);
  the_input_steps_truth_tree->SetBranchAddress("step_py_f", &the_input_root_data_structure->step_py_f);
  the_input_steps_truth_tree->SetBranchAddress("step_pz_f", &the_input_root_data_structure->step_pz_f);
  the_input_steps_truth_tree->SetBranchAddress("step_initial_energy", &the_input_root_data_structure->step_initial_energy);
  the_input_steps_truth_tree->SetBranchAddress("step_final_energy", &the_input_root_data_structure->step_final_energy);
  the_input_steps_truth_tree->SetBranchAddress("step_process_name", &the_input_root_data_structure->step_process_name);
  
  if (DEBUG_PRINT_LEVEL > 1)
  {
    cout << output_message;
  }
  
  return initialization_succeeded;
}

int input_output_manager::get_number_of_events_to_process()
{
  return number_of_events_to_process;
}

bool input_output_manager::record_detailed_mc_data(input_root_data_structure* the_input_root_data_structure)
{
  // cout << "Recording MC detailed data." << endl;
  
  return true;
}

void input_output_manager::record_event()
{
  // this is a temporary debugging tool
  // should be commented out in the production version
  //
  
  int dummy_counter;
  
  out_file << "event" << endl;
  
  out_file << "event_class " << the_global_data_dispatcher->get_event_reco_output_data_structure()->event_class_by_cluster_multiplicity << endl;
  
  vector<pair<double, int>> order_plaque_8;
  vector<pair<double, int>> order_plaque_11_12;
  vector<pair<double, int>> order_plaque_19_20;

  order_plaque_8.clear();
  order_plaque_11_12.clear();
  order_plaque_19_20.clear();
  
  XYZVector dummy_point;
  
  for (int i = 0; i < the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.size(); i++)
  {
    if (the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.size() > 0)
    {
      dummy_counter = 0;
      
      for (int j = 0; j < the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.size(); j++)
      {
        // if (the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.at(j) > 3)
          //&& the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.at(j) < 9)
        {
          dummy_counter++;
        }
      }
      
      out_file << "particle " << the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).track_id << " "
      << the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).particle_code << " "
      << the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).ekin_vertex << " "
      << dummy_counter << " ";
      for (int j = 0; j < the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.size(); j++)
      {
        //if (the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.at(j) > 3)
          //&& the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.at(j) < 9)
        {
          out_file << the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.at(j) << " "
          << the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).position.at(j).X() << " "
          << the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).position.at(j).Y() << " "
          << the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).position.at(j).Z() << " ";
        }
      }
      out_file << endl;
      
      
      
      //if (the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.size() == 3)
      {
        //cout << "particle " << the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).particle_code << " "
                            //<< the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).ekin_vertex << " ";
        for (int j = 0; j < the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.size(); j++)
        {
          //if (the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.at(j) > 3)
          //&& the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.at(j) < 9)
          if (the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.at(j) > 8)
          {
            //cout << the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.at(j) << " "
            //<< the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).position.at(j).X() << " ";
            //<< the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).position.at(j).Y() << " ";
               //<< the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).position.at(j).Z() << " ";
          }
          
          if (the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.at(j) == 7)
          {
            order_plaque_8.push_back(make_pair(the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).position.at(j).Y(), i));
          }
          else if (the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.at(j) == 9 ||
                   the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.at(j) == 10)
          {
            order_plaque_11_12.push_back(make_pair(the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).position.at(j).Y(), i));
          }
          else if (the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.at(j) == 17 ||
                   the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).step_plate_index.at(j) == 18)
          {
            order_plaque_19_20.push_back(make_pair(the_global_data_dispatcher->get_event_mc_output_data_structure()->the_charged_tracks.at(i).position.at(j).Y(), i));
          }

        }
        //cout << endl;
      }
      
      
      
      
    }
  }
  
  sort(order_plaque_8.begin(), order_plaque_8.end());
  sort(order_plaque_11_12.begin(), order_plaque_11_12.end());
  sort(order_plaque_19_20.begin(), order_plaque_19_20.end());

  for (int i = 0; i < order_plaque_8.size(); i++)
  {
    //cout << order_plaque_8.at(i).second << " ";
  }
  //cout << endl;
  
  for (int i = 0; i < order_plaque_11_12.size(); i++)
  {
    //cout << order_plaque_11_12.at(i).second << " ";
  }
  //cout << endl;

  for (int i = 0; i < order_plaque_19_20.size(); i++)
  {
    //cout << order_plaque_19_20.at(i).second << " ";
  }
  //cout << endl;

  
  out_file << "number_of_clusters ";
  for (int i = 0; i < 16; i++)
  {
    out_file << the_global_data_dispatcher->get_event_reco_output_data_structure()->cluster_lines[i].size() << " ";
  }
  out_file << endl;
  
  out_file << "mc_vertex " << the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex.X() << " " << the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex.Y() << " " << the_global_data_dispatcher->get_event_mc_output_data_structure()->mc_vertex.Z() << endl;
  
  for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); i++)
  {
    out_file << "clusters " << i << " " << the_global_data_dispatcher->get_event_reco_output_data_structure()->cluster_lines[i].size() << " ";
    
    for (int j = 0; j < the_global_data_dispatcher->get_event_reco_output_data_structure()->cluster_lines[i].size(); j++)
    {
      out_file << the_global_data_dispatcher->get_event_reco_output_data_structure()->cluster_lines[i].at(j).point.X() << " "
               << the_global_data_dispatcher->get_event_reco_output_data_structure()->cluster_lines[i].at(j).point.Y() << " "
               << the_global_data_dispatcher->get_event_reco_output_data_structure()->cluster_lines[i].at(j).point.Z() << " "
               << the_global_data_dispatcher->get_event_reco_output_data_structure()->cluster_lines[i].at(j).direction.X() << " "
               << the_global_data_dispatcher->get_event_reco_output_data_structure()->cluster_lines[i].at(j).direction.Y() << " "
               << the_global_data_dispatcher->get_event_reco_output_data_structure()->cluster_lines[i].at(j).direction.Z() << " ";
    }
    out_file << endl;
  }
  
  out_file << "pattern " << the_global_data_dispatcher->get_event_reco_output_data_structure()->track_directions.size() << " " << the_global_data_dispatcher->get_event_reco_output_data_structure()->global_fit_quality << endl;
  out_file << "vertex " << the_global_data_dispatcher->get_event_reco_output_data_structure()->vertex.X() << " " << the_global_data_dispatcher->get_event_reco_output_data_structure()->vertex.Y() << " " << the_global_data_dispatcher->get_event_reco_output_data_structure()->vertex.Z() << endl;
  
  out_file << "tracks ";
  for (int i = 0; i < the_global_data_dispatcher->get_event_reco_output_data_structure()->track_directions.size(); i++)
  {
    out_file << i << " " << the_global_data_dispatcher->get_event_reco_output_data_structure()->track_directions.at(i).X() << " " << the_global_data_dispatcher->get_event_reco_output_data_structure()->track_directions.at(i).Y() << " " << the_global_data_dispatcher->get_event_reco_output_data_structure()->track_directions.at(i).Z() << " ";
  }
  out_file << endl;
  
  for (int i = 0; i < the_global_data_dispatcher->get_event_reco_output_data_structure()->reconstructed_tracks.size(); i++)
  {
    for (int j = 0; j < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); j++)
    {
      out_file << the_global_data_dispatcher->get_event_reco_output_data_structure()->reconstructed_tracks.at(i).track_plate_intersection_points[j].X() << " " << the_global_data_dispatcher->get_event_reco_output_data_structure()->reconstructed_tracks.at(i).track_plate_intersection_points[j].Y() << " "
      << the_global_data_dispatcher->get_event_reco_output_data_structure()->reconstructed_tracks.at(i).track_plate_intersection_points[j].Z() << " ";
      
    }
    out_file << endl;
  }
  
  out_file << "upstream_track_hits " << the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates() << endl;
  for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_upstream_plates(); i++)
  {
    dummy_point = the_global_data_dispatcher->get_event_reco_output_data_structure()->the_upstream_track.track_plate_intersection_points[i];
    out_file << i << " " << dummy_point.X() << " " << dummy_point.Y() << " " << dummy_point.Z() << " ";
  }
  out_file << endl;
  
  out_file << "endevent" << endl;

/*
cout << "^^^^^^^^" << endl;
cout << the_global_data_dispatcher->get_event_reco_output_data_structure()->reconstructed_tracks.size() << endl;
if (the_global_data_dispatcher->get_event_reco_output_data_structure()->reconstructed_tracks.size() == 1)
{
cout << the_global_data_dispatcher->get_event_reco_output_data_structure()->reconstructed_tracks[0].track_points.size() << endl;
}
*/
  
  return;
}

void input_output_manager::record_run()
{
  
  return;
}
