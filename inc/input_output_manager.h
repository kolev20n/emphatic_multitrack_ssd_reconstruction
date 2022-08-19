// This class manages the input and output communications.
//
// Blame: Nikolay Kolev, kolev20n@uregina.ca
//
// 2022
//

#ifndef input_output_manager_h
#define input_output_manager_h

#include <string>
#include <fstream>

#include "common_data_structures_and_functions.h"
#include "matej_minimum_data_structures_and_functions.h"

class TFile;
class TTree;
class partial_geometry;

//class histogram_manager;

class input_output_manager
{
public:
  input_output_manager();
  ~input_output_manager();
  
  // the expected input is the single configuration file, which
  // contains the information for other files (geometry, input root,
  // output root, etc.), cuts and other configuration information
  bool initialize_files(int argc, char** argv);
  
  // reads from the geometry json file and initializes the geometry class
  bool initialize_geometry(partial_geometry* the_geometry);
  
  // prepares the root input trees for reading later on by reconstruction_manager;
  // all three trees are connected to the input_root_data_structure, which has the
  // data variables for all of them
  bool initialize_root_trees(TFile*& an_input_root_file, TTree*& an_input_ssd_energy_tree, TTree*& an_input_particle_truth_tree, TTree*& an_input_steps_truth_tree, input_root_data_structure*& the_input_root_data_structure);
  
  int get_number_of_events_to_process();
  
  // int get_number_of_plates();
  
  // provides access to the reconstruction options (i.e. cuts and tunings for the
  // reconstruction algorithm)
  reconstruction_options get_reconstruction_options();

  // dump mc data from input
  bool record_detailed_mc_data(input_root_data_structure* the_input_root_data_structure);
  
  void record_event();
  void record_run();
  
private:
  std::string output_message;
  
  std::string magnetic_field_file_name;
  std::string geometry_file_name;
  std::string input_root_file_name;
  TFile* the_input_root_file;
  TTree* the_input_ssd_energy_tree;
  TTree* the_input_particle_truth_tree;
  TTree* the_input_steps_truth_tree;
  
  std::string output_root_file_name;
  TFile* the_output_root_file;
  TTree* the_run_output_tree;
  TTree* the_event_mc_output_tree;
  TTree* the_event_reco_output_tree;

  std::string input_root_ssd_energy_tree_name;
  std::string input_root_particle_truth_tree_name;
  std::string input_root_steps_truth_tree_name;

  std::string run_output_tree_name = "run_output";
  std::string event_mc_output_tree_name = "event_mc_output";
  std::string event_reco_output_tree_name = "event_reco_output";

  int number_of_events_to_process;
  
  reconstruction_options* the_reconstruction_options;
  reconstruction_algorithm* the_algorithm;
  MagneticField* the_magnetic_field;
  
  std::ofstream out_file;
};

#endif
