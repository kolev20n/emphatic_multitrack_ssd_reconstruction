// This class is the engine of the silicon strip reconstruction. It
// contains the main loop, controls the algorithm, performs the 
// reconstruction steps (depending on the case) and controls the 
// output (general root or specific other output). It mostly contains
// pointers to data structures, to the geometry and different managers,
// which are the classes that actually perform the various reconstruction
// tasks; most notably it drives the track finding manager, which is the 
// class that performs the majority of reconstruction tasks.
//
// Blame: Nikolay Kolev, kolev20n@uregina.ca
//
// 2022
//

#ifndef reconstruction_manager_h
#define reconstruction_manager_h

#include "common_data_structures_and_functions.h"
#include "partial_geometry.h"

class TTree;
class TFile;

class input_output_manager;
class partial_geometry;
class cluster_manager;
class transformation_manager;
class track_finding_manager;

class reconstruction_manager
{
public:
  reconstruction_manager(input_output_manager* an_input_output_manager);
  void process_events();

private:
  input_output_manager* the_input_output_manager;
  cluster_manager* the_cluster_manager;
  transformation_manager* the_transformation_manager;
  track_finding_manager* the_track_finding_manager;

  partial_geometry* the_geometry;
  
  TFile* the_input_root_file;
  TTree* the_input_ssd_energy_tree;
  TTree* the_input_particle_truth_tree;
  TTree* the_input_steps_truth_tree;
  input_root_data_structure* the_input_root_data_structure;

  run_mc_output_data_structure* the_run_mc_output_data_structure;
  run_reco_output_data_structure* the_run_reco_output_data_structure;
  event_mc_output_data_structure* the_event_mc_output_data_structure;
  event_reco_output_data_structure* the_event_reco_output_data_structure;
  
  void calculate_particle_data();
  
  //track_finding_manager* the_track_finding_manager;
};

#endif
