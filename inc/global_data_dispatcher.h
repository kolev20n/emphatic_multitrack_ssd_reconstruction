// This is the global_data_dispatcher class. It serves as a
// dispatcher for data-exchange classes, as well as some commonly
// used classes such as partial_geometry. The main reason for its
// existence is that during development the chain of function calls
// and which class is using what data (normally with several layers
// of calles separating them) changes quite often, so it is much more
// convenient to provide data access to any class or function needing
// it using the global dispather.
// This may change in the production version, because global data access
// os not good OOP.
// All items follow the same pattern: once a data structure is created,
// it must be registered with its own register method. A getter function
// must also be provided that returns a pointer to the data item.
//
// Blame: Nikolay Kolev, kolev20n@uregina.ca
//
// 2022
//

#ifndef global_data_dispatcher_h
#define global_data_dispatcher_h

class partial_geometry;
class input_root_data_structure;
class reconstruction_options;
class reconstruction_algorithm;
class run_mc_output_data_structure;
class run_reco_output_data_structure;
class event_mc_output_data_structure;
class event_reco_output_data_structure;
class MagneticField;

class global_data_dispatcher
{
public:
  partial_geometry* get_partial_geometry();
  void register_geometry(partial_geometry* a_geometry);
  
  MagneticField* get_magnetic_field();
  void register_magnetic_field(MagneticField* a_magnetic_field);
  
  input_root_data_structure* get_input_root_data_structure();
  void register_input_root_data_structure(input_root_data_structure* an_input_root_data_structure);
  
  reconstruction_options* get_reconstruction_options();
  void register_reconstruction_options(reconstruction_options* a_reconstruction_options);
  
  reconstruction_algorithm* get_reconstruction_algorithm();
  void register_reconstruction_algorithm(reconstruction_algorithm* a_reconstruction_algorithm);

  run_mc_output_data_structure* get_run_mc_output_data_structure();
  void register_run_mc_output_data_structure(run_mc_output_data_structure* a_run_mc_output_data_structure);
  
  run_reco_output_data_structure* get_run_reco_output_data_structure();
  void register_run_reco_output_data_structure(run_reco_output_data_structure* a_run_reco_output_data_structure);
  
  event_mc_output_data_structure* get_event_mc_output_data_structure();
  void register_event_mc_output_data_structure(event_mc_output_data_structure* an_event_mc_output_data_structure);
  
  event_reco_output_data_structure* get_event_reco_output_data_structure();
  void register_event_reco_output_data_structure(event_reco_output_data_structure* an_event_reco_output_data_structure);
  
private:
  partial_geometry* the_geometry;
  input_root_data_structure* the_input_root_data_structure;
  reconstruction_options* the_reconstruction_options; 
  reconstruction_algorithm* the_reconstruction_algorithm;
  run_mc_output_data_structure* the_run_mc_output_data_structure;
  run_reco_output_data_structure* the_run_reco_output_data_structure;
  event_mc_output_data_structure* the_event_mc_output_data_structure;
  event_reco_output_data_structure* the_event_reco_output_data_structure;
  MagneticField* the_magnetic_field;
};

#endif
