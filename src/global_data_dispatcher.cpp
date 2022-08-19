#include "global_data_dispatcher.h"
#include "common_data_structures_and_functions.h"
#include "partial_geometry.h"

using namespace std;
using namespace ROOT::Math;

partial_geometry* global_data_dispatcher::get_partial_geometry()
{
  return the_geometry;
}

void global_data_dispatcher::register_geometry(partial_geometry* a_geometry)
{
  the_geometry = a_geometry;
}

MagneticField* global_data_dispatcher::get_magnetic_field()
{
  return the_magnetic_field;
}

void global_data_dispatcher::register_magnetic_field(MagneticField* a_magnetic_field)
{
  the_magnetic_field = a_magnetic_field;
}

input_root_data_structure* global_data_dispatcher::get_input_root_data_structure()
{
  return the_input_root_data_structure;
}

void global_data_dispatcher::register_input_root_data_structure(input_root_data_structure* an_input_root_data_structure)
{
  the_input_root_data_structure = an_input_root_data_structure;
}

reconstruction_options* global_data_dispatcher::get_reconstruction_options()
{
  return the_reconstruction_options;
}

void global_data_dispatcher::register_reconstruction_options(reconstruction_options* a_reconstruction_options)
{
  the_reconstruction_options = a_reconstruction_options;
}

reconstruction_algorithm* global_data_dispatcher::get_reconstruction_algorithm()
{
  return the_reconstruction_algorithm;
}

void global_data_dispatcher::register_reconstruction_algorithm(reconstruction_algorithm* a_reconstruction_algorithm)
{
  the_reconstruction_algorithm = a_reconstruction_algorithm;
}

run_mc_output_data_structure* global_data_dispatcher::get_run_mc_output_data_structure()
{
  return the_run_mc_output_data_structure;
}

void global_data_dispatcher::register_run_mc_output_data_structure(run_mc_output_data_structure* a_run_mc_output_data_structure)
{
  the_run_mc_output_data_structure = a_run_mc_output_data_structure;
}

run_reco_output_data_structure* global_data_dispatcher::get_run_reco_output_data_structure()
{
  return the_run_reco_output_data_structure;
}

void global_data_dispatcher::register_run_reco_output_data_structure(run_reco_output_data_structure* a_run_reco_output_data_structure)
{
  the_run_reco_output_data_structure = a_run_reco_output_data_structure;
}

event_mc_output_data_structure* global_data_dispatcher::get_event_mc_output_data_structure()
{
  return the_event_mc_output_data_structure;
}

void global_data_dispatcher::register_event_mc_output_data_structure(event_mc_output_data_structure* an_event_mc_output_data_structure)
{
  the_event_mc_output_data_structure = an_event_mc_output_data_structure;
}

event_reco_output_data_structure* global_data_dispatcher::get_event_reco_output_data_structure()
{
  return the_event_reco_output_data_structure;
}

void global_data_dispatcher::register_event_reco_output_data_structure(event_reco_output_data_structure* an_event_reco_output_data_structure)
{
  the_event_reco_output_data_structure = an_event_reco_output_data_structure;
}




