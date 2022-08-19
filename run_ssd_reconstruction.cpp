// This is the main function of the track reconstruction code
// for the silicon strip detectors of the EMPHATIC experiment.
//
// Blame: Nikolay Kolev, kolev20n@uregina.ca
//
// 2022
//

#include <iostream>

#include "input_output_manager.h"
#include "reconstruction_manager.h"
#include "global_data_dispatcher.h"
#include "partial_geometry.h"

using namespace std;

// The global data dispatcher is a global solution for data exchange.
// Every data structure used for data exchange between classes and/or functions
// can be registered with it and then it will provide a pointer to it from
// anywhere in the code.
// This should be removed in the production version, where the data flow will be
// finalized.
//
global_data_dispatcher* the_global_data_dispatcher;

int main(int argc, char **argv)
{
  if (DEBUG_PRINT_LEVEL > 0)
  {
    cout << "EMPHATIC ssd reconstruction starting." << endl;
  }
  
  the_global_data_dispatcher = new global_data_dispatcher;

  input_output_manager* the_input_output_manager = new input_output_manager();

  if (!the_input_output_manager->initialize_files(argc, argv))
  {
    cout << "Error 1001: The initialize_files operation of the_input_output_manager failed. Quitting..." << endl;
    return 1001;
  }
  
  partial_geometry* the_geometry = new partial_geometry();
  
  the_global_data_dispatcher->register_geometry(the_geometry);
  
  if (!the_input_output_manager->initialize_geometry(the_geometry))
  {
    cout << "Error 1002: The initialize_geometry operation of the_input_output_manager failed. Quitting..." << endl;
    return 1002;
  }
  
  reconstruction_manager the_reconstruction_manager(the_input_output_manager);
  
  // The processing part starts here. The process_events() function contains the main
  // event loop.
  //
  if (DEBUG_PRINT_LEVEL > 0)
  {
    cout << "Starting event loop." << endl;
  }
  the_reconstruction_manager.process_events();
  
  // Since output is written during the processing of each event, there is not much to do here.
  // Maybe the overall statistics can be printed (from the run_reco_output_data_structure via input_output_manager)
  // or something of that sort.
  //
  
  if (DEBUG_PRINT_LEVEL > 0)
  {
    cout << "EMPHATIC ssd reconstruction ending." << endl;
  }
  
  return 0;
}
