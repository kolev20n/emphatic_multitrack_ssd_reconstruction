// This is the cluster manager class. Its purpose is to combine
// neighbouring strips that are fired into a single cluster of strips.
// This is done separately for all plates.
//
// Blame: Nikolay Kolev, kolev20n@uregina.ca
//
// 2022
//

#ifndef cluster_manager_h
#define cluster_manager_h

#include <vector>
#include <string>

#include "common_data_structures_and_functions.h"

class hit_pair;

//class input_output_manager;
//class TTree;
//class TFile;
//class transformation_manager;
//class cluster_interpretation_manager;

class cluster_manager
{
public:
  cluster_manager();
  
  // hit_pair_list is an array (indexing each plate) of vectors of type hit_pair,
  // where hit_pair contains the strip number and the energy value (ADC value in
  // experimental data)
  // the clusters array (indexing plates again) is the output
  // the cluster interpretation method determines how the strips will be combined:
  // -> plain average value (not implemented, but easy to add)
  // -> average value weighted by the energy value (default)
  // -> other (not implemented)
  bool cluster(std::vector<hit_pair> hit_pair_list[], std::vector<double> clusters[], int actual_number_of_plaques, std::string cluster_interpretation_method);
private:

};

#endif
