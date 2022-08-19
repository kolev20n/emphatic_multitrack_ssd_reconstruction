// This class follows Munuit 2 instructions to fit several tracks
// as straight line through a common vertex and a set of skew lines,
// each track with its own set assigned to it.
//
// Blame: Nikolay Kolev, kolev20n@uregina.ca
//
// 2022
//

#ifndef my_global_fit_no_field_fcn_h
#define my_global_fit_no_field_fcn_h

#include "Minuit2/FCNBase.h"
#include <TMath.h>
#include <Math/Vector3D.h>

#include <vector>

#include "common_data_structures_and_functions.h"
#include "partial_geometry.h"

namespace ROOT
{
  namespace Minuit2
  {
    class my_global_fit_no_field_fcn : public FCNBase
    {
    public:
      my_global_fit_no_field_fcn(const pattern& a_pattern, const track& an_upstream_track)
      : the_pattern(a_pattern), the_upstream_track(an_upstream_track), fErrorDef(1.) {}
      
      my_global_fit_no_field_fcn() {}
      
      virtual double Up() const {return fErrorDef;}
      virtual double operator()(const std::vector<double>&) const;
      
      void SetErrorDef(double def) {fErrorDef = def;}
      
    private:
      pattern the_pattern;
      track the_upstream_track;
      
      double fErrorDef;
    };
    
  }
  
}

#endif
