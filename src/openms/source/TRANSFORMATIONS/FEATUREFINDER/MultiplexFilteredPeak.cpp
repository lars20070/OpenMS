// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexSatellite.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexSatelliteProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredPeak.h>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{
  MultiplexFilteredPeak::MultiplexFilteredPeak(double mz, double rt, size_t mz_idx, size_t rt_idx) :
    mz_(mz), rt_(rt), mz_idx_(mz_idx), rt_idx_(rt_idx)
  {
  }

  double MultiplexFilteredPeak::getMZ() const
  {
    return mz_;
  }

  double MultiplexFilteredPeak::getRT() const
  {
    return rt_;
  }

  size_t MultiplexFilteredPeak::getMZidx() const
  {
    return mz_idx_;
  }

  size_t MultiplexFilteredPeak::getRTidx() const
  {
    return rt_idx_;
  }

  void MultiplexFilteredPeak::addSatellite(size_t rt_idx, size_t mz_idx, size_t pattern_idx)
  {
    satellites_.insert(std::make_pair(pattern_idx, MultiplexSatellite(rt_idx, mz_idx)));
  }
  
  void MultiplexFilteredPeak::addSatellite(MultiplexSatellite satellite, size_t pattern_idx)
  {
    satellites_.insert(std::make_pair(pattern_idx, satellite));
  }
  
  const std::multimap<size_t, MultiplexSatellite >& MultiplexFilteredPeak::getSatellites() const
  {
    return satellites_;
  }
  
  size_t MultiplexFilteredPeak::size() const
  {
    return satellites_.size();
  }
  
  size_t MultiplexFilteredPeak::sizeProfile() const
  {
    if (satellites_.size() == 0)
    {
      return 0;
    }
    
    // The size of the mz_ and intensity_ vectors of each satellite are identical.
    // Here we return the size of the first mz_ vector.
    return satellites_.begin()->second.getMZ().size();
  }
  
  void MultiplexFilteredPeak::updateCandidates(const MSExperiment& exp_picked, double mz_shift, std::vector<SplineSpectrum::Navigator>& navigators)
  {
    // construct the set of spline-interpolated satellites for this specific mz_profile
    for (std::multimap<size_t, MultiplexSatellite >::iterator satellite_it = satellites_.begin(); satellite_it != satellites_.end(); ++satellite_it)
    {
      // find indices of the peak
      size_t rt_idx = (satellite_it->second).getRTidx();
      size_t mz_idx = (satellite_it->second).getMZidx();
      
      // find peak itself
      MSExperiment::ConstIterator it_rt = exp_picked.begin();
      std::advance(it_rt, rt_idx);
      MSSpectrum<Peak1D>::ConstIterator it_mz = it_rt->begin();
      std::advance(it_mz, mz_idx);

      double rt_satellite = it_rt->getRT();
      double mz_satellite = it_mz->getMZ();
      
      // determine m/z and corresponding intensity
      double mz = mz_satellite + mz_shift;
      double intensity = navigators[rt_idx].eval(mz);

      (satellite_it->second).setMZTemp(mz);
      (satellite_it->second).setIntensityTemp(intensity);   
    }
  }
  
  void MultiplexFilteredPeak::pushPeakToResults(const MSExperiment& exp_picked)
  {
    // loop over satellites
    for (std::multimap<size_t, MultiplexSatellite >::iterator it = satellites_.begin(); it != satellites_.end(); ++it)
    {
      // find indices of the peak
      size_t rt_idx = (it->second).getRTidx();
      size_t mz_idx = (it->second).getMZidx();
        
      // find peak itself
      MSExperiment::ConstIterator it_rt = exp_picked.begin();
      std::advance(it_rt, rt_idx);
      MSSpectrum<Peak1D>::ConstIterator it_mz = it_rt->begin();
      std::advance(it_mz, mz_idx);
      
      // push to the result vectors
      (it->second).addMZ(it_mz->getMZ());
      (it->second).addIntensity(it_mz->getIntensity());
    }
  }
  
  void MultiplexFilteredPeak::pushDataPointToResults()
  {
    // loop over satellites
    for (std::multimap<size_t, MultiplexSatellite >::iterator it = satellites_.begin(); it != satellites_.end(); ++it)
    {
      // push candidate to the result vectors
      (it->second).addMZ((it->second).getMZTemp());
      (it->second).addIntensity((it->second).getIntensityTemp());
    }
  }
  
}
