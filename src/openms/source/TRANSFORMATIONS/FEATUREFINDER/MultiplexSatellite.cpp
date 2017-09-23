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
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexSatellite.h>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{
  MultiplexSatellite::MultiplexSatellite(size_t rt_idx, size_t mz_idx) :
    rt_idx_(rt_idx), mz_idx_(mz_idx)
  {
  }
  
  MultiplexSatellite::MultiplexSatellite(size_t rt_idx, size_t mz_idx, std::vector<double> mz, std::vector<double> intensity) :
    rt_idx_(rt_idx), mz_idx_(mz_idx), mz_(mz), intensity_(intensity)
  {
  }
  
  size_t MultiplexSatellite::getMZidx() const
  {
    return mz_idx_;
  }

  size_t MultiplexSatellite::getRTidx() const
  {
    return rt_idx_;
  }
  
  void MultiplexSatellite::setMZTemp(double mz_temp)
  {
    mz_temp_ = mz_temp;
  }
  
  void MultiplexSatellite::setIntensityTemp(double intensity_temp)
  {
    intensity_temp_ = intensity_temp;
  }
  
  double MultiplexSatellite::getMZTemp() const
  {
    return mz_temp_;
  }
  
  double MultiplexSatellite::getIntensityTemp() const
  {
    return intensity_temp_;
  }
  
  void MultiplexSatellite::addMZ(double mz)
  {
    mz_.push_back(mz);
  }
  
  void MultiplexSatellite::addIntensity(double intensity)
  {
    intensity_.push_back(intensity);
  }
  
  std::vector<double> MultiplexSatellite::getMZ() const
  {
    return mz_;
  }

  std::vector<double> MultiplexSatellite::getIntensity() const
  {
    return intensity_;
  }
  
}
