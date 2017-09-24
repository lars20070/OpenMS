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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXFILTEREDPEAK_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXFILTEREDPEAK_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexSatellite.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexSatelliteProfile.h>

#include <map>
#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief data structure storing a single peak that passed all filters
   * 
   * Each filter result corresponds to a successful search for a particular
   * peak pattern in the centroided data. The actual m/z shifts seen in the filter
   * result might differ from the theoretical shifts listed in the peak pattern.
   * 
   * Each MultiplexFilteredPeak consists of a primary peak and a set of satellite peaks.
   * The primary peak is a peak in the mono-isotopic masstrace of the lightest peptide
   * in the multiplet. The satellite peaks are peaks that form the m/z shift pattern
   * relative to the primary peak within a retention time range rt_band_.
   * 
   * @see MultiplexPeakPattern
   */
  class OPENMS_DLLAPI MultiplexFilteredPeak
  {
    public:

    /**
     * @brief constructor
     */
    MultiplexFilteredPeak(double mz, double rt, size_t mz_idx, size_t rt_idx);

    /**
     * @brief returns m/z of the peak
     */
    double getMZ() const;
     
    /**
     * @brief returns RT of the peak
     */
    double getRT() const;
    
    /**
     * @brief returns the m/z index of the peak
     */
    size_t getMZidx() const;
     
    /**
     * @brief returns the RT index of the peak
     */
    size_t getRTidx() const;
    
    /**
     * @brief add a satellite peak
     */
    void addSatellite(size_t rt_idx, size_t mz_idx, size_t pattern_idx);
    
    void addSatellite(MultiplexSatellite satellite, size_t pattern_idx);
    
    /**
     * @brief return all satellite peaks
     */
    const std::multimap<size_t, MultiplexSatellite >& getSatellites() const;
    
    /**
     * @brief return number of satellite peaks
     */
    size_t size() const;
    
    /**
     * @brief return number of satellite data points
     */
    size_t sizeProfile() const;
    
    /**
     * @brief update the temporary data points in each satellite
     * 
     * These (m/z, intensity) pairs will be used for the subsequent averagine and peptide correlation filters.
     */
    void updateCandidate(const MSExperiment& exp_picked, double mz_shift, std::vector<SplineSpectrum::Navigator>& navigators);
    
    /**
     * @brief push peak to result vector (used in centroid mode)
     * 
     * push peak m/z and peak intensity of each satellite to the result vectors
     * (later on we will construct peptide features from them)
     */
    void pushPeakToResults(const MSExperiment& exp_picked);
    
    /**
     * @brief push candidate data point to result vector (used in profile mode)
     * 
     * push m/z and corresponding spline-interpolated intensity of a candidate data to the result vectors
     * (later on we will construct peptide features from them)
     */
    void pushDataPointToResults();
    
    private:
    /**
     * @brief position of the primary peak
     * 
     * Position of the primary peak in the m/z-RT plane in [Th, sec].
     * It is the input for the subsequent clustering step. 
     */
    double mz_;
    double rt_;
    
    /**
     * @brief indices of the primary peak position in the centroided experiment
     * 
     * Spectral index and peak index within the spectrum of the primary peak.
     * The indices are used to check the blacklist.
     */
    size_t mz_idx_;
    size_t rt_idx_;
    
    /**
     * @brief set of satellites
     * 
     * Mapping from a pattern index i.e. a specific mass trace to all peaks forming
     * the pattern. The primary peak is part of the satellite peak set.
     * 
     * pattern_idx -> (rt_idx, mz_idx)
     * 
     * Note that we store only indices, not iterators or pointers. We filter
     * 'white' experiments, but all indices refer to the original experiment.
     * White experiments are temporary (for each pattern), but the original
     * <exp_picked_> experiment is permanent.
     */
    std::multimap<size_t, MultiplexSatellite > satellites_;
    
  };
  
}

#endif /* MULTIPLEXFILTEREDPEAK_H */
