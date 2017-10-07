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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXFILTERING_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXFILTERING_H

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResult.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredPeak.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief base class for filtering centroided and profile data for peak patterns
   *
   * The algorithm searches for patterns of multiple peptides in the data.
   * The peptides appear as characteristic patterns of isotopic peaks in
   * MS1 spectra. We first search the centroided data, and optionally in
   * a second step the spline interpolated profile data. For each
   * peak pattern the algorithm generates a filter result.
   *
   * The algorithm differs slightly for centroided and profile input data.
   * This base class comprises code common to both. The two child classes
   * MultiplexFilteringCentroided and MultiplexFilteringProfile contain
   * specific functions and the primary filter() method.
   *
   * @see MultiplexIsotopicPeakPattern
   * @see MultiplexFilterResult
   * @see MultiplexFilteringCentroided
   * @see MultiplexFilteringProfile
   */
  class OPENMS_DLLAPI MultiplexFiltering :
    public ProgressLogger
  {
public:
    /**
     * @brief type for peak blacklisting
     * 
     * white    white in this and subsequent patterns
     * grey     white in this pattern and black in subsequent patterns
     * black    black in this and in subsequent patterns
     * 
     * We assume that one peak cannot belong to two or more patterns
     * i.e. peptides at the same time.
     */
    enum BlacklistEntry
    {
      white,
      grey,
      black
    };
    
    /**
     * @brief index mapping from a 'white' experiment to its original experiment
     * 
     * An MSExperiment contains a set of spectra each containing a number of peaks.
     * In the course of the filtering, some peaks are blacklisted since they are
     * identified to belong to a certain pattern i.e. peptide. An experiment in
     * which blacklisted peaks are removed is called 'white'. White spectra
     * contain fewer peaks than their corresponding primary spectra. Consequently,
     * their indices are shifted. The type maps a peak index in a 'white'
     * spectrum back to its original spectrum.
     */
    typedef std::vector<std::map<int, int> > White2Original;

    /**
     * @brief constructor
     *
     * @param exp_picked    experimental data in centroid mode
     * @param patterns    patterns of isotopic peaks to be searched for
     * @param isotopes_per_peptide_min    minimum number of isotopic peaks in peptides
     * @param isotopes_per_peptide_max    maximum number of isotopic peaks in peptides
     * @param intensity_cutoff    intensity cutoff
     * @param rt_band    RT range used for filtering
     * @param mz_tolerance    error margin in m/z for matching expected patterns to experimental data
     * @param mz_tolerance_unit    unit for mz_tolerance, ppm (true), Da (false)
     * @param peptide_similarity    similarity score for two peptides in the same multiplet
     * @param averagine_similarity    similarity score for peptide isotope pattern and averagine model
     * @param averagine_similarity_scaling    scaling factor x for the averagine similarity parameter p when detecting peptide singlets. With p' = p + x(1-p). 
     */
    MultiplexFiltering(const MSExperiment& exp_picked, const std::vector<MultiplexIsotopicPeakPattern> patterns, int isotopes_per_peptide_min, int isotopes_per_peptide_max, double intensity_cutoff, double rt_band, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, String averagine_type="peptide");

protected:
    /**
     * @brief construct an MS experiment from exp_picked_ containing
     * peaks which have not previously blacklisted in blacklist_
     * 
     * @param mapping    index mapping of 'white' peak positions to their position in the corresponding, original spectrum 
     */
    MSExperiment getWhiteMSExperiment_(White2Original& mapping); 

    /**
     * @brief check for significant peak
     *
     * @param mz    position where the peak is expected
     * @param mz_tolerance    m/z tolerance within the peak may lie
     * @param it_rt    pointer to the spectrum
     * @param intensity_first_peak    intensity to compare to
     *
     * @return boolean if there is a significant peak
     */
    bool checkForSignificantPeak_(double mz, double mz_tolerance, MSExperiment::ConstIterator& it_rt, double intensity_first_peak) const;

    /**
     * @brief check if there are enough peaks in the RT band to form the pattern
     *
     * Checks if there are peaks at m/z positions corresponding to the pattern
     * and that the primary peak position is not blacklisted.
     *
     * @param it_mz    m/z iterator of the primary
     * @param it_rt_begin    RT iterator of the very first spectrum of the experiment (needed to determine indices)
     * @param it_rt_band_begin    RT iterator of the first spectrum in the RT band
     * @param it_rt_band_end    RT iterator of the last spectrum in the RT band
     * @param pattern    m/z pattern to search for
     * @param peak    filter result output
     *
     * @return boolean if this filter was passed i.e. there are <isotopes_per_peptide_min_> or more mass traces which form the pattern.
     */
    bool filterPeakPositions_(const MSSpectrum<Peak1D>::ConstIterator& it_mz, White2Original& index_mapping, const MSExperiment::ConstIterator& it_rt_begin, const MSExperiment::ConstIterator& it_rt_band_begin, const MSExperiment::ConstIterator& it_rt_band_end, const MultiplexIsotopicPeakPattern& pattern, MultiplexFilteredPeak& peak) const;

    /**
     * @brief blacklist this peak
     * 
     * Blacklist all satellites associated with this peak.
     * 
     * @param peak    peak to be blacklisted
     */
    void blacklistPeak_(const MultiplexFilteredPeak& peak);
    
    /**
     * @brief blacklist this peak
     * 
     * Each of the satellites is associated with a specific mass trace. We blacklist
     * all peaks in these mass traces (even if they are not a satellite) extending them
     * by a margin <rt_band_>.  
     * 
     * @param peak    peak to be blacklisted
     * @param pattern_idx    index of the pattern in <patterns_>
     */
    void blacklistPeak_(const MultiplexFilteredPeak& peak, unsigned pattern_idx);

    /**
     * @brief turn grey blacklist_ entries into black ones
     * 
     * Grey entries function as white in the current pattern but black in subsequent patterns,
     * i.e. at the end of a pattern these entries need to be turned black.
     */
    void ungreyBlacklist_();

    /**
    * @brief centroided experimental data
    */
    MSExperiment exp_picked_;

    /**
    * @brief auxiliary structs for blacklisting
    */
    std::vector<std::vector<BlacklistEntry> > blacklist_;

    /**
     * @brief list of peak patterns
     */
    std::vector<MultiplexIsotopicPeakPattern> patterns_;

    /**
     * @brief minimum number of isotopic peaks per peptide
     */
    size_t isotopes_per_peptide_min_;

    /**
     * @brief maximum number of isotopic peaks per peptide
     */
    size_t isotopes_per_peptide_max_;

    /**
     * @brief intensity cutoff
     */
    double intensity_cutoff_;

    /**
     * @brief RT range used for filtering
     */
    double rt_band_;
    
    /**
     * @brief m/z shift tolerance
     */
    double mz_tolerance_;

    /**
     * @brief unit for m/z shift tolerance (ppm - true, Da - false)
     */
    bool mz_tolerance_unit_;

    /**
      * @brief peptide similarity
      */
    double peptide_similarity_;

    /**
     * @brief averagine similarity
     */
    double averagine_similarity_;

    /**
     * @brief averagine similarity scaling
     */
    double averagine_similarity_scaling_;

    /**
     * @brief type of averagine to use
     */
    String averagine_type_;

  };

}

#endif /* MULTIPLEXFILTERING_H_ */
