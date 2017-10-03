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
// $Authors: Lars Nilse, Timo Sachsenberg, Samuel Wein, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/FORMAT/MzQuantMLFile.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>

#include <OpenMS/METADATA/MSQuantifications.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzQuantMLFile.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexClustering.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexSatellite.h>
#include <OpenMS/COMPARISON/CLUSTERING/GridBasedCluster.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

//Contrib includes
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <QDir>

//std includes
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <limits>
#include <locale>
#include <iomanip>

using namespace std;
using namespace OpenMS;
using namespace boost::math;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_FeatureFinderMultiplex FeatureFinderMultiplex

  @brief Detects peptide pairs in LC-MS data and determines their relative abundance.

<CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FeatureFinderMultiplex \f$ \longrightarrow \f$</td>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_IDMapper</td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileFilter </td>
    </tr>
  </table>
</CENTER>

  FeatureFinderMultiplex is a tool for the fully automated analysis of quantitative proteomics data. It detects pairs of isotopic envelopes with fixed m/z separation. It requires no prior sequence identification of the peptides. In what follows we outline the algorithm.

  <b>Algorithm</b>

  The algorithm is divided into three parts: filtering, clustering and linear fitting, see Fig. (d), (e) and (f). In the following discussion let us consider a particular mass spectrum at retention time 1350 s, see Fig. (a). It contains a peptide of mass 1492 Da and its 6 Da heavier labelled counterpart. Both are doubly charged in this instance. Their isotopic envelopes therefore appear at 746 and 749 in the spectrum. The isotopic peaks within each envelope are separated by 0.5. The spectrum was recorded at finite intervals. In order to read accurate intensities at arbitrary m/z we spline-fit over the data, see Fig. (b).

  We would like to search for such peptide pairs in our LC-MS data set. As a warm-up let us consider a standard intensity cut-off filter, see Fig. (c). Scanning through the entire m/z range (red dot) only data points with intensities above a certain threshold pass the filter. Unlike such a local filter, the filter used in our algorithm takes intensities at a range of m/z positions into account, see Fig. (d). A data point (red dot) passes if
  - all six intensities at m/z, m/z+0.5, m/z+1, m/z+3, m/z+3.5 and m/z+4 lie above a certain threshold,
  - the intensity profiles in neighbourhoods around all six m/z positions show a good correlation and
  - the relative intensity ratios within a peptide agree up to a factor with the ratios of a theoretic averagine model.

  Let us now filter not only a single spectrum but all spectra in our data set. Data points that pass the filter form clusters in the t-m/z plane, see Fig. (e). Each cluster corresponds to the mono-isotopic mass trace of the lightest peptide of a SILAC pattern. We now use hierarchical clustering methods to assign each data point to a specific cluster. The optimum number of clusters is determined by maximizing the silhouette width of the partitioning. Each data point in a cluster corresponds to three pairs of intensities (at [m/z, m/z+3], [m/z+0.5, m/z+3.5] and [m/z+1, m/z+4]). A plot of all intensity pairs in a cluster shows a clear linear correlation, see Fig. (f). Using linear regression we can determine the relative amounts of labelled and unlabelled peptides in the sample.

  @image html SILACAnalyzer_algorithm.png

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_FeatureFinderMultiplex.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_FeatureFinderMultiplex.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderMultiplex :
  public TOPPBase
{
private:

  // input and output files
  String in_;
  String out_;
  String out_features_;
  String out_mzq_;

  // section "algorithm"
  String labels_;
  std::vector<std::vector<String> > samples_labels_;
  unsigned charge_min_;
  unsigned charge_max_;
  unsigned missed_cleavages_;
  unsigned isotopes_per_peptide_min_;
  unsigned isotopes_per_peptide_max_;
  double rt_typical_;
  double rt_band_;
  double rt_min_;
  double mz_tolerance_;
  bool mz_unit_; // ppm (true), Da (false)
  double intensity_cutoff_;
  double peptide_similarity_;
  double averagine_similarity_;
  double averagine_similarity_scaling_;
  bool knock_out_;
  String spectrum_type_;
  String averagine_type_;

  // section "labels"
  map<String, double> label_mass_shift_;
  
  // experimental data
  MSExperiment exp_profile_;
  MSExperiment exp_centroid_;
  std::vector<SplineSpectrum> exp_spline_profile_;
  std::vector<SplineSpectrum::Navigator> navigators_;
  
public:
  TOPPFeatureFinderMultiplex() :
    TOPPBase("FeatureFinderMultiplex", "Determination of peak ratios in LC-MS data", true),
    charge_min_(1), charge_max_(1), missed_cleavages_(0), isotopes_per_peptide_min_(1), isotopes_per_peptide_max_(1), rt_typical_(0.0), rt_band_(0.0), rt_min_(0.0),
    mz_tolerance_(0.0), mz_unit_(true), intensity_cutoff_(0.0), peptide_similarity_(0.0), averagine_similarity_(0.0), averagine_similarity_scaling_(0.0), knock_out_(false)
  {
  }

  typedef std::vector<double> MassPattern; // list of mass shifts

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "LC-MS dataset in centroid or profile mode");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Set of all identified peptide groups (i.e. peptide pairs or triplets or singlets or ..). The m/z-RT positions correspond to the lightest peptide in each group.", false);
    setValidFormats_("out", ListUtils::create<String>("consensusXML"));
    registerOutputFile_("out_features", "<file>", "", "Optional output file containing the individual peptide features in \'out\'.", false, true);
    setValidFormats_("out_features", ListUtils::create<String>("featureXML"));
    registerOutputFile_("out_mzq", "<file>", "", "Optional output file of MzQuantML.", false, true);
    setValidFormats_("out_mzq", ListUtils::create<String>("mzq"));

    registerSubsection_("algorithm", "Parameters for the algorithm.");
    registerSubsection_("labels", "Isotopic labels that can be specified in section \'algorithm:labels\'.");
  }

  // create parameters for sections (set default values and restrictions)
  Param getSubsectionDefaults_(const String& section) const
  {
    Param defaults;

    if (section == "algorithm")
    {
      defaults.setValue("labels", "[][Lys8,Arg10]", "Labels used for labelling the samples. [...] specifies the labels for a single sample. For example\n\n[][Lys8,Arg10]        ... SILAC\n[][Lys4,Arg6][Lys8,Arg10]        ... triple-SILAC\n[Dimethyl0][Dimethyl6]        ... Dimethyl\n[Dimethyl0][Dimethyl4][Dimethyl8]        ... triple Dimethyl\n[ICPL0][ICPL4][ICPL6][ICPL10]        ... ICPL");
      defaults.setValue("charge", "1:4", "Range of charge states in the sample, i.e. min charge : max charge.");
      defaults.setValue("isotopes_per_peptide", "3:6", "Range of isotopes per peptide in the sample. For example 3:6, if isotopic peptide patterns in the sample consist of either three, four, five or six isotopic peaks. ", ListUtils::create<String>("advanced"));
      defaults.setValue("rt_typical", 40.0, "Typical retention time [s] over which a characteristic peptide elutes. (This is not an upper bound. Peptides that elute for longer will be reported.)");
      defaults.setMinFloat("rt_typical", 0.0);
      defaults.setValue("rt_band", 10.0, "RT band which is taken into considerations when filtering.TODO docu");
      defaults.setMinFloat("rt_band", 0.0);
      defaults.setValue("rt_min", 2.0, "Lower bound for the retention time [s]. (Any peptides seen for a shorter time period are not reported.)");
      defaults.setMinFloat("rt_min", 0.0);
      defaults.setValue("mz_tolerance", 6.0, "m/z tolerance for search of peak patterns.");
      defaults.setMinFloat("mz_tolerance", 0.0);
      defaults.setValue("mz_unit", "ppm", "Unit of the 'mz_tolerance' parameter.");
      defaults.setValidStrings("mz_unit", ListUtils::create<String>("Da,ppm"));
      defaults.setValue("intensity_cutoff", 1000.0, "Lower bound for the intensity of isotopic peaks.");
      defaults.setMinFloat("intensity_cutoff", 0.0);
      defaults.setValue("peptide_similarity", 0.5, "Two peptides in a multiplet are expected to have the same isotopic pattern. This parameter is a lower bound on their similarity.");
      defaults.setMinFloat("peptide_similarity", -1.0);
      defaults.setMaxFloat("peptide_similarity", 1.0);
      defaults.setValue("averagine_similarity", 0.4, "The isotopic pattern of a peptide should resemble the averagine model at this m/z position. This parameter is a lower bound on similarity between measured isotopic pattern and the averagine model.");
      defaults.setMinFloat("averagine_similarity", -1.0);
      defaults.setMaxFloat("averagine_similarity", 1.0);
      defaults.setValue("averagine_similarity_scaling", 0.75, "Let x denote this scaling factor, and p the averagine similarity parameter. For the detection of single peptides, the averagine parameter p is replaced by p' = p + x(1-p), i.e. x = 0 -> p' = p and x = 1 -> p' = 1. (For knock_out = true, peptide doublets and singlets are detected simulataneously. For singlets, the peptide similarity filter is irreleavant. In order to compensate for this 'missing filter', the averagine parameter p is replaced by the more restrictive p' when searching for singlets.)", ListUtils::create<String>("advanced"));
      defaults.setMinFloat("averagine_similarity_scaling", 0.0);
      defaults.setMaxFloat("averagine_similarity_scaling", 1.0);
      defaults.setValue("missed_cleavages", 0, "Maximum number of missed cleavages due to incomplete digestion. (Only relevant if enzymatic cutting site coincides with labelling site. For example, Arg/Lys in the case of trypsin digestion and SILAC labelling.)");
      defaults.setMinInt("missed_cleavages", 0);
      defaults.setValue("knock_out", "false", "Is it likely that knock-outs are present? (Supported for doublex, triplex and quadruplex experiments only.)", ListUtils::create<String>("advanced"));
      defaults.setValidStrings("knock_out", ListUtils::create<String>("true,false"));
      defaults.setValue("spectrum_type", "automatic", "Type of MS1 spectra in input mzML file. 'automatic' determines the spectrum type directly from the input mzML file.", ListUtils::create<String>("advanced"));
      defaults.setValidStrings("spectrum_type", ListUtils::create<String>("profile,centroid,automatic"));
      defaults.setValue("averagine_type","peptide","The type of averagine to use, currently RNA, DNA or peptide", ListUtils::create<String>("advanced"));
      defaults.setValidStrings("averagine_type", ListUtils::create<String>("peptide,RNA,DNA"));
    }

    if (section == "labels")
    {
      MultiplexDeltaMassesGenerator generator;
      Param p = generator.getParameters();
      
      for (Param::ParamIterator it = p.begin(); it != p.end(); ++it)
      {
        defaults.setValue(it->name, it->value, it->description, ListUtils::create<String>("advanced"));
        defaults.setMinFloat(it->name, 0.0);
      }
    }

    return defaults;
  }

  /**
   * @brief process parameters of 'input/output' section
   */
  void getParameters_in_out_()
  {
    in_ = getStringOption_("in");
    out_ = getStringOption_("out");
    out_features_ = getStringOption_("out_features");
    out_mzq_ = getStringOption_("out_mzq");
  }

  /**
   * @brief process parameters of 'algorithm' section
   */
  void getParameters_algorithm_()
  {
    // get selected labels
    labels_ = getParam_().getValue("algorithm:labels");
    samples_labels_ = splitLabelString_();

    // get selected charge range
    String charge_string = getParam_().getValue("algorithm:charge");
    double charge_min_temp, charge_max_temp;
    parseRange_(charge_string, charge_min_temp, charge_max_temp);
    charge_min_ = charge_min_temp;
    charge_max_ = charge_max_temp;
    if (charge_min_ > charge_max_)
    {
      swap(charge_min_, charge_max_);
    }

    // get isotopes per peptide range
    String isotopes_per_peptide_string = getParam_().getValue("algorithm:isotopes_per_peptide");
    double isotopes_per_peptide_min_temp, isotopes_per_peptide_max_temp;
    parseRange_(isotopes_per_peptide_string, isotopes_per_peptide_min_temp, isotopes_per_peptide_max_temp);
    isotopes_per_peptide_min_ = isotopes_per_peptide_min_temp;
    isotopes_per_peptide_max_ = isotopes_per_peptide_max_temp;
    if (isotopes_per_peptide_min_ > isotopes_per_peptide_max_)
    {
      swap(isotopes_per_peptide_min_, isotopes_per_peptide_max_);
    }

    rt_typical_ = getParam_().getValue("algorithm:rt_typical");
    rt_band_ = getParam_().getValue("algorithm:rt_band");
    rt_min_ = getParam_().getValue("algorithm:rt_min");
    mz_tolerance_ = getParam_().getValue("algorithm:mz_tolerance");
    mz_unit_ = (getParam_().getValue("algorithm:mz_unit") == "ppm");
    intensity_cutoff_ = getParam_().getValue("algorithm:intensity_cutoff");
    peptide_similarity_ = getParam_().getValue("algorithm:peptide_similarity");
    averagine_similarity_ = getParam_().getValue("algorithm:averagine_similarity");
    averagine_similarity_scaling_ = getParam_().getValue("algorithm:averagine_similarity_scaling");
    missed_cleavages_ = getParam_().getValue("algorithm:missed_cleavages");
    knock_out_ = (getParam_().getValue("algorithm:knock_out") == "true");
    spectrum_type_ = getParam_().getValue("algorithm:spectrum_type");
    averagine_type_ = getParam_().getValue("algorithm:averagine_type");
  }

  /**
   * @brief process parameters of 'labels' section
   */
  void getParameters_labels_()
  {
    Param p = getParam_();
    
    // create map of pairs (label as string, mass shift as double)
    for (Param::ParamIterator it = p.begin(); it != p.end(); ++it)
    {
      label_mass_shift_.insert(make_pair(it->name, it->value));
    }    
  }

  /**
   * @brief split labels string
   *
   * @return list of samples containing lists of corresponding labels
   */
  std::vector<std::vector<String> > splitLabelString_()
  {
    std::vector<std::vector<String> > samples_labels;
    std::vector<String> temp_samples;
    
    String labels(labels_);
    boost::replace_all(labels, "[]", "no_label");
    boost::replace_all(labels, "()", "no_label");
    boost::replace_all(labels, "{}", "no_label");
    boost::split(temp_samples, labels, boost::is_any_of("[](){}")); // any bracket allowed to separate samples
    
    for (unsigned i = 0; i < temp_samples.size(); ++i)
    {
      if (!temp_samples[i].empty())
      {
        if (temp_samples[i]=="no_label")
        {
          vector<String> temp_labels;
          temp_labels.push_back("no_label");
          samples_labels.push_back(temp_labels);
        }
        else
        {
          vector<String> temp_labels;
          boost::split(temp_labels, temp_samples[i], boost::is_any_of(",;: ")); // various separators allowed to separate labels
          samples_labels.push_back(temp_labels);
        }
      }
    }
    
    if (samples_labels.empty())
    {
      vector<String> temp_labels;
      temp_labels.push_back("no_label");
      samples_labels.push_back(temp_labels);
    }

    return samples_labels;
  }

  /**
   * @brief order of charge states
   *
   * 2+ 3+ 4+ 1+ 5+ 6+ ...
   *
   * Order charge states by the likelihood of their occurrence, i.e. we search for the most likely charge states first.
   */
  static size_t order_charge(int charge)
  {
    if ((1 < charge) && (charge < 5))
    {
      return (charge - 1);
    }
    else if (charge == 1)
    {
      return 4;
    }
    else
    {
      return charge;
    }
  }
  
  /**
   * @brief comparator of peak patterns
   *
   * The comperator determines in which order the peak patterns are searched for.
   * First we check the number of mass shifts (triplets before doublets before singlets). 
   * Then we check the first mass shift (for example 6 Da before 12 Da i.e. misscleavage).
   * Finally we check for charges (2+ before 1+, most likely first).
   *
   * @param pattern1    first peak pattern
   * @param pattern2    second peak pattern
   *
   * @return true if pattern1 should be searched before pattern2
   */
  static bool less_pattern(const MultiplexIsotopicPeakPattern& pattern1, const MultiplexIsotopicPeakPattern& pattern2)
  {
    if (pattern1.getMassShiftCount() == pattern2.getMassShiftCount())
    {
      // The first mass shift is by definition always zero.
      if ((pattern1.getMassShiftCount() > 1) && (pattern2.getMassShiftCount() > 1))
      {
        if (pattern1.getMassShiftAt(1) == pattern2.getMassShiftAt(1))
        {
          // 2+ before 3+ before 4+ before 1+ before 5+ before 6+ etc.
          return order_charge(pattern1.getCharge()) < order_charge(pattern2.getCharge());
        }
        else
        {
          return pattern1.getMassShiftAt(1) < pattern2.getMassShiftAt(1);
        }
      }
      else
      {
        // 2+ before 3+ before 4+ before 1+ before 5+ before 6+ etc.
        return order_charge(pattern1.getCharge()) < order_charge(pattern2.getCharge());
      }
    }
    else
    {
      // triplets before doublets before singlets
      return pattern1.getMassShiftCount() > pattern2.getMassShiftCount();
    }
  }

  /**
   * @brief generate list of m/z shifts
   *
   * @param charge_min    minimum charge
   * @param charge_max    maximum charge
   * @param peaks_per_peptide_max    maximum number of isotopes in peptide
   * @param mass_pattern_list    mass shifts due to labelling
   *
   * @return list of m/z shifts
   */
  std::vector<MultiplexIsotopicPeakPattern> generatePeakPatterns_(int charge_min, int charge_max, int peaks_per_peptide_max, std::vector<MultiplexDeltaMasses> mass_pattern_list)
  {
    std::vector<MultiplexIsotopicPeakPattern> list;

    // iterate over all charge states
    for (int c = charge_max; c >= charge_min; --c)
    {
      // iterate over all mass shifts
      for (unsigned i = 0; i < mass_pattern_list.size(); ++i)
      {
        MultiplexIsotopicPeakPattern pattern(c, peaks_per_peptide_max, mass_pattern_list[i], i);
        list.push_back(pattern);
      }
    }
    
    sort(list.begin(), list.end(), less_pattern);
    
    // debug output
    /*for (int i = 0; i < list.size(); ++i)
    {
      std::cout << "charge = " << list[i].getCharge() << "+    shift = " << list[i].getMassShiftAt(1) << " Da\n";
    }*/
    
    return list;
  }
  
  /**
   * @brief calculate peptide intensities
   *
   * @param pattern
   * @param satellites
   *
   * @return vector with intensities for each of the peptides
   */
  std::vector<double> determinePeptideIntensities_(MultiplexIsotopicPeakPattern& pattern, std::multimap<size_t, MultiplexSatellite >& satellites, bool centroided)
  {
    // determine RT shift between the peptides
    // i.e. first determine the RT centre of mass for each peptide
    std::vector<double> rt_peptide;
    std::vector<double> intensity_peptide;
    // loop over peptides
    for (size_t peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      // coordinates of the peptide feature
      // RT is the intensity-average of all satellites peaks of all (!) mass traces
      double rt(0);
      double intensity_sum(0);
      
      // loop over isotopes i.e. mass traces of the peptide
      for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
      {
        // find satellites for this isotope i.e. mass trace
        size_t idx = peptide * isotopes_per_peptide_max_ + isotope;
        std::pair<std::multimap<size_t, MultiplexSatellite >::const_iterator, std::multimap<size_t, MultiplexSatellite >::const_iterator> satellites_isotope;
        satellites_isotope = satellites.equal_range(idx);
        
        // loop over satellites for this isotope i.e. mass trace
        for (std::multimap<size_t, MultiplexSatellite >::const_iterator satellite_it = satellites_isotope.first; satellite_it != satellites_isotope.second; ++satellite_it)
        {
          // find indices of the peak
          size_t rt_idx = (satellite_it->second).getRTidx();
          size_t mz_idx = (satellite_it->second).getMZidx();
          
          // find peak itself
          MSExperiment::ConstIterator it_rt = exp_centroid_.begin();
          std::advance(it_rt, rt_idx);
          MSSpectrum<Peak1D>::ConstIterator it_mz = it_rt->begin();
          std::advance(it_mz, mz_idx);
          
          // get profile vectors
          // In the centroid case, a single entry.
          // In the profile case, all data points from scanning over the peak profile.
          std::vector<double> mz_profile = (satellite_it->second).getMZ();
          std::vector<double> intensity_profile = (satellite_it->second).getIntensity();
          
          // loop over profile data points
          for (std::vector<double>::const_iterator it_intensity = intensity_profile.begin(); it_intensity < intensity_profile.end(); ++it_intensity)
          {
            rt += it_rt->getRT() * (*it_intensity);
            intensity_sum += (*it_intensity);
          }
        }
      }
      
      if (intensity_sum <= 0)
      {
        std::ostringstream strs;
        strs << intensity_sum;
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The total peptide intensity is not positive. Aborting.", strs.str());
      }
      
      rt /= intensity_sum;
      rt_peptide.push_back(rt);
      intensity_peptide.push_back(intensity_sum);
    }
    
    // determine the fold changes between the lightest peptide and the remaining ones
    // TODO Replace the lightest peptide by the highest intensity peptide in the multiplet
    std::vector<double> ratio_peptide;
    ratio_peptide.push_back(1.0);
    // loop over other peptides
    for (size_t peptide = 1; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      // fill the vectors from which the Pearson correlation for the ratio i.e. fold change will be calculated
      std::vector<double> intensities_light;
      std::vector<double> intensities_other;
      // loop over isotopes i.e. mass traces of the peptide
      for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
      {
        // find satellites for this isotope in the light peptide
        std::pair<std::multimap<size_t, MultiplexSatellite >::const_iterator, std::multimap<size_t, MultiplexSatellite >::const_iterator> satellites_isotope_1;
        satellites_isotope_1 = satellites.equal_range(isotope);
        
        // find satellites for this isotope in the second peptide
        std::pair<std::multimap<size_t, MultiplexSatellite >::const_iterator, std::multimap<size_t, MultiplexSatellite >::const_iterator> satellites_isotope_2;
        satellites_isotope_2 = satellites.equal_range(peptide * isotopes_per_peptide_max_ + isotope);
        
        // if the satellite set for either isotope is empty, we can move on
        if ((satellites_isotope_1.first == satellites_isotope_1.second) || (satellites_isotope_2.first == satellites_isotope_2.second))
        {
          continue;
        }
        
        // loop over satellites for this isotope in the light peptide
        for (std::multimap<size_t, MultiplexSatellite >::const_iterator satellite_it_1 = satellites_isotope_1.first; satellite_it_1 != satellites_isotope_1.second; ++satellite_it_1)
        {
          // find indices of the peak
          size_t rt_idx_1 = (satellite_it_1->second).getRTidx();
          size_t mz_idx_1 = (satellite_it_1->second).getMZidx();
          
          // find peak itself
          MSExperiment::ConstIterator it_rt_1 = exp_centroid_.begin();
          std::advance(it_rt_1, rt_idx_1);
          MSSpectrum<Peak1D>::ConstIterator it_mz_1 = it_rt_1->begin();
          std::advance(it_mz_1, mz_idx_1);
          
          // get profile vectors
          // In the centroid case, a single entry.
          // In the profile case, all data points from scanning over the peak profile.
          std::vector<double> mz_profile_1 = (satellite_it_1->second).getMZ();
          std::vector<double> intensity_profile_1 = (satellite_it_1->second).getIntensity();
          
          // find corresponding spectrum
          double mz_1 = it_mz_1->getMZ();
          double rt_1 = it_rt_1->getRT();
          double rt_2_target = rt_1 + rt_peptide[peptide] - rt_peptide[0];
          
          // find the spectra which bracket <rt_2_target>
          int rt_idx_2_before = -1;
          int rt_idx_2_after = -1;
          double rt_2_before = -1.0;
          double rt_2_after = -1.0;
          double mz_2_before = -1.0;
          double mz_2_after = -1.0;
          std::vector<double> mz_profile_2_before;
          std::vector<double> mz_profile_2_after;
          std::vector<double> intensity_profile_2_before;
          std::vector<double> intensity_profile_2_after;
          // loop over satellites for this isotope in the second peptide
          for (std::multimap<size_t, MultiplexSatellite >::const_iterator satellite_it_2 = satellites_isotope_2.first; satellite_it_2 != satellites_isotope_2.second; ++satellite_it_2)
          {
            // find RT of the peak
            size_t rt_idx_temp = (satellite_it_2->second).getRTidx();
            size_t mz_idx_temp = (satellite_it_2->second).getMZidx();
            MSExperiment::ConstIterator it_rt_temp = exp_centroid_.begin();
            std::advance(it_rt_temp, rt_idx_temp);
            MSSpectrum<Peak1D>::ConstIterator it_mz_temp = it_rt_temp->begin();
            std::advance(it_mz_temp, mz_idx_temp);
            double rt_temp = it_rt_temp->getRT();
            double mz_temp = it_mz_temp->getMZ();
            
            // a better rt_2_before
            if ((rt_temp <= rt_2_target) && ((rt_2_before < rt_temp) || (rt_2_before == -1.0)))
            {
              rt_idx_2_before = rt_idx_temp;
              rt_2_before = rt_temp;
              mz_2_before = mz_temp;
              mz_profile_2_before = (satellite_it_2->second).getMZ();
              intensity_profile_2_before = (satellite_it_2->second).getMZ();
            }
            
            // a better rt_2_after
            if ((rt_temp >= rt_2_target) && ((rt_temp < rt_2_after) || (rt_2_after == -1.0)))
            {
              rt_idx_2_after = rt_idx_temp;
              rt_2_after = rt_temp;
              mz_2_after = mz_temp;
              mz_profile_2_after = (satellite_it_2->second).getMZ();
              intensity_profile_2_after = (satellite_it_2->second).getMZ();
            }
          }
          
          // No lower bracket found. (No satellites for the exact RT shift. But let us use this satellite as next best fit.)
          if ((rt_2_before == -1.0) && (rt_2_after >= 0.0))
          {
            rt_idx_2_before = rt_idx_2_after;
            rt_2_before = rt_2_after;
            mz_2_before = mz_2_after;
            mz_profile_2_before = mz_profile_2_after;
            intensity_profile_2_before = intensity_profile_2_after;
          }
          
          // No upper bracket found. (No satellites for the exact RT shift. But let us use this satellite as next best fit.)
          if ((rt_2_after == -1.0) && (rt_2_before >= 0.0))
          {
            rt_idx_2_after = rt_idx_2_before;
            rt_2_after = rt_2_before;
            mz_2_after = mz_2_before;
            mz_profile_2_after = mz_profile_2_before;
            intensity_profile_2_after = intensity_profile_2_before;
          }
          
          // fill intensity vector (centroid mode)
          if (centroided)
          {
            intensities_light.push_back(intensity_profile_1[0]);
            if (rt_idx_2_before == rt_idx_2_after)
            {
              intensities_other.push_back(intensity_profile_2_before[0]);
            }
            else
            {
              // There is no spectrum at rt_2_target. So we interpolate linearly between spactra.
              intensities_other.push_back(intensity_profile_2_before[0] + (intensity_profile_2_after[0] - intensity_profile_2_before[0]) * (rt_2_target - rt_2_before) / (rt_2_after - rt_2_before));
            }
          }
          // fill intensity vector (profile mode)
          else
          {
            std::cout << "mz_profile_1 size = " << mz_profile_1.size() << "    intensity_profile_1 size = " << intensity_profile_1.size() << "\n";
            
            // loop over m/z around the first satellite
            for (int i = 0; i < mz_profile_1.size(); ++i)
            {
              std::cout << "    m/z = " << mz_profile_1[i] << "    m/z difference = " << (mz_profile_1[i] - mz_1) << "    intensity = " << intensity_profile_1[i] << "\n";
              
              double intensity_before_temp = navigators_[rt_idx_2_before].eval(mz_profile_1[i] - mz_1 + mz_2_before);
              double intensity_after_temp = navigators_[rt_idx_2_after].eval(mz_profile_1[i] - mz_1 + mz_2_after);
              
              if ((intensity_before_temp > intensity_cutoff_) && (intensity_after_temp > intensity_cutoff_))
              {
                if (rt_idx_2_before == rt_idx_2_after)
                {
                  intensities_other.push_back(intensity_before_temp);
                }
                else
                {
                  // There is no spectrum at rt_2_target. So we interpolate linearly between spactra.
                  intensities_other.push_back(intensity_before_temp + (intensity_after_temp - intensity_before_temp) * (rt_2_target - rt_2_before) / (rt_2_after - rt_2_before));
                }
              }
            }
          }
          
          
          
          
          
          
          
          
          
          
          // loop over satellites for this isotope in the second peptide
          /*double rt_earlier = -1;
          std::vector<double> intensity_earlier(0);
          double rt_later = -1;
          std::vector<double> intensity_later(0);
          for (std::multimap<size_t, MultiplexSatellite >::const_iterator satellite_it_2 = satellites_isotope_2.first; satellite_it_2 != satellites_isotope_2.second; ++satellite_it_2)
          {
            // find indices of the peak
            size_t rt_idx_2 = (satellite_it_2->second).getRTidx();
            size_t mz_idx_2 = (satellite_it_2->second).getMZidx();
            
            // find peak itself
            MSExperiment::ConstIterator it_rt_2 = exp_centroid_.begin();
            std::advance(it_rt_2, rt_idx_2);
            MSSpectrum<Peak1D>::ConstIterator it_mz_2 = it_rt_2->begin();
            std::advance(it_mz_2, mz_idx_2);
 
            // get profile vectors
            // In the centroid case, a single entry.
            // In the profile case, all data points from scanning over the peak profile. Note in the profile case the vectors can be of any size including empty.
            std::vector<double> mz_profile_2 = (satellite_it_2->second).getMZ();
            std::vector<double> intensity_profile_2 = (satellite_it_2->second).getIntensity();
           
            if (it_rt_2->getRT() <= rt_2 && (std::abs(it_rt_2->getRT() - rt_2) < std::abs(rt_earlier - rt_2)))
            {
              rt_earlier = it_rt_2->getRT();
              intensity_earlier = intensity_profile_2;
            }
            
            if (it_rt_2->getRT() >= rt_2 && (std::abs(it_rt_2->getRT() - rt_2) < std::abs(rt_later - rt_2)))
            {
              rt_later = it_rt_2->getRT();
              intensity_later = intensity_profile_2;
            }
            
          }
          
          if (rt_earlier == -1 || rt_later == -1)
          {
            continue;
          }*/
          
          
          
          
          // DEBUG OUTPUT
          //std::cout << "size intensity_profile_1 = " << intensity_profile_1.size() << "    size intensity_earlier = " << intensity_earlier.size() << "    size intensity_later = " << intensity_later.size() << "\n";
          
          /*if (intensity_profile_1.size() == intensity_earlier.size())
          {
            //std::cout << "SAME.    size (profile 1) = " << intensity_profile_1.size() << "    size (earlier) = " << intensity_earlier.size() << "\n";
          }
          else
          {
            //std::cout << "    DIFFERENT.    size (profile 1) = " << intensity_profile_1.size() << "    size (earlier) = " << intensity_earlier.size() << "\n";
            std::cout << "    DIFFERENT.    difference = " << std::abs((int) intensity_profile_1.size() - (int) intensity_earlier.size()) << "\n";
          }*/
          
          
          
          
          
          // Our target lies on or between two satellites of the 'other' peptide.
          /*if ((rt_earlier > 0) && (rt_later > 0))
          {
            // linearly interpolated intensity
            // TODO: Maybe linear interpolation between spectra does more bad than good. Perhaps simply pick the nearest spectrum?
            std::vector<double> intensity_other;
            if ((rt_2 == rt_earlier) || (rt_later == rt_earlier))
            {
              intensity_other = intensity_earlier;
            }
            else
            {
              // loop simultaneously over earlier and later intensity vectors
              std::vector<double>::const_iterator it_earlier;
              std::vector<double>::const_iterator it_later;
              for (it_earlier = intensity_earlier.begin(), it_later = intensity_later.begin();
                   it_earlier != intensity_earlier.end(), it_later != intensity_later.end();
                   ++it_earlier, ++it_later)
              {
                intensity_other.push_back(*it_earlier + (*it_later - *it_earlier)*(rt_2 - rt_earlier)/(rt_later - rt_earlier));
              }
            }
            
            intensities_light.insert(intensities_light.end(), intensity_profile_1.begin(), intensity_profile_1.end());
            intensities_other.insert(intensities_other.end(), intensity_other.begin(), intensity_other.end());
          }*/
        
          // END OF THE LINE // END OF THE LINE // END OF THE LINE // END OF THE LINE // END OF THE LINE // END OF THE LINE // 
        }
        
      }
      
      // If less than three matches are found, we cannot reliably calculate the intensity ratio (aka slope) and report the uncorrected intensities.
      if ((intensities_light.size() < 3) || (intensities_other.size() < 3))
      {
        return intensity_peptide;
      }
      
      /*if (intensities_light.size() != intensities_other.size())
      {
        std::cout << "Vector size are not equal.      intensities_light.size() = " << intensities_light.size() << "  intensities_other.size() = " << intensities_other.size() << "\n";
      }*/
      
      // determine ratios through linear regression of all corresponding intensities
      LinearRegressionWithoutIntercept linreg;
      linreg.addData(intensities_light, intensities_other);
      double slope = linreg.getSlope();
      
      // DEBUG OUTPUT
      /*if (slope < 0)
      {
        std::cout << "slope = " << slope << "    size (intensity light) = " << intensities_light.size() << "    size (intensity other) = " << intensities_other.size() << "\n";
      }*/
      
      ratio_peptide.push_back(slope);
    }
    
    // correct peptide intensities
    // The peptide ratios are calculated as linear regression of (spline-interpolated) profile intensities, @see linreg
    // The individual peptide intensities are the sum of the same profile intensities. But the quotient of these peptide intensities
    // is not necessarily the same as the independently calculated ratio from the linear regression. Since the peptide ratio
    // from linear regression is the more accurate one, we correct the two peptide intensities by projecting them onto the ratio.
    // In the end, both peptide ratio from linear regression and the quotient of the peptide intensities are identical.
    std::vector<double> intensity_peptide_corrected;
    if (intensity_peptide.size() == 2)
    {
      double intensity1 = (intensity_peptide[0] + ratio_peptide[1] * intensity_peptide[1]) / (1 + ratio_peptide[1] * ratio_peptide[1]);
      double intensity2 = ratio_peptide[1] * intensity1;
      
      // DEBUG OUTPUT
      /*if (intensity1 <= 0)
      {
        std::cout << "peptide intensity 1 = " << intensity1 << "    intensity_peptide[0] = " << intensity_peptide[0] << "    intensity_peptide[1] = " << intensity_peptide[1] << "    ratio_peptide[1] = " << ratio_peptide[1] << "\n";
      }
      if (intensity2 <= 0)
      {
        std::cout << "peptide intensity 2 = " << intensity2 << "    ratio_peptide[1] = " << ratio_peptide[1] << "\n";
      }*/
      
      intensity_peptide_corrected.push_back(intensity1);
      intensity_peptide_corrected.push_back(intensity2);
    }
    else if (intensity_peptide.size() > 2)
    {
      // Now with n instead of two peptide intensities, one needs to project the peptide intensities onto the hyperplane defined
      // by the set of all peptide ratios (TODO). Instead, it is simpler to keep the lightest peptide intensity fixed, and correct
      // only the remaining ones. The correct peptide ratio (from linear regression) is reported on both cases.
      intensity_peptide_corrected.push_back(intensity_peptide[0]);
      for (unsigned i = 1; i < intensity_peptide.size(); ++i)
      {
        intensity_peptide_corrected.push_back(ratio_peptide[i] * intensity_peptide[0]);
      }
    }
    else
    {
      // For simple feature detection (singlets) the intensities remain unchanged.
      intensity_peptide_corrected.push_back(intensity_peptide[0]);
    }
    
    return intensity_peptide_corrected;
  }
  
  /**
   * @brief generates consensus and feature maps containing all peptide multiplets
   *
   * @param patterns    patterns of isotopic peaks we have been searching for
   * @param filter_results    filter results for each of the patterns
   * @param cluster_results    clusters of filter results
   * @param consensus_map    consensus map with peptide multiplets (to be filled)
   * @param feature_map    feature map with peptides (to be filled)
   */
  
  
  void generateMaps_(std::vector<MultiplexIsotopicPeakPattern> patterns, std::vector<MultiplexFilteredMSExperiment> filter_results, std::vector<std::map<int, GridBasedCluster> > cluster_results, ConsensusMap& consensus_map, FeatureMap& feature_map, bool centroided)
  {
    // loop over peak patterns
    for (unsigned pattern = 0; pattern < patterns.size(); ++pattern)
    {
      // loop over clusters
      size_t cluster_idx(0);
      for (std::map<int, GridBasedCluster>::const_iterator cluster_it = cluster_results[pattern].begin(); cluster_it != cluster_results[pattern].end(); ++cluster_it)
      {
        GridBasedCluster cluster = cluster_it->second;
        std::vector<int> points = cluster.getPoints();
        
        // Construct a satellite set for the complete peptide multiplet
        // Make sure there are no duplicates, i.e. the same satellite from different filtered peaks.
        std::multimap<size_t, MultiplexSatellite > satellites;
        // loop over points in cluster
        for (std::vector<int>::const_iterator point_it = points.begin(); point_it != points.end(); ++point_it)
        {
          MultiplexFilteredPeak peak = filter_results[pattern].getPeak(*point_it);
          // loop over satellites of the peak
          for (std::multimap<size_t, MultiplexSatellite >::const_iterator satellite_it = peak.getSatellites().begin(); satellite_it != peak.getSatellites().end(); ++satellite_it)
          {
            // check if this satellite (i.e. these indices) are already in the set
            bool satellite_in_set = false;
            for (std::multimap<size_t, MultiplexSatellite >::const_iterator satellite_it_2 = satellites.begin(); satellite_it_2 != satellites.end(); ++satellite_it_2)
            {
              if ((satellite_it_2->second.getRTidx() == satellite_it->second.getRTidx()) && (satellite_it_2->second.getMZidx() == satellite_it->second.getMZidx()))
              {
                satellite_in_set = true;
                continue;
              }
            }
            if (satellite_in_set)
            {
              continue;
            }
            
            satellites.insert(std::make_pair(satellite_it->first, MultiplexSatellite(satellite_it->second.getRTidx(), satellite_it->second.getMZidx(), satellite_it->second.getMZ(), satellite_it->second.getIntensity())));
          }
        }
        
        // determine peptide intensities
        std::vector<double> peptide_intensities = determinePeptideIntensities_(patterns[pattern], satellites, centroided);
        
        // If no reliable peptide intensity can be determined, we do not report the peptide multiplet.
        if (peptide_intensities[0] == -1)
        {
          continue;
        }
        
        std::vector<Feature> features;
        ConsensusFeature consensus;
        bool abort = false;
        
        // construct the feature and consensus maps
        // loop over peptides
        for (size_t peptide = 0; (peptide < patterns[pattern].getMassShiftCount() && !abort); ++peptide)
        {
          // coordinates of the peptide feature
          // RT is the intensity-average of all satellites peaks of the mono-isotopic mass trace
          // m/z is the intensity-average of all satellites peaks of the mono-isotopic mass trace
          Feature feature;
          double rt(0);
          double mz(0);
          double intensity_sum(0);
          
          // loop over isotopes i.e. mass traces of the peptide
          for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
          {
            // find satellites for this isotope i.e. mass trace
            size_t idx = peptide * isotopes_per_peptide_max_ + isotope;
            std::pair<std::multimap<size_t, MultiplexSatellite >::const_iterator, std::multimap<size_t, MultiplexSatellite >::const_iterator> satellites_isotope;
            satellites_isotope = satellites.equal_range(idx);
            
            DBoundingBox<2> mass_trace;
            
            // loop over satellites for this isotope i.e. mass trace
            for (std::multimap<size_t, MultiplexSatellite >::const_iterator satellite_it = satellites_isotope.first; satellite_it != satellites_isotope.second; ++satellite_it)
            {
              // find indices of the peak
              size_t rt_idx = (satellite_it->second).getRTidx();
              size_t mz_idx = (satellite_it->second).getMZidx();
              
              // find peak itself
              MSExperiment::ConstIterator it_rt = exp_centroid_.begin();
              std::advance(it_rt, rt_idx);
              MSSpectrum<Peak1D>::ConstIterator it_mz = it_rt->begin();
              std::advance(it_mz, mz_idx);
              
              // get profile vectors
              // In the centroid case, a single entry.
              // In the profile case, all data points from scanning over the peak profile.
              std::vector<double> mz_profile = (satellite_it->second).getMZ();
              std::vector<double> intensity_profile = (satellite_it->second).getIntensity();
              
              // loop over profile data points
              std::vector<double>::const_iterator it_mz_profile;
              std::vector<double>::const_iterator it_intensity_profile;
              for (it_mz_profile = mz_profile.begin(), it_intensity_profile = intensity_profile.begin();
                   it_mz_profile != mz_profile.end(), it_intensity_profile != intensity_profile.end();
                   ++it_mz_profile, ++it_intensity_profile)
              {
                if (isotope == 0)
                {
                  rt += it_rt->getRT() * (*it_intensity_profile);
                  mz += it_mz->getMZ() * (*it_intensity_profile);
                  intensity_sum += (*it_intensity_profile);
                }
                
                mass_trace.enlarge(it_rt->getRT(), (*it_mz_profile));
                
              }
   
            }
            
            if ((mass_trace.width() == 0) || (mass_trace.height() == 0))
            {
              // The mass trace contains only a single point. Add a small margin around
              // the point, otherwise the mass trace is considered empty and not drawn.
              // TODO: Remove the magic number for the margin.
              mass_trace.enlarge(mass_trace.minX() - 0.01, mass_trace.minY() - 0.01);
              mass_trace.enlarge(mass_trace.maxX() + 0.01, mass_trace.maxY() + 0.01);
            }
            
            if (!(mass_trace.isEmpty()))
            {
              ConvexHull2D hull;
              hull.addPoint(DPosition<2>(mass_trace.minX(), mass_trace.minY()));
              hull.addPoint(DPosition<2>(mass_trace.minX(), mass_trace.maxY()));
              hull.addPoint(DPosition<2>(mass_trace.maxX(), mass_trace.minY()));
              hull.addPoint(DPosition<2>(mass_trace.maxX(), mass_trace.maxY()));
              feature.getConvexHulls().push_back(hull);
            }
          }
          
          if (intensity_sum <= 0)
          {
            std::ostringstream strs;
            strs << intensity_sum;
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The total intensity of the mono-isotopic mass trace is non-positive. Aborting.", strs.str());
          }
          
          rt /= intensity_sum;
          mz /= intensity_sum;
          
          feature.setRT(rt);
          feature.setMZ(mz);
          feature.setIntensity(peptide_intensities[peptide]);
          feature.setCharge(patterns[pattern].getCharge());
          feature.setOverallQuality(1.0);
          
          // Check that the feature eluted long enough.
          DBoundingBox<2> box = feature.getConvexHull().getBoundingBox();
          if (box.maxX() - box.minX() < rt_min_)
          {
            abort = true;
            continue;
          }
          
          features.push_back(feature);
          
          if (peptide == 0)
          {
            // The first/lightest peptide acts as anchor of the peptide multiplet consensus.
            // All peptide feature handles are connected to this point.
            consensus.setRT(rt);
            consensus.setMZ(mz);
            consensus.setIntensity(peptide_intensities[peptide]);
            consensus.setCharge(patterns[pattern].getCharge());
            consensus.setQuality(1.0);
          }
          
          FeatureHandle feature_handle;
          feature_handle.setRT(rt);
          feature_handle.setMZ(mz);
          feature_handle.setIntensity(peptide_intensities[peptide]);
          feature_handle.setCharge(patterns[pattern].getCharge());
          feature_handle.setMapIndex(peptide);
          //feature_handle.setUniqueId(&UniqueIdInterface::setUniqueId);    // TODO: Do we need to set unique ID?
          consensus.insert(feature_handle);
          consensus_map.getFileDescriptions()[peptide].size++;
        }
        
        if (!abort)
        {
          consensus_map.push_back(consensus);
          for(std::vector<Feature>::iterator it = features.begin(); it != features.end(); ++it)
          {
            feature_map.push_back(*it);
          }
        }

      ++cluster_idx;        
      }
      
    }
    
  }
  
  /**
   * @brief generates the data structure for mzQuantML output
   *
   * @param exp    experimental data
   * @param consensus_map    consensus map with complete quantitative information
   * @param quantifications    MSQuantifications data structure for writing mzQuantML (mzq)
   */
  void generateMSQuantifications(MSExperiment& exp, ConsensusMap& consensus_map, MSQuantifications& quantifications)
  {
    // generate the labels
    // (for each sample a list of (label string, mass shift) pairs)
    // for example triple-SILAC: [(none,0)][(Lys4,4.0251),(Arg6,6.0201)][Lys8,8.0141)(Arg10,10.0082)]
    std::vector<std::vector<std::pair<String, double> > > labels;
    
    for (unsigned sample = 0; sample < samples_labels_.size(); ++sample)
    {
      // The labels are required to be ordered in mass shift.
      std::map<double, String> single_label_map;
      std::vector<std::pair<String, double> > single_label;
      for (unsigned label = 0; label < samples_labels_[sample].size(); ++label)
      {
        String label_string = samples_labels_[sample][label];
        double shift;
        if (label_string == "")
        {
          label_string = "none";
          shift = 0;
        }
        else
        {
          shift = label_mass_shift_[label_string];
        }

        single_label_map[shift] = label_string;
      }
      for (std::map<double, String>::const_iterator it = single_label_map.begin(); it != single_label_map.end(); ++it)
      {
        std::pair<String, double> label_shift(it->second, it->first);
        single_label.push_back(label_shift);
      }
      labels.push_back(single_label);
    }

    quantifications.registerExperiment(exp, labels);
    quantifications.assignUIDs();

    MSQuantifications::QUANT_TYPES quant_type = MSQuantifications::MS1LABEL;
    quantifications.setAnalysisSummaryQuantType(quant_type);

    // add results from  analysis
    LOG_DEBUG << "Generating output mzQuantML file..." << endl;
    ConsensusMap numap(consensus_map);
    
    //calculate ratios
    for (ConsensusMap::iterator cit = numap.begin(); cit != numap.end(); ++cit)
    {
      // make ratio templates
      std::vector<ConsensusFeature::Ratio> rts;
      for (std::vector<MSQuantifications::Assay>::const_iterator ait = quantifications.getAssays().begin() + 1; ait != quantifications.getAssays().end(); ++ait)
      {
        ConsensusFeature::Ratio r;
        r.numerator_ref_ = String(quantifications.getAssays().begin()->uid_);
        r.denominator_ref_ = String(ait->uid_);
        r.description_.push_back("Simple ratio calc");
        r.description_.push_back("light to medium/.../heavy");
        rts.push_back(r);
      }

      const ConsensusFeature::HandleSetType& feature_handles = cit->getFeatures();
      if (feature_handles.size() > 1)
      {
        std::set<FeatureHandle, FeatureHandle::IndexLess>::const_iterator fit = feature_handles.begin(); // this is unlabeled
        ++fit;
        for (; fit != feature_handles.end(); ++fit)
        {
          Size ri = std::distance(feature_handles.begin(), fit);
          rts[ri - 1].ratio_value_ =  feature_handles.begin()->getIntensity() / fit->getIntensity(); // a proper algo should never have 0-intensities so no 0devison ...
        }
      }

      cit->setRatios(rts);
    }
    quantifications.addConsensusMap(numap); //add FeatureFinderMultiplex result

  }

  /**
   * @brief Write consensus map to consensusXML file.
   *
   * @param filename    name of consensusXML file
   * @param map    consensus map for output
   */
  void writeConsensusMap_(const String& filename, ConsensusMap& map) const
  {
    map.sortByPosition();
    map.applyMemberFunction(&UniqueIdInterface::setUniqueId);
    map.setExperimentType("labeled_MS1");

    // annotate maps
    for (unsigned i = 0; i < samples_labels_.size(); ++i)
    {
      ConsensusMap::FileDescription& desc = map.getFileDescriptions()[i];
      desc.filename = filename;

      if (knock_out_)
      {
        // With knock-outs present, the correct labels can only be determined during ID mapping.
        // For now, we simply store a unique identifier.
        std::stringstream stream;
        stream << "label " << i;
        desc.label = stream.str();
      }
      else
      {
        String label_string;
        for (unsigned j = 0; j < samples_labels_[i].size(); ++j)
        {
          label_string.append(samples_labels_[i][j]);
        }
        desc.label = label_string;
      }
    }

    ConsensusXMLFile file;
    file.store(filename, map);
  }

  /**
   * @brief Write feature map to featureXML file.
   *
   * @param filename    name of featureXML file
   * @param map    feature map for output
   */
  void writeFeatureMap_(const String& filename, FeatureMap& map) const
  {
    map.sortByPosition();
    map.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    FeatureXMLFile file;
    file.store(filename, map);
  }

  /**
   * @brief Write MS quantification map to mzq file.
   *
   * @param filename    name of mzq file
   * @param map    MS quantification map for output
   */
  void writeMSQuantifications(const String& filename, MSQuantifications& msq) const
  {
    MzQuantMLFile file;
    file.store(filename, msq);
  }

  /**
  * @brief simple linear regression through the origin
  *
  * TODO: combine with OpenMS/MATH/STATISTICS/LinearRegression
  */
  class LinearRegressionWithoutIntercept
  {
public:
    /**
     * @brief constructor
     */
    LinearRegressionWithoutIntercept() :
      sum_xx_(0), sum_xy_(0), n_(0)
    {
    }

    /**
     * @brief adds an observation (x,y) to the regression data set.
     *
     * @param x    independent variable value
     * @param y    dependent variable value
     */
    void addData(double x, double y)
    {
      sum_xx_ += x * x;
      sum_xy_ += x * y;

      ++n_;
    }

    /**
     * @brief adds observations (x,y) to the regression data set.
     *
     * @param x    vector of independent variable values
     * @param y    vector of dependent variable values
     */
    void addData(std::vector<double>& x, std::vector<double>& y)
    {
      for (unsigned i = 0; i < x.size(); ++i)
      {
        addData(x[i], y[i]);
      }
    }

    /**
     * @brief returns the slope of the estimated regression line.
     */
    double getSlope()
    {
      if (n_ < 2)
      {
        return std::numeric_limits<double>::quiet_NaN(); // not enough data
      }
      
      return sum_xy_ / sum_xx_;
    }

private:
    /**
     * @brief total variation in x
     */
    double sum_xx_;

    /**
     * @brief sum of products
     */
    double sum_xy_;

    /**
     * @brief number of observations
     */
    int n_;

  };


  ExitCodes main_(int, const char**)
  {

    /**
     * handle parameters
     */
    getParameters_in_out_();
    getParameters_labels_();
    getParameters_algorithm_();
    
    /**
     * load input
     */
    MzMLFile file;
    MSExperiment exp;

    // only read MS1 spectra
    std::vector<int> levels;
    levels.push_back(1);
    file.getOptions().setMSLevels(levels);

    LOG_DEBUG << "Loading input..." << endl;
    file.setLogType(log_type_);
    file.load(in_, exp);

    if (exp.getSpectra().empty())
    {
      throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__, "Error: No MS1 spectra in input file.");
    }

    // update m/z and RT ranges
    exp.updateRanges();

    // sort according to RT and MZ
    exp.sortSpectra();

    // determine type of spectral data (profile or centroided)
    SpectrumSettings::SpectrumType spectrum_type = exp[0].getType();
    if (spectrum_type == SpectrumSettings::UNKNOWN)
    {
      spectrum_type = PeakTypeEstimator().estimateType(exp[0].begin(), exp[0].end());
    }
    
    bool centroided;
    if (spectrum_type_=="automatic")
    {
      centroided = spectrum_type == SpectrumSettings::PEAKS;
    }
    else if (spectrum_type_=="centroid")
    {
      centroided = true;
    }
    else  // "profile"
    {
      centroided = false;
    }
    
    if (centroided)
    {
      exp_centroid_ = exp;
    }
    else
    {
      exp_profile_ = exp;
      
      // spline interpolate the profile data
      for (MSExperiment::Iterator it = exp_profile_.begin(); it < exp_profile_.end(); ++it)
      {
        SplineSpectrum spline(*it);
        exp_spline_profile_.push_back(spline);
        navigators_.push_back(spline.getNavigator());
      }
    }
    // TODO: Explicitly destruct <exp>?

    /**
     * pick peaks
     */
    //MSExperiment exp_picked;
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_s; // peak boundaries for spectra
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_c; // peak boundaries for chromatograms

    if (!centroided)
    {
      PeakPickerHiRes picker;
      Param param = picker.getParameters();
      picker.setLogType(log_type_);
      param.setValue("ms_levels", ListUtils::create<Int>("1"));
      param.setValue("signal_to_noise", 0.0); // signal-to-noise estimation switched off
      picker.setParameters(param);

      picker.pickExperiment(exp_profile_, exp_centroid_, boundaries_exp_s, boundaries_exp_c);
    }

    /**
     * filter for peak patterns
     */
    MultiplexDeltaMassesGenerator generator = MultiplexDeltaMassesGenerator(labels_, missed_cleavages_, label_mass_shift_);
    if (knock_out_)
    {
      generator.generateKnockoutDeltaMasses();
    }
    generator.printSamplesLabelsList();
    generator.printDeltaMassesList();
    
    std::vector<MultiplexDeltaMasses> masses = generator.getDeltaMassesList();
    std::vector<MultiplexIsotopicPeakPattern> patterns = generatePeakPatterns_(charge_min_, charge_max_, isotopes_per_peptide_max_, masses);

    std::vector<MultiplexFilteredMSExperiment> filter_results;
    if (centroided)
    {
      // centroided data
      MultiplexFilteringCentroided filtering(exp_centroid_, patterns, isotopes_per_peptide_min_, isotopes_per_peptide_max_, intensity_cutoff_, rt_band_, mz_tolerance_, mz_unit_, peptide_similarity_, averagine_similarity_, averagine_similarity_scaling_, averagine_type_);
      filtering.setLogType(log_type_);
      filter_results = filtering.filter();
    }
    else
    {
      // profile data
      MultiplexFilteringProfile filtering(exp_profile_, exp_centroid_, boundaries_exp_s, patterns, isotopes_per_peptide_min_, isotopes_per_peptide_max_, intensity_cutoff_, rt_band_, mz_tolerance_, mz_unit_, peptide_similarity_, averagine_similarity_, averagine_similarity_scaling_, averagine_type_);
      filtering.setLogType(log_type_);
      filter_results = filtering.filter();
    }

    /**
     * cluster filter results
     */
    std::vector<std::map<int, GridBasedCluster> > cluster_results;
    if (centroided)
    {
      // centroided data
      MultiplexClustering clustering(exp_centroid_, mz_tolerance_, mz_unit_, rt_typical_, rt_min_);
      clustering.setLogType(log_type_);
      cluster_results = clustering.cluster(filter_results);
    }
    else
    {
      // profile data
      MultiplexClustering clustering(exp_profile_, exp_centroid_, boundaries_exp_s, rt_typical_, rt_min_);
      clustering.setLogType(log_type_);
      cluster_results = clustering.cluster(filter_results);
    }

    /**
     * write to output
     */
    ConsensusMap consensus_map;
    FeatureMap feature_map;
    
    if (centroided)
    {
      consensus_map.setPrimaryMSRunPath(exp_centroid_.getPrimaryMSRunPath());
      feature_map.setPrimaryMSRunPath(exp_centroid_.getPrimaryMSRunPath());
    }
    else
    {
      consensus_map.setPrimaryMSRunPath(exp_profile_.getPrimaryMSRunPath());
      feature_map.setPrimaryMSRunPath(exp_profile_.getPrimaryMSRunPath());
    }
    
    generateMaps_(patterns, filter_results, cluster_results, consensus_map, feature_map, centroided);
    
    if (out_ != "")
    {
      writeConsensusMap_(out_, consensus_map);
    }
    if (out_features_ != "")
    {
      writeFeatureMap_(out_features_, feature_map);
    }

    /*if (out_mzq_ != "")
    {
      MSQuantifications quantifications;
      generateMSQuantifications(exp, consensus_map, quantifications);
      writeMSQuantifications(out_mzq_, quantifications);
    }*/

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPFeatureFinderMultiplex tool;
  return tool.main(argc, argv);
}

//@endcond
