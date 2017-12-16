// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Douglas McCloskey $
// $Authors: Douglas McCloskey $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>

//Analysis classes
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(AbsoluteQuantitationMethod, "$Id$")

/////////////////////////////////////////////////////////////

AbsoluteQuantitationMethod* ptr = 0;
AbsoluteQuantitationMethod* nullPointer = 0;
START_SECTION((AbsoluteQuantitationMethod()))
	ptr = new AbsoluteQuantitationMethod();
	TEST_NOT_EQUAL(ptr, nullPointer);
END_SECTION

START_SECTION((~AbsoluteQuantitationMethod()))
	delete ptr;
END_SECTION

START_SECTION((bool checkLOD(const double & value)))

  AbsoluteQuantitationMethod aqm;
  double value = 2.0;

  // tests
  aqm.setLLOD(0.0);
  aqm.setULOD(4.0);
  TEST_EQUAL(aqm.checkLOD(value),true);
  aqm.setLLOD(0.0);
  aqm.setULOD(1.0);
  TEST_EQUAL(aqm.checkLOD(value),false);
  aqm.setLLOD(3.0);
  aqm.setULOD(4.0);
  TEST_EQUAL(aqm.checkLOD(value),false);
END_SECTION

START_SECTION((bool checkLOQ(const double & value)))

  AbsoluteQuantitationMethod aqm;
  double value = 2.0;

  // tests
  aqm.setLLOQ(0.0);
  aqm.setULOQ(4.0);
  TEST_EQUAL(aqm.checkLOQ(value),true);
  aqm.setLLOQ(0.0);
  aqm.setULOQ(1.0);
  TEST_EQUAL(aqm.checkLOQ(value),false);
  aqm.setLLOQ(3.0);
  aqm.setULOQ(4.0);
  TEST_EQUAL(aqm.checkLOQ(value),false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST