//* This file is part of potassium
//* https://github.com/idaholab/potassium
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/potassium/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseObjectUnitTest.h"
#include "PotassiumTwoPhaseFluidProperties.h"

class PotassiumTwoPhaseFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  PotassiumTwoPhaseFluidPropertiesTest() : MooseObjectUnitTest("PotassiumApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("PotassiumTwoPhaseFluidProperties");
    _fe_problem->addUserObject("PotassiumTwoPhaseFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<PotassiumTwoPhaseFluidProperties>("fp");
  }

  const PotassiumTwoPhaseFluidProperties * _fp;
};
