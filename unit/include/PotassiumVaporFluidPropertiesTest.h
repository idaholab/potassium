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
#include "PotassiumVaporFluidProperties.h"

class PotassiumVaporFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  PotassiumVaporFluidPropertiesTest() : MooseObjectUnitTest("PotassiumApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("PotassiumVaporFluidProperties");
    _fe_problem->addUserObject("PotassiumVaporFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<PotassiumVaporFluidProperties>("fp");
  }

  const PotassiumVaporFluidProperties * _fp;
};
