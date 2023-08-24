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
#include "PotassiumLiquidFluidProperties.h"

class PotassiumLiquidFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  PotassiumLiquidFluidPropertiesTest() : MooseObjectUnitTest("PotassiumApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("PotassiumLiquidFluidProperties");
    _fe_problem->addUserObject("PotassiumLiquidFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<PotassiumLiquidFluidProperties>("fp");
  }

  const PotassiumLiquidFluidProperties * _fp;
};
