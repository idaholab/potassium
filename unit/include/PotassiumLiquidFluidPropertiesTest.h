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
