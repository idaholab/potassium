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
    _fp = &_fe_problem->getUserObjectTempl<PotassiumVaporFluidProperties>("fp");
  }

  const PotassiumVaporFluidProperties * _fp;
};
