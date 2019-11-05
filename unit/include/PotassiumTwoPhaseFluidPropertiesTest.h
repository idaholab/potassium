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
    _fp = &_fe_problem->getUserObjectTempl<PotassiumTwoPhaseFluidProperties>("fp");
  }

  const PotassiumTwoPhaseFluidProperties * _fp;
};
