#pragma once

#include "MooseObjectUnitTest.h"
#include "Potassium7EqnFluidProperties.h"

class Potassium7EqnFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  Potassium7EqnFluidPropertiesTest() : MooseObjectUnitTest("PotassiumApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("Potassium7EqnFluidProperties");
    _fe_problem->addUserObject("Potassium7EqnFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObjectTempl<Potassium7EqnFluidProperties>("fp");
  }

  const Potassium7EqnFluidProperties * _fp;
};
