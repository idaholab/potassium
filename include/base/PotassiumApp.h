//* This file is part of potassium
//* https://github.com/idaholab/potassium
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/potassium/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseApp.h"

class Factory;

class PotassiumApp : public MooseApp
{
public:
  PotassiumApp(const InputParameters & parameters);

public:
  static InputParameters validParams();
  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
  static void registerObjects(Factory & factory);
};
