#pragma once

#include "MooseApp.h"

class Factory;

class PotassiumApp : public MooseApp
{
public:
  PotassiumApp(InputParameters parameters);

public:
  static InputParameters validParams();
  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
  static void registerObjects(Factory & factory);
};
