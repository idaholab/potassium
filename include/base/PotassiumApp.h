#pragma once

#include "MooseApp.h"

class Factory;
class PotassiumApp;

template <>
InputParameters validParams<PotassiumApp>();

class PotassiumApp : public MooseApp
{
public:
  PotassiumApp(InputParameters parameters);

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
  static void registerObjects(Factory & factory);

protected:
};
