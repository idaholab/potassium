//* This file is part of potassium
//* https://github.com/idaholab/potassium
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/potassium/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PotassiumApp.h"
#include "PotassiumRevision.h"
#include "MooseSyntax.h"
#include "AppFactory.h"

// Modules
#include "FluidPropertiesApp.h"

InputParameters
PotassiumApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_output_syntax") = false;
  return params;
}

registerKnownLabel("PotassiumApp");

PotassiumApp::PotassiumApp(InputParameters parameters)
  : MooseApp(parameters)
{
  PotassiumApp::registerAll(_factory, _action_factory, _syntax);
}

// External entry point for dynamic application loading
extern "C" void
PotassiumApp__registerApps()
{
  PotassiumApp::registerApps();
}

void
PotassiumApp::registerApps()
{
  registerApp(PotassiumApp);
}

extern "C" void
PotassiumApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  PotassiumApp::registerAll(f, af, s);
}

void
PotassiumApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  Registry::registerObjectsTo(f, {"PotassiumApp"});
  Registry::registerActionsTo(af, {"PotassiumApp"});

  FluidPropertiesApp::registerAll(f, af, s);
}

void
PotassiumApp::registerObjects(Factory & factory)
{
  mooseDeprecated("PotassiumApp: use registerAll instead of registerObjects");
  Registry::registerObjectsTo(factory, {"PotassiumApp"});
}
