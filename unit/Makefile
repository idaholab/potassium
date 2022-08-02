###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
# HERD_TRUNK_DIR   - Location of the HERD repository
# FRAMEWORK_DIR    - Location of the MOOSE framework
#
###############################################################################
POTASSIUM_DIR        ?= $(shell dirname `pwd`)
MOOSE_SUBMODULE    := $(POTASSIUM_DIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname $(POTASSIUM_DIR))/moose
endif
# check that MOOSE is available
MOOSE_CONTENT      := $(shell ls $(MOOSE_DIR) 2> /dev/null)
ifeq ($(MOOSE_CONTENT),)
  $(error MOOSE framework does not seem to be available. Make sure that either the submodule is checked out or that your MOOSE_DIR points to the correct location)
endif

MOOSE_DIR          ?= $(POTASSIUM_DIR)/moose
FRAMEWORK_DIR      ?= $(MOOSE_DIR)/framework
###############################################################################
CURRENT_DIR        := $(shell pwd)

# Extra stuff for GTEST
ADDITIONAL_INCLUDES += -I$(FRAMEWORK_DIR)/contrib/gtest
ADDITIONAL_LIBS     += $(FRAMEWORK_DIR)/contrib/gtest/libgtest.la

# framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk
################################## MODULES ####################################
FLUID_PROPERTIES  := yes
include           $(MOOSE_DIR)/modules/modules.mk
###############################################################################

include $(POTASSIUM_DIR)/libPotassiumProperties.mk

# dep apps
APPLICATION_DIR    := $(POTASSIUM_DIR)
APPLICATION_NAME   := potassium
include            $(FRAMEWORK_DIR)/app.mk

APPLICATION_DIR    := $(CURRENT_DIR)
APPLICATION_NAME   := potassium-unit
BUILD_EXEC         := yes
app_BASE_DIR       :=      # Intentionally blank
DEP_APPS           ?= $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include $(FRAMEWORK_DIR)/app.mk

# Find all the potassium unit test source files and include their dependencies.
potassium-unit_srcfiles := $(shell find $(CURRENT_DIR) -name "*.C")
potassium-unit_deps := $(patsubst %.C, %.$(obj-suffix).d, $(potassium-unit_srcfiles))
-include $(potassium-unit_deps)

###############################################################################
# Additional special case targets should be added here
