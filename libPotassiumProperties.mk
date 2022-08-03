#
# build libPotassiumProperties using libMesh's build system
#
LIBPOTASSIUM_DIR       := $(POTASSIUM_DIR)/contrib/libPotassiumProperties

LIBPOTASSIUM_srcfiles  :=
LIBPOTASSIUM_srcfiles  += $(LIBPOTASSIUM_DIR)/K_Golden.cpp

LIBPOTASSIUM_objects   := $(patsubst %.cpp, %.$(obj-suffix), $(LIBPOTASSIUM_srcfiles))
LIBPOTASSIUM_deps      := $(patsubst %.$(obj-suffix), %.$(obj-suffix).d, $(LIBPOTASSIUM_objects))
LIBPOTASSIUM_LIB       := $(LIBPOTASSIUM_DIR)/libPotassiumProperties-$(METHOD).la

app_INCLUDES += -I$(POTASSIUM_DIR)
app_LIBS += $(LIBPOTASSIUM_LIB)

$(LIBPOTASSIUM_LIB): $(LIBPOTASSIUM_objects)
	@echo "Linking Library "$@"..."
	@$(libmesh_LIBTOOL) --tag=CC $(LIBTOOLFLAGS) --mode=link --quiet \
	  $(libmesh_CC) $(libmesh_CFLAGS) -o $@ $(LIBPOTASSIUM_objects) $(libmesh_LDFLAGS) $(EXTERNAL_FLAGS) -rpath $(LIBPOTASSIUM_DIR)
	@$(libmesh_LIBTOOL) --mode=install --quiet install -c $(LIBPOTASSIUM_LIB) $(LIBPOTASSIUM_DIR)

$(app_EXEC): $(LIBPOTASSIUM_LIB)

-include $(LIBPOTASSIUM_deps)

cleanlibpotassium:
	@echo "Cleaning libPotassiumProperties"
	@rm -f $(LIBPOTASSIUM_objects)
	@rm -f $(LIBPOTASSIUM_deps)
	@rm -f $(LIBPOTASSIUM_LIB)
	@rm -f $(LIBPOTASSIUM_DIR)/libPotassiumProperties-$(METHOD)*.dylib
	@rm -f $(LIBPOTASSIUM_DIR)/libPotassiumProperties-$(METHOD)*.so*
	@rm -f $(LIBPOTASSIUM_DIR)/libPotassiumProperties-$(METHOD)*.a
