SUBDIRS := gpumcml mcml

all: $(SUBDIRS) copy
$(SUBDIRS):
		$(MAKE) -C $@

copy: $(SUBDIRS)
	cp gpumcml/gpumcml build/pymcml_gpumcml_executable
	cp mcml/mcml build/pymcml_mcml_executable

.PHONY: all $(SUBDIRS)
