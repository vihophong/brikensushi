PACKAGE = myanalysis

DECODER_PACKAGE = pdecoder

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --glibs)

LIBPATH = /home/phong/projects/rcdaqnew/install/bin/newanalaysis

CXXFLAGS = -I.  $(ROOTFLAGS) -I$(ONLINE_MAIN)/include -I$(OFFLINE_MAIN)/include
RCFLAGS = -I.  -I$(ONLINE_MAIN)/include -I$(OFFLINE_MAIN)/include -I/usr/lib
AFLAGS = -I. -I/usr/lib

LDFLAGS = -Wl,--no-as-needed  -L$(ONLINE_MAIN)/lib  -L$(OFFLINE_MAIN)/lib -lpmonitor -lEvent -lNoRootEvent -lmessage $(ROOTLIBS) -fPIC 

ALDFLAGS = -Wl,--no-as-needed $(ROOTLIBS) -fPIC


HDRFILES = $(DECODER_PACKAGE).h
LINKFILE = $(DECODER_PACKAGE)LinkDef.h

APKG = DataStruct
AHDRFILES = libDataStruct.h
ALINKFILE = libDataStructLinkDef.h

ADDITIONAL_SOURCES = libDataStruct.h dpp.cpp dpp.h
ADDITIONAL_LIBS = -L$(PWD) -$(APKG)

SO = lib$(PACKAGE).so

DATASO = lib$(APKG).so

$(SO) : $(DECODER_PACKAGE).cc $(PACKAGE).cc $(DECODER_PACKAGE)_dict.C $(ADDITIONAL_SOURCES) $(LINKFILE) $(APKG)_dict.C $(DATASO)
	$(CXX) $(CXXFLAGS) -o $@ -shared  $<  $(ADDITIONAL_SOURCES) $(DECODER_PACKAGE)_dict.C $(APKG)_dict.C $(LDFLAGS) $(ADDITIONAL_LIBS)


$(DATASO) : $(AHDRFILES) $(APKG)_dict.C
	$(CXX) -I. $(ROOTFLAGS) -o $@ -shared $^ $(ALDFLAGS)

$(DECODER_PACKAGE)_dict.C : $(HDRFILES) $(LINKFILE)
	rootcint -f $@  -c $(RCFLAGS) $^

$(APKG)_dict.C: $(AHDRFILES) $(ALINKFILE)
	rootcint -f $@  -c $(AFLAGS) $^

 
.PHONY: clean

clean: 
	rm -f $(SO) $(DATASO) $(DECODER_PACKAGE)_dict.C $(DECODER_PACKAGE)_dict.h $(APKG)_dict.C $(APKG)_dict.h 

