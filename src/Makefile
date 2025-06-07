# Makefile for compiling required oct files (parallel version)

# List all target oct files
OCTFILES = editDistance.oct libsvmread.oct libsvmwrite.oct svmpredict.oct \
           svmtrain.oct fcnntrain.oct fcnnpredict.oct

# Default target - depends on all oct files (enables parallel build)
all: $(OCTFILES)

# Individual targets for each oct file with dependencies
editDistance.oct: editDistance.cc
	$(MKOCTFILE) editDistance.cc

libsvmread.oct: libsvmread.cc
	$(MKOCTFILE) libsvmread.cc

libsvmwrite.oct: libsvmwrite.cc
	$(MKOCTFILE) libsvmwrite.cc

svmpredict.oct: svmpredict.cc svm.cpp svm_model_octave.cc
	$(MKOCTFILE) svmpredict.cc svm.cpp svm_model_octave.cc

svmtrain.oct: svmtrain.cc svm.cpp svm_model_octave.cc
	$(MKOCTFILE) svmtrain.cc svm.cpp svm_model_octave.cc

fcnntrain.oct: fcnntrain.cc
	$(MKOCTFILE) fcnntrain.cc

fcnnpredict.oct: fcnnpredict.cc
	$(MKOCTFILE) fcnnpredict.cc

# Clean target
clean:
	rm -f $(OCTFILES)

# Declare targets that don't create files
.PHONY: all clean

