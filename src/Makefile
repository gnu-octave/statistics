# Makefile for compiling required oct files (parallel version)

# List all target oct files
OCTFILES = editDistance.oct libsvmread.oct libsvmwrite.oct svmpredict.oct \
           svmtrain.oct fcnntrain.oct fcnnpredict.oct \
           gamtrain.oct gampredict.oct __lbfgs__.oct

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

# fcnn.cpp is #included by both, so it is a prerequisite of both: without it
# an edit to the network code leaves the oct-files stale and the build silent.
fcnntrain.oct: fcnntrain.cc fcnn.cpp
	$(MKOCTFILE) fcnntrain.cc

fcnnpredict.oct: fcnnpredict.cc fcnn.cpp
	$(MKOCTFILE) fcnnpredict.cc

# gam.cpp is #included by both, exactly as fcnn.cpp is, and for the same
# reason: the engine is shared between the classification and the regression
# learner, and an edit to it must rebuild both oct-files.
gamtrain.oct: gamtrain.cc gam.cpp
	$(MKOCTFILE) gamtrain.cc

gampredict.oct: gampredict.cc gam.cpp
	$(MKOCTFILE) gampredict.cc

# lbfgs.h is the engine both this wrapper and any compiled learner include, so
# an edit to it must rebuild the oct-file.
__lbfgs__.oct: __lbfgs__.cc lbfgs.h
	$(MKOCTFILE) __lbfgs__.cc

# Clean target
clean:
	rm -f $(OCTFILES)

# Declare targets that don't create files
.PHONY: all clean

