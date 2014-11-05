#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=MinGW-Windows
CND_DLIB_EXT=dll
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/1804899502/gtest-all.o \
	${OBJECTDIR}/_ext/1804899502/gtest_main.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f1

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libgtest.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libgtest.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libgtest.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libgtest.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libgtest.a

${OBJECTDIR}/_ext/1804899502/gtest-all.o: ../gtest-1.7.0/src/gtest-all.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1804899502
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../gtest-1.7.0 -I../gtest-1.7.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1804899502/gtest-all.o ../gtest-1.7.0/src/gtest-all.cc

${OBJECTDIR}/_ext/1804899502/gtest_main.o: ../gtest-1.7.0/src/gtest_main.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1804899502
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../gtest-1.7.0 -I../gtest-1.7.0/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1804899502/gtest_main.o ../gtest-1.7.0/src/gtest_main.cc

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}
${TESTDIR}/TestFiles/f1: ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS} ../FemBase/dist/Debug/MinGW-Windows/libfembase.a ../dist/Debug/MinGW-Windows/libfem.a 


${OBJECTDIR}/_ext/1804899502/gtest-all_nomain.o: ${OBJECTDIR}/_ext/1804899502/gtest-all.o ../gtest-1.7.0/src/gtest-all.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1804899502
	@NMOUTPUT=`${NM} ${OBJECTDIR}/_ext/1804899502/gtest-all.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -I../gtest-1.7.0 -I../gtest-1.7.0/include -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1804899502/gtest-all_nomain.o ../gtest-1.7.0/src/gtest-all.cc;\
	else  \
	    ${CP} ${OBJECTDIR}/_ext/1804899502/gtest-all.o ${OBJECTDIR}/_ext/1804899502/gtest-all_nomain.o;\
	fi

${OBJECTDIR}/_ext/1804899502/gtest_main_nomain.o: ${OBJECTDIR}/_ext/1804899502/gtest_main.o ../gtest-1.7.0/src/gtest_main.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1804899502
	@NMOUTPUT=`${NM} ${OBJECTDIR}/_ext/1804899502/gtest_main.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -I../gtest-1.7.0 -I../gtest-1.7.0/include -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1804899502/gtest_main_nomain.o ../gtest-1.7.0/src/gtest_main.cc;\
	else  \
	    ${CP} ${OBJECTDIR}/_ext/1804899502/gtest_main.o ${OBJECTDIR}/_ext/1804899502/gtest_main_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/f1 || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/libgtest.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
