include_directories("${PROJECT_SOURCE_DIR}/include")
#include_directories("${PROJECT_SOURCE_DIR}/tests")

set(CMAKE_MODULE_PATH "${CMAKE_BINARY_DIR}/Modules/" CACHE PATH "CMake module path" FORCE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_LIBSTDCXX_LINK_FLAGS "-lstdc++")

option(BUILD_TESTS "Build unit tests" OFF)
option(COLOR_TERMINAL "Allow colorful SOFT output" ON)
option(DEBUG "Compile with debug symbols and no optimizations" OFF)
option(OPTIMIZE_NATIVE "Enable native optimizations (non-debug mode only)" ON)
option(PROFILING "Enable gprof profiling hooks" OFF)
option(WITH_MPI "Enable MPI support" OFF)

set(init
    "${PROJECT_SOURCE_DIR}/src/Init/Init.cpp"
    "${PROJECT_SOURCE_DIR}/src/Init/Config.cpp"
    "${PROJECT_SOURCE_DIR}/src/Init/Distribution.cpp"
    "${PROJECT_SOURCE_DIR}/src/Init/Global.cpp"
    "${PROJECT_SOURCE_DIR}/src/Init/MagneticField.cpp"
    "${PROJECT_SOURCE_DIR}/src/Init/ParticlePusher.cpp"
    "${PROJECT_SOURCE_DIR}/src/Init/RadialProfile.cpp"
    "${PROJECT_SOURCE_DIR}/src/Init/Tools.cpp"
    "${PROJECT_SOURCE_DIR}/src/Init/Local/ParticlePusher.cpp"
)
set(main
	"${PROJECT_SOURCE_DIR}/src/Main.cpp"
)
set(main_tests
	"${PROJECT_SOURCE_DIR}/tests/deep/runtest.cpp"
	"${PROJECT_SOURCE_DIR}/tests/deep/UnitTest.cpp"
)
set(orbit
    "${PROJECT_SOURCE_DIR}/src/Orbit/GuidingCenterEquation.cpp"
    "${PROJECT_SOURCE_DIR}/src/Orbit/Orbit.cpp"
    "${PROJECT_SOURCE_DIR}/src/Orbit/ParticleEquation.cpp"
    "${PROJECT_SOURCE_DIR}/src/Orbit/ParticlePusher.cpp"
    "${PROJECT_SOURCE_DIR}/src/Orbit/SOFTEquation.cpp"
)
set(orbit_tests
    "${PROJECT_SOURCE_DIR}/tests/deep/Orbit/ParticlePusher/Test_ParticlePusher.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Orbit/ParticlePusher/GCMagneticAxis.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Orbit/ParticlePusher/GuidingCenter.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Orbit/ParticlePusher/orbits.cpp"
)
set(phasespace
	"${PROJECT_SOURCE_DIR}/src/PhaseSpace/Particle.cpp"
	"${PROJECT_SOURCE_DIR}/src/PhaseSpace/Particle.momentum.cpp"
	"${PROJECT_SOURCE_DIR}/src/PhaseSpace/ParticleGenerator.cpp"
)
set(phasespace_tests
	"${PROJECT_SOURCE_DIR}/tests/deep/PhaseSpace/Test_Particle.cpp"
	"${PROJECT_SOURCE_DIR}/tests/deep/PhaseSpace/Test_ParticleGenerator.cpp"
)
set(soft
    "${PROJECT_SOURCE_DIR}/src/MemoryManager.cpp"
    "${PROJECT_SOURCE_DIR}/src/SOFT/Info.cpp"
    "${PROJECT_SOURCE_DIR}/src/SOFT/Run.cpp"
    "${PROJECT_SOURCE_DIR}/src/SOFT/SOFT.cpp"
    "${PROJECT_SOURCE_DIR}/src/SOFT/SOFTLocal.cpp"
	"${PROJECT_SOURCE_DIR}/src/softexit.cpp"
)
set(tools
    "${PROJECT_SOURCE_DIR}/src/Tools/Tool.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/ToolHandler.cpp"
)
set(tool_orbits
    "${PROJECT_SOURCE_DIR}/src/Tools/Orbits/Orbits.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Orbits/Configure.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Orbits/Output.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Orbits/Welcome.cpp"
)
set(tool_radiation
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Configure.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Detector.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/ToroidalIntegral/LocateSurfaceOfVisibility.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/ToroidalIntegral/sov_helpers.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/ToroidalIntegral/Trapz.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/ToroidalIntegral/TrapzImproved.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/AngularDistribution/Configure.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/AngularDistribution/Handle.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/AngularDistribution/Emission/ADEmission.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/AngularDistribution/Emission/Synchrotron/AngularDistribution.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/AngularDistribution/Emission/Synchrotron/AngularDistributionDrifts.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/AngularDistribution/Emission/Synchrotron/AngularSpectralDistribution.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/AngularDistribution/Emission/Synchrotron/SynchrotronEmission.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/AngularDistribution/Quadrature2D/ADQuadrature2D.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/AngularDistribution/Quadrature2D/Eval2D.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/AngularDistribution/Quadrature2D/Simpson2D.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/AngularDistribution/Quadrature2D/Trapz2D.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/Cone/Configure.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/Cone/EdgeCondition.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/Cone/Handle.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/Cone/Emission/ConeEmission.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/Cone/Emission/Bremsstrahlung/ConeBremsstrahlungEmission.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/Cone/Emission/Bremsstrahlung/ConeBremsstrahlungScreenedEmission.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/Cone/Emission/Bremsstrahlung/dilog_func.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/Cone/Emission/Synchrotron/ConeSynchrotronEmission.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/synchrotron_func.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/Cone/Emission/Unit/ConeUnitEmission.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/Cone/Projection/Original.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/Cone/Projection/Reverse.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/Isotropic/Configure.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Models/Isotropic/Isotropic.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Optics/Korger.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/RadiationOutput.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Green/Configure.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Green/Green.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Green/Handle.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Green/Output.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Green/Welcome.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Image/Configure.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Image/Generate.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Image/Handle.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Image/Image.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Image/Welcome.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/SoVVolume/Configure.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/SoVVolume/Generate.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/SoVVolume/Handle.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/SoVVolume/SoVVolume.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/SoVVolume/Welcome.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Space3D/Configure.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Space3D/Generate.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Space3D/Handle.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Space3D/Space3D.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Space3D/Welcome.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Spectrum/Configure.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Spectrum/Generate.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Spectrum/Handle.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Spectrum/Spectrum.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Spectrum/Welcome.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Topview/Configure.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Topview/Generate.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Topview/Handle.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Topview/Topview.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Output/Topview/Welcome.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Radiation.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/RadiationParticle.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Welcome.cpp"
)
set(tool_radiation_tests
    #"${PROJECT_SOURCE_DIR}/tests/deep/Tools/Radiation/Models/AngularDistribution/Test_AngularDistribution.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Tools/Radiation/Models/AngularDistribution/GOAMomentum.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Tools/Radiation/Models/AngularDistribution/Output_AngularDistributionDrifts.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Tools/Radiation/Models/AngularDistribution/Test_ADQuadrature2D.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Tools/Radiation/Models/AngularDistribution/Test_ADQuadrature2D_Emission.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Tools/Radiation/Models/AngularDistribution/Test_ADSynchrotronEmission.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Tools/Radiation/Models/AngularDistribution/Test_AngularDistributionDrifts.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Tools/Radiation/Models/Isotropic/Test_Isotropic.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Tools/Radiation/Models/Cone/Test_ConeProjection.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Tools/Radiation/Models/Cone/Test_SynchrotronEmission.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Tools/Radiation/Models/Cone/Bremsstrahlung/Test_BremsstrahlungScreenedEmission.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Tools/Radiation/Models/Cone/Bremsstrahlung/Test_BremsstrahlungEmission.cpp"
    "${PROJECT_SOURCE_DIR}/tests/deep/Tools/Radiation/Models/Cone/Reverse/Test_ConeProjectionReverse.cpp"
)

set(mpi
    "${PROJECT_SOURCE_DIR}/src/SMPI.cpp"
)

if (DEBUG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wno-format-security -O0 -g -D_FILE_OFFSET_BITS=64")
else (DEBUG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wno-format-security -O3 -D_FILE_OFFSET_BITS=64")
    if (OPTIMIZE_NATIVE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif (OPTIMIZE_NATIVE)
endif (DEBUG)

if (PROFILING)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg")
endif (PROFILING)

set(all_tools ${tool_orbits} ${tool_radiation})
set(core ${init} ${orbit} ${phasespace} ${soft} ${tools} ${all_tools})
set(target ${main} ${core})

if (WITH_MPI)
    set(target ${target} ${mpi})
endif (WITH_MPI)

add_executable(soft ${target})

# Set LIBSTDC++ link flags
get_target_property(STDCXX_LINK_FLAGS soft LINK_FLAGS)
if (STDCXX_LINK_FLAGS STREQUAL "STDCXX_LINK_FLAGS-NOTFOUND")
    set_target_properties(soft PROPERTIES LINK_FLAGS "${CMAKE_LIBSTDCXX_LINK_FLAGS}")
else ()
    set_target_properties(soft PROPERTIES LINK_FLAGS "${CMAKE_LIBSTDCXX_LINK_FLAGS} ${STDCXX_LINK_FLAGS}")
endif ()

if (BUILD_TESTS)
	set(core_tests ${main_tests} ${orbit_tests} ${phasespace_tests} ${soft} ${tool_radiation_tests})
	add_executable(soft_tests ${core} ${core_tests})
endif (BUILD_TESTS)

# Get git commit hash
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# Find OpenMP
find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
    target_link_libraries(soft PUBLIC OpenMP::OpenMP_CXX)

    if (BUILD_TESTS)
        target_link_libraries(soft_tests PUBLIC OpenMP::OpenMP_CXX)
    endif (BUILD_TESTS)
endif ()

# Find SOFTLib
find_package(SOFTLIB REQUIRED)
if (SOFTLIB_FOUND)
    include_directories(${SOFTLIB_INCLUDE_DIRS})
    target_link_libraries(soft PUBLIC ${SOFTLIB_LIBRARIES})
	
	if (BUILD_TESTS)
		target_link_libraries(soft_tests PUBLIC ${SOFTLIB_LIBRARIES})
	endif (BUILD_TESTS)
endif (SOFTLIB_FOUND)

# Find GSL!
find_package(GSL REQUIRED)
if (GSL_FOUND)
	if ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
		include_directories(${GSL_INCLUDE_DIRS})
		message(STATUS ${GSL_LIBRARIES})
		target_link_libraries(soft PUBLIC ${GSL_LIBRARIES})
		target_link_libraries(soft PUBLIC ${GSL_CBLAS_LIBRARY})

		if (BUILD_TESTS)
			target_link_libraries(soft_tests PUBLIC ${GSL_LIBRARIES})
			target_link_libraries(soft_tests PUBLIC ${GSL_CBLAS_LIBRARY})
		endif (BUILD_TESTS)
	else ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
		message(FATAL_ERROR "SOFT needs GSL >= 2.0")
	endif ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
endif (GSL_FOUND)

# Compile with official Matlab OR HDF5 support?
if (OFFICIAL_MATLAB)
	find_package(Matlab COMPONENTS MAT_LIBRARY MX_LIBRARY)
	if (Matlab_FOUND)
		include_directories(${Matlab_INCLUDE_DIRS})
		target_link_libraries(soft PUBLIC ${Matlab_MAT_LIBRARY} ${Matlab_MX_LIBRARY})

		if (BUILD_TESTS)
			target_link_libraries(soft_tests PUBLIC ${Matlab_MAT_LIBRARY} ${Matlab_MX_LIBRARY})
		endif (BUILD_TESTS)
	else (Matlab_FOUND)
		message(FATAL_ERROR "Official MATLAB support was requested, but no MATLAB installation was found")
	endif (Matlab_FOUND)
else (OFFICIAL_MATLAB)
	find_package(HDF5 COMPONENTS CXX)
	if (HDF5_FOUND)
		include_directories(${HDF5_INCLUDE_DIRS})
		target_link_libraries(soft PUBLIC ${HDF5_LIBRARIES})

		if (BUILD_TESTS)
			target_link_libraries(soft_tests PUBLIC ${HDF5_LIBRARIES})
		endif (BUILD_TESTS)
	else (HDF5_FOUND)
		message(FATAL_ERROR "HDF5 support was requested, but no HDF5 library was found on this system")
	endif (HDF5_FOUND)
endif (OFFICIAL_MATLAB)

# Compile with MPI support?
if (WITH_MPI)
    find_package(MPI COMPONENTS CXX)
    if (MPI_CXX_FOUND)
        include_directories(${MPI_CXX_INCLUDE_DIRS})

        set_target_properties(soft PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_OPTIONS} ${MPI_CXX_COMPILE_DEFINITIONS}")
        set_target_properties(soft PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
        target_link_libraries(soft PUBLIC "${MPI_CXX_LIBRARIES}")
    else (MPI_CXX_FOUND)
        message(FATAL_ERROR "MPI support was requested, but no MPI library was found.")
    endif (MPI_CXX_FOUND)
endif (WITH_MPI)

configure_file(
	"${PROJECT_SOURCE_DIR}/include/config.h.in"
	"${PROJECT_SOURCE_DIR}/include/config.h"
)

