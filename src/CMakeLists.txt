
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
set(orbit
    "${PROJECT_SOURCE_DIR}/src/Orbit/GuidingCenterEquation.cpp"
    "${PROJECT_SOURCE_DIR}/src/Orbit/Orbit.cpp"
    "${PROJECT_SOURCE_DIR}/src/Orbit/ParticleEquation.cpp"
    "${PROJECT_SOURCE_DIR}/src/Orbit/ParticlePusher.cpp"
    "${PROJECT_SOURCE_DIR}/src/Orbit/SOFTEquation.cpp"
)
set(phasespace
	"${PROJECT_SOURCE_DIR}/src/PhaseSpace/Particle.cpp"
	"${PROJECT_SOURCE_DIR}/src/PhaseSpace/Particle.momentum.cpp"
	"${PROJECT_SOURCE_DIR}/src/PhaseSpace/ParticleGenerator.cpp"
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
    "${PROJECT_SOURCE_DIR}/src/Tools/OutputModule.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Tool.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/ToolHandler.cpp"
)
set(tool_integrator
    "${PROJECT_SOURCE_DIR}/src/Tools/Integrator/Integrator.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Integrator/Welcome.cpp"
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
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/TimeIntegral.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/ToroidalIntegral/LocateSurfaceOfVisibility.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/ToroidalIntegral/sov_helpers.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/ToroidalIntegral/Trapz.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/ToroidalIntegral/TrapzImproved.cpp"
    "${PROJECT_SOURCE_DIR}/src/Tools/Radiation/Welcome.cpp"
)
set(mpi
    "${PROJECT_SOURCE_DIR}/src/SMPI.cpp"
)

set(all_tools ${tool_integrator} ${tool_orbits} ${tool_radiation})
set(core ${init} ${orbit} ${phasespace} ${soft} ${tools} ${all_tools})
set(target ${main} ${core})

if (WITH_MPI)
    set(target ${target} ${mpi})
endif (WITH_MPI)

add_executable(soft ${target})
target_include_directories(soft PUBLIC "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/include")

# Require C++17
set_target_properties(soft PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# Set LIBSTDC++ link flags
set(CMAKE_LIBSTDCXX_LINK_FLAGS "-lstdc++")
get_target_property(STDCXX_LINK_FLAGS soft LINK_FLAGS)
if (STDCXX_LINK_FLAGS STREQUAL "STDCXX_LINK_FLAGS-NOTFOUND")
    set_target_properties(soft PROPERTIES LINK_FLAGS "${CMAKE_LIBSTDCXX_LINK_FLAGS}")
else ()
    set_target_properties(soft PROPERTIES LINK_FLAGS "${CMAKE_LIBSTDCXX_LINK_FLAGS} ${STDCXX_LINK_FLAGS}")
endif ()

# Get git commit hash
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

####################
# MODULES
####################
# Find OpenMP
find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
    target_link_libraries(soft PUBLIC OpenMP::OpenMP_CXX)
endif ()

# Find SOFTLib
find_package(SOFTLIB REQUIRED)
if (SOFTLIB_FOUND)
    target_include_directories(soft PUBLIC ${SOFTLIB_INCLUDE_DIRS})
    target_link_libraries(soft PUBLIC ${SOFTLIB_LIBRARIES})
endif (SOFTLIB_FOUND)

# Find GSL!
find_package(GSL REQUIRED)
if (GSL_FOUND)
	if ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
		target_include_directories(soft PUBLIC ${GSL_INCLUDE_DIRS})
		message(STATUS ${GSL_LIBRARIES})
		target_link_libraries(soft PUBLIC ${GSL_LIBRARIES})
		target_link_libraries(soft PUBLIC ${GSL_CBLAS_LIBRARY})
	else ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
		message(FATAL_ERROR "SOFT needs GSL >= 2.0")
	endif ((GSL_VERSION VERSION_GREATER "2.0") OR (GSL_VERSION VERSION_EQUAL "2.0"))
endif (GSL_FOUND)

# Find libpng (optional dependency of 'softlib')
find_package(PNG)
if (PNG_FOUND)
    target_include_directories(soft PUBLIC ${PNG_INCLUDE_DIRS})
    target_link_libraries(soft PUBLIC ${PNG_LIBRARIES})
    target_compile_definitions(soft PUBLIC HAS_LIBPNG)
else (PNG_FOUND)
    message(INFO "No libpng was found. Configuring without PNG image support...")
endif (PNG_FOUND)

# Compile with official Matlab OR HDF5 support?
if (OFFICIAL_MATLAB)
	find_package(Matlab COMPONENTS MAT_LIBRARY MX_LIBRARY)
	if (Matlab_FOUND)
		target_include_directories(soft PUBLIC ${Matlab_INCLUDE_DIRS})
		target_link_libraries(soft PUBLIC ${Matlab_MAT_LIBRARY} ${Matlab_MX_LIBRARY})
	else (Matlab_FOUND)
		message(FATAL_ERROR "Official MATLAB support was requested, but no MATLAB installation was found")
	endif (Matlab_FOUND)
else (OFFICIAL_MATLAB)
	find_package(HDF5 COMPONENTS CXX)
	if (HDF5_FOUND)
		target_include_directories(soft PUBLIC ${HDF5_INCLUDE_DIRS})
		target_link_libraries(soft PUBLIC ${HDF5_LIBRARIES})
	else (HDF5_FOUND)
		message(FATAL_ERROR "HDF5 support was requested, but no HDF5 library was found on this system")
	endif (HDF5_FOUND)
endif (OFFICIAL_MATLAB)

# Compile with MPI support?
if (WITH_MPI)
    find_package(MPI COMPONENTS CXX)
    if (MPI_CXX_FOUND)
        target_include_directories(soft PUBLIC ${MPI_CXX_INCLUDE_DIRS})

        set_target_properties(soft PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_OPTIONS} ${MPI_CXX_COMPILE_DEFINITIONS}")
        set_target_properties(soft PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
        target_link_libraries(soft PUBLIC "${MPI_CXX_LIBRARIES}")
    else (MPI_CXX_FOUND)
        message(FATAL_ERROR "MPI support was requested, but no MPI library was found.")
    endif (MPI_CXX_FOUND)
endif (WITH_MPI)

