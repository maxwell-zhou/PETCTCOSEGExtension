cmake_minimum_required(VERSION 2.8.9)

project(PETCTCOSEGExtension)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/maxwell-zhou/PETCTCOSEGExtension")
set(EXTENSION_CATEGORY "Segmentation")
set(EXTENSION_CONTRIBUTORS "Qi Song, Junjie Bai (University of Iowa)")
set(EXTENSION_DESCRIPTION "The PETCTCOSEG extension implementes a co-segmentation of tumor in PET and CT simultaneously.")
set(EXTENSION_ICONURL "https://raw.github.com/maxwell-zhou/PETCTCOSEGExtension/master/PETCTCOSEGExtension.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.github.com/llasion/PETCTCOSEGExtension/master/PETCTCOSEGExtension.png")
set(EXTENSION_DEPENDS "ITK") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(PETCTCOSEG)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
