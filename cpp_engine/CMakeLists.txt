cmake_minimum_required(VERSION 3.20)
project(CHEMSI LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ============================================================
# Build options
# ============================================================
# The visualizer pulls GUI dependencies (OpenGL/GLFW/ImGui/ImPlot). Make it
# optional so a "core simulation" build works out-of-the-box.
option(CHEMSI_BUILD_VIS "Build CHEMSI visualizer (GLFW/OpenGL/ImGui)" OFF)
option(CHEMSI_ENABLE_GRPC "Enable Unity integration via gRPC + Protobuf" ON)
option(CHEMSI_USE_SYSTEM_GRPC "Prefer system-installed gRPC/Protobuf instead of FetchContent" ON)

# ============================================================
# Core library
# ============================================================
set(CHEMSI_CORE_SRCS
  src/chemistry.cpp
  src/LiIonRunaway.cpp
  src/ObjectModel.cpp
  src/MechanicsSim.cpp
  src/Reactor.cpp
  src/Simulation.cpp
  src/Aerodynamics.cpp
  src/Suppression.cpp
  src/Ventilation.cpp
)

# Only compile the gRPC server implementation when enabled.
if (CHEMSI_ENABLE_GRPC)
  list(APPEND CHEMSI_CORE_SRCS src/GrpcSimServer.cpp)
endif()

add_library(chemsi ${CHEMSI_CORE_SRCS})
target_include_directories(chemsi PUBLIC ${CMAKE_SOURCE_DIR}/include)

# ============================================================
# gRPC + Protobuf (Unity integration)
# ============================================================
if (CHEMSI_ENABLE_GRPC)

  # Proto layout (repo root has /proto next to /cpp_engine)
  set(CHEMSI_PROTO_DIR "${CMAKE_SOURCE_DIR}/../proto")
  set(CHEMSI_PROTO_GEN_DIR "${CMAKE_BINARY_DIR}/generated")

  set(CHEMSI_PROTO_FILES
    "${CHEMSI_PROTO_DIR}/command_v1.proto"
    "${CHEMSI_PROTO_DIR}/telemetry_frame_v1.proto"
    "${CHEMSI_PROTO_DIR}/vfep_types_v1.proto"
    "${CHEMSI_PROTO_DIR}/world_snapshot_v1.proto"
    "${CHEMSI_PROTO_DIR}/vfep_sim_service_v1.proto"
  )

  foreach(_p IN LISTS CHEMSI_PROTO_FILES)
    if (NOT EXISTS "${_p}")
      message(FATAL_ERROR "Missing proto file: ${_p}")
    endif()
  endforeach()

  # Map proto -> generated C++ sources/headers
  set(CHEMSI_PROTO_SRCS "")
  set(CHEMSI_PROTO_HDRS "")
  foreach(_p IN LISTS CHEMSI_PROTO_FILES)
    get_filename_component(_name "${_p}" NAME_WE)
    list(APPEND CHEMSI_PROTO_SRCS "${CHEMSI_PROTO_GEN_DIR}/${_name}.pb.cc")
    list(APPEND CHEMSI_PROTO_HDRS "${CHEMSI_PROTO_GEN_DIR}/${_name}.pb.h")
  endforeach()

  # gRPC stubs for the service proto
  set(CHEMSI_GRPC_SRCS "${CHEMSI_PROTO_GEN_DIR}/vfep_sim_service_v1.grpc.pb.cc")
  set(CHEMSI_GRPC_HDRS "${CHEMSI_PROTO_GEN_DIR}/vfep_sim_service_v1.grpc.pb.h")

  # ============================================================
  # Dependency selection
  # ============================================================
  if (CHEMSI_USE_SYSTEM_GRPC)
    # Prefer MSYS2 /mingw64 packaged deps
    if (MINGW OR MSYS)
      list(PREPEND CMAKE_PREFIX_PATH "/mingw64")
    endif()

    find_package(Protobuf CONFIG REQUIRED)
    find_package(gRPC CONFIG REQUIRED)

    # Resolve protoc path (MSYS2 packages export targets; we also support PATH fallback)
    set(CHEMSI_PROTOC_EXE "")
    if (TARGET protobuf::protoc)
      get_target_property(CHEMSI_PROTOC_EXE protobuf::protoc IMPORTED_LOCATION)
    endif()
    if (NOT CHEMSI_PROTOC_EXE)
      find_program(CHEMSI_PROTOC_EXE protoc REQUIRED)
    endif()

    # Resolve grpc_cpp_plugin path
    set(CHEMSI_GRPC_CPP_PLUGIN_EXE "")
    if (TARGET gRPC::grpc_cpp_plugin)
      get_target_property(CHEMSI_GRPC_CPP_PLUGIN_EXE gRPC::grpc_cpp_plugin IMPORTED_LOCATION)
    endif()
    if (NOT CHEMSI_GRPC_CPP_PLUGIN_EXE)
      find_program(CHEMSI_GRPC_CPP_PLUGIN_EXE grpc_cpp_plugin REQUIRED)
    endif()

    message(STATUS "Using system Protobuf/gRPC")
    message(STATUS "  protoc: ${CHEMSI_PROTOC_EXE}")
    message(STATUS "  grpc_cpp_plugin: ${CHEMSI_GRPC_CPP_PLUGIN_EXE}")

  else()
    # ----------------------------
    # FetchContent path (downloads/builds gRPC + Protobuf)
    # ----------------------------
    include(FetchContent)

    set(gRPC_INSTALL OFF CACHE BOOL "" FORCE)
    set(gRPC_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(protobuf_BUILD_TESTS OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
      grpc
      GIT_REPOSITORY https://github.com/grpc/grpc.git
      GIT_TAG        v1.62.2
    )
    FetchContent_MakeAvailable(grpc)

    if (NOT TARGET protobuf::protoc)
      message(FATAL_ERROR "FetchContent gRPC did not provide protobuf::protoc target.")
    endif()
    if (NOT TARGET gRPC::grpc_cpp_plugin)
      message(FATAL_ERROR "FetchContent gRPC did not provide gRPC::grpc_cpp_plugin target.")
    endif()

    set(CHEMSI_PROTOC_EXE "$<TARGET_FILE:protobuf::protoc>")
    set(CHEMSI_GRPC_CPP_PLUGIN_EXE "$<TARGET_FILE:gRPC::grpc_cpp_plugin>")

    message(STATUS "Using FetchContent Protobuf/gRPC")
  endif()

  # ============================================================
  # Proto generation
  # ============================================================
  file(MAKE_DIRECTORY "${CHEMSI_PROTO_GEN_DIR}")

  add_custom_command(
    OUTPUT
      ${CHEMSI_PROTO_SRCS} ${CHEMSI_PROTO_HDRS}
      ${CHEMSI_GRPC_SRCS}  ${CHEMSI_GRPC_HDRS}
    COMMAND "${CHEMSI_PROTOC_EXE}"
    ARGS
      -I "${CHEMSI_PROTO_DIR}"
      --cpp_out=generated
      --grpc_out=generated
      "--plugin=protoc-gen-grpc=${CHEMSI_GRPC_CPP_PLUGIN_EXE}"
      ${CHEMSI_PROTO_FILES}
    DEPENDS
      ${CHEMSI_PROTO_FILES}
    COMMENT "Generating protobuf/gRPC sources"
    VERBATIM
  )

  add_custom_target(chemsi_proto_gen
    DEPENDS
      ${CHEMSI_PROTO_SRCS} ${CHEMSI_PROTO_HDRS}
      ${CHEMSI_GRPC_SRCS}  ${CHEMSI_GRPC_HDRS}
  )

  add_library(chemsi_grpc_proto
    ${CHEMSI_PROTO_SRCS}
    ${CHEMSI_GRPC_SRCS}
  )
  add_dependencies(chemsi_grpc_proto chemsi_proto_gen)

  target_include_directories(chemsi_grpc_proto PUBLIC "${CHEMSI_PROTO_GEN_DIR}")
  target_link_libraries(chemsi_grpc_proto
    PUBLIC
      protobuf::libprotobuf
      gRPC::grpc++
  )

  target_compile_definitions(chemsi PUBLIC CHEMSI_ENABLE_GRPC=1)
  target_include_directories(chemsi PUBLIC "${CHEMSI_PROTO_GEN_DIR}")
  target_link_libraries(chemsi PUBLIC chemsi_grpc_proto)

endif()

# ============================================================
# Main simulation (console)
# ============================================================
add_executable(VFEP_Sim src/main.cpp)
target_link_libraries(VFEP_Sim PRIVATE chemsi)

add_executable(VFEP_GrpcClient src/grpc_client.cpp)
target_link_libraries(VFEP_GrpcClient PRIVATE chemsi)

# ============================================================
# Tests
# ============================================================
enable_testing()

add_executable(TestAtom tests/TestAtom.cpp)
target_link_libraries(TestAtom PRIVATE chemsi)
add_test(NAME TestAtom COMMAND TestAtom)

add_executable(NumericIntegrity tests/TestNumericIntegrity.cpp)
target_link_libraries(NumericIntegrity PRIVATE chemsi)
add_test(NAME NumericIntegrity COMMAND NumericIntegrity)

# MSVC Debug stack overflow fix for NumericIntegrity
if (MSVC)
  target_link_options(NumericIntegrity PRIVATE "$<$<CONFIG:Debug>:/STACK:8388608>")
endif()

# ============================================================
# Visualization executable: GLFW + OpenGL + ImGui (+ ImPlot)
# ============================================================
if (CHEMSI_BUILD_VIS)
  include(FetchContent)

  find_package(OpenGL REQUIRED)

  set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
  set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)

  FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG        3.4
  )
  FetchContent_MakeAvailable(glfw)

  # --- Dear ImGui ---
  FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        v1.91.1
  )
  FetchContent_MakeAvailable(imgui)

  # --- ImPlot ---
  FetchContent_Declare(
    implot
    GIT_REPOSITORY https://github.com/epezent/implot.git
    GIT_TAG        v0.16
  )
  FetchContent_MakeAvailable(implot)

  add_library(imgui_lib STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp

    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp

    ${implot_SOURCE_DIR}/implot.cpp
    ${implot_SOURCE_DIR}/implot_items.cpp
  )

  target_include_directories(imgui_lib PUBLIC
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
    ${implot_SOURCE_DIR}
  )

  target_link_libraries(imgui_lib PUBLIC
    glfw
    OpenGL::GL
  )

  add_executable(VFEP_Vis vis/main_vis.cpp)
  set_target_properties(VFEP_Vis PROPERTIES OUTPUT_NAME "VFEP")

  target_link_libraries(VFEP_Vis PRIVATE chemsi imgui_lib)
  target_include_directories(VFEP_Vis PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/vis
  )
endif()
