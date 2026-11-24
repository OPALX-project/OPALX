# -----------------------------------------------------------------------------
# CheckRunSideBoundary.cmake
#
# Checks the boundary between the input layer (parser, input objects, OPTION globals) and the
# run side (everything else under src/ except Main.cpp). A run-side source file must not
#   (a) use Options:: or OpalData outside a comment or string literal, or
#   (b) include an input-layer header; only the build bridge below may.
# With OPALX_LIBRARY_SOURCES, the file the build writes with the sources of the libraries opalx_run
# and opalx, it also checks that
#   (c) opalx_run compiles exactly the run-side .cpp files outside the build bridge, and opalx
#       exactly the input-layer and bridge ones,
# so that the library split and these path rules cannot drift apart.
# Prints every violation as <file>:<line> (or <file> for (c)) and fails if there is one.
#
# Usage: cmake -DOPALX_SOURCE_DIR=<repository>/src
#              [-DOPALX_LIBRARY_SOURCES=<build>/unit_tests/OPALXLibrarySources.cmake]
#              -P CheckRunSideBoundary.cmake
# -----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.20)

if(NOT IS_DIRECTORY "${OPALX_SOURCE_DIR}")
  message(FATAL_ERROR "OPALX_SOURCE_DIR='${OPALX_SOURCE_DIR}' is not a directory.")
endif()

# Rule (b) exceptions: the build bridge turns a run description into runtime objects through
# value-only static functions of the command classes.
set(BUILD_BRIDGE
    # build() and registerModel() of every element and time-dependence keyword class.
    "Track/BuildElements.cpp"
    # EmissionSource and FieldSolverCmd value readers for the samplers and the solver banner.
    "Track/BuildTracking.cpp"
    # BinningCmd, EmissionSource and FieldSolverCmd value readers for the space-charge setup.
    "SpaceCharge/SpaceChargeConfigBuilder.cpp")

# Input layer, as paths relative to src/.
function(is_input_layer path result)
  set(_input FALSE)
  if(path MATCHES "^(Attributes|Expressions|OpalParser|BasicActions|ValueDefinitions|OpalConfigure|Lines)/")
    set(_input TRUE)
  elseif(path MATCHES "^AbstractObjects/" AND NOT path MATCHES "^AbstractObjects/OpalParticle\\.")
    # OpalParticle is a runtime record of one particle.
    set(_input TRUE)
  elseif(path MATCHES "^Elements/" AND NOT path MATCHES "^Elements/(OpalBeamline|PlacementResolver)\\.")
    set(_input TRUE)
  elseif(path MATCHES "^Structure/(Beam|BinningCmd|FieldSolverCmd|EmissionSource)\\.[a-z]+$"
         OR path MATCHES "^Distribution/Distribution\\.[a-z]+$"
         OR path MATCHES "^Track/(TrackCmd|TrackRun|TrackEnd|FreezeRun)\\.[a-z]+$"
         OR path MATCHES "^Utilities/(Options|ParseError)\\.[a-z]+$")
    set(_input TRUE)
  endif()
  set(${result} ${_input} PARENT_SCOPE)
endfunction()

# CMake lists split at ';', keep ';' inside '[...]' together and treat '\;' as no separator, so
# these characters are replaced by control characters before the source is split into tokens.
string(ASCII 1 SEMICOLON)
string(ASCII 2 OPEN_BRACKET)
string(ASCII 3 CLOSE_BRACKET)
string(ASCII 4 BACKSLASH)
set(STRING_LITERAL "\"([^\"${BACKSLASH}\n]|${BACKSLASH}.)*\"")
set(CHAR_LITERAL "'([^'${BACKSLASH}\n]|${BACKSLASH}.)*'")
set(LINE_COMMENT "//[^\n]*")
set(BLOCK_COMMENT "/[*]([^*]|[*]+[^*/])*[*]+/")
# Tokens in order of preference; every other character sequence is code.
set(TOKEN "${STRING_LITERAL}|${CHAR_LITERAL}|${LINE_COMMENT}|${BLOCK_COMMENT}|[^\"'/]+|/")
set(GLOBAL_USE "(^|[^A-Za-z0-9_])(Options::|OpalData[^A-Za-z0-9_]|OpalData$)")

file(GLOB_RECURSE SOURCES RELATIVE "${OPALX_SOURCE_DIR}" "${OPALX_SOURCE_DIR}/*.cpp"
     "${OPALX_SOURCE_DIR}/*.h" "${OPALX_SOURCE_DIR}/*.hpp" "${OPALX_SOURCE_DIR}/*.tpp")
list(SORT SOURCES)

set(VIOLATIONS 0)
set(CHECKED 0)
foreach(path IN LISTS SOURCES)
  is_input_layer("${path}" input)
  if(input OR path STREQUAL "Main.cpp")
    continue()
  endif()
  math(EXPR CHECKED "${CHECKED} + 1")
  list(FIND BUILD_BRIDGE "${path}" bridge)
  get_filename_component(directory "${OPALX_SOURCE_DIR}/${path}" DIRECTORY)

  file(READ "${OPALX_SOURCE_DIR}/${path}" content)
  string(REPLACE "\\" "${BACKSLASH}" content "${content}")
  string(REPLACE ";" "${SEMICOLON}" content "${content}")
  string(REPLACE "[" "${OPEN_BRACKET}" content "${content}")
  string(REPLACE "]" "${CLOSE_BRACKET}" content "${content}")
  string(REGEX MATCHALL "${TOKEN}" tokens "${content}")

  set(line 1)
  set(previous "")
  foreach(token IN LISTS tokens)
    string(SUBSTRING "${token}" 0 2 head)
    if(head MATCHES "^\"")
      # The path of an #include "..." directive.
      if(previous MATCHES "#[ \t]*include[ \t]*$" AND bridge EQUAL -1)
        string(REGEX REPLACE "^\"(.*)\"$" "\\1" included "${token}")
        set(target "")
        if(EXISTS "${OPALX_SOURCE_DIR}/${included}")
          set(target "${OPALX_SOURCE_DIR}/${included}")
        elseif(EXISTS "${directory}/${included}")
          get_filename_component(target "${directory}/${included}" ABSOLUTE)
        endif()
        if(target)
          file(RELATIVE_PATH target "${OPALX_SOURCE_DIR}" "${target}")
          is_input_layer("${target}" input)
          if(input)
            message("src/${path}:${line}: includes the input-layer header \"${included}\"")
            math(EXPR VIOLATIONS "${VIOLATIONS} + 1")
          endif()
        endif()
      endif()
    elseif(NOT head MATCHES "^('|//|/[*])")
      set(previous "${token}")
      if(token MATCHES "${GLOBAL_USE}")
        # Find the lines of the code token that use a global.
        string(REPLACE "\n" ";" code_lines "${token}")
        set(offset 0)
        foreach(code_line IN LISTS code_lines)
          if(code_line MATCHES "${GLOBAL_USE}")
            set(global "OpalData")
            if(CMAKE_MATCH_2 STREQUAL "Options::")
              set(global "Options::")
            endif()
            math(EXPR at "${line} + ${offset}")
            string(STRIP "${code_line}" shown)
            string(REPLACE "${BACKSLASH}" "\\" shown "${shown}")
            string(REPLACE "${SEMICOLON}" ";" shown "${shown}")
            string(REPLACE "${OPEN_BRACKET}" "[" shown "${shown}")
            string(REPLACE "${CLOSE_BRACKET}" "]" shown "${shown}")
            message("src/${path}:${at}: uses ${global} outside a comment: ${shown}")
            math(EXPR VIOLATIONS "${VIOLATIONS} + 1")
          endif()
          math(EXPR offset "${offset} + 1")
        endforeach()
      endif()
    endif()
    string(REPLACE "\n" "" without_newlines "${token}")
    string(LENGTH "${token}" length)
    string(LENGTH "${without_newlines}" length_without_newlines)
    math(EXPR line "${line} + ${length} - ${length_without_newlines}")
  endforeach()
endforeach()

# Rule (c): the library of every .cpp file under src/ agrees with the rules above.
if(DEFINED OPALX_LIBRARY_SOURCES)
  include("${OPALX_LIBRARY_SOURCES}")
  if(NOT DEFINED opalx_run_SOURCES OR NOT DEFINED opalx_SOURCES)
    message(FATAL_ERROR "OPALX_LIBRARY_SOURCES='${OPALX_LIBRARY_SOURCES}' does not set "
                        "opalx_run_SOURCES and opalx_SOURCES.")
  endif()
  foreach(path IN LISTS SOURCES)
    if(NOT path MATCHES "\\.cpp$" OR path STREQUAL "Main.cpp")
      continue()
    endif()
    is_input_layer("${path}" input)
    list(FIND BUILD_BRIDGE "${path}" bridge)
    if(input OR NOT bridge EQUAL -1)
      set(expected "opalx")
    else()
      set(expected "opalx_run")
    endif()
    set(libraries "")
    foreach(library IN ITEMS opalx_run opalx)
      list(FIND ${library}_SOURCES "${path}" found)
      if(NOT found EQUAL -1)
        list(APPEND libraries "${library}")
      endif()
    endforeach()
    if(NOT libraries STREQUAL expected)
      if(libraries STREQUAL "")
        set(libraries "no library")
      endif()
      string(REPLACE ";" " and " libraries "${libraries}")
      message("src/${path}: compiled into ${libraries}, but the rules above put it into ${expected}")
      math(EXPR VIOLATIONS "${VIOLATIONS} + 1")
    endif()
  endforeach()
endif()

if(VIOLATIONS GREATER 0)
  message(FATAL_ERROR "${VIOLATIONS} run-side boundary violation(s) in ${CHECKED} run-side files.")
endif()
if(DEFINED OPALX_LIBRARY_SOURCES)
  message(STATUS "Run-side boundary: ${CHECKED} run-side files, no violations; opalx_run and opalx "
                 "agree with the rules.")
else()
  message(STATUS "Run-side boundary: ${CHECKED} run-side files, no violations.")
endif()
