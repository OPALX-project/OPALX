"""Build-system regression tests using small header fixtures, without OPALX execution."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


MODULE = Path(__file__).resolve().parents[3] / "cmake/PreferFetchedKokkosHeaders.cmake"


@unittest.skipUnless(shutil.which("cmake") and shutil.which("c++"), "CMake and C++ are required")
class PreferFetchedKokkosHeadersTest(unittest.TestCase):
    def test_conflicting_headers_build_only_with_precedence_fix_and_export_cleanly(self):
        with tempfile.TemporaryDirectory(prefix="opalx header ordering ") as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            for kind in ("fetched", "old-uenv"):
                for header in ("core/Kokkos_Core.hpp", "core/desul/atomics/Config.hpp",
                               "desul/desul/atomics.hpp", "mdspan/mdspan.hpp"):
                    path = source / kind / header
                    path.parent.mkdir(parents=True, exist_ok=True)
                    if kind == "old-uenv":
                        path.write_text('#error "stale uenv header selected"\n')
                    elif header.endswith("Kokkos_Core.hpp"):
                        path.write_text('#include <desul/atomics/Config.hpp>\n'
                                        '#include <desul/atomics.hpp>\n#include <mdspan.hpp>\n')
                    else:
                        path.write_text("// Matching fetched dependency header.\n")
            # The old uenv exports one broad include root containing all packages.
            old_include = source / "old-uenv/include"
            old_include.mkdir()
            for subdir in ("core", "desul", "mdspan"):
                shutil.copytree(source / "old-uenv" / subdir, old_include, dirs_exist_ok=True)
            (source / "probe.cpp").write_text('#include <Kokkos_Core.hpp>\nint probe() { return 0; }\n')
            (source / "consumer.cpp").write_text('#include <Kokkos_Core.hpp>\nint probe();\n'
                                                 'int main() { return probe(); }\n')
            (source / "CMakeLists.txt").write_text('''cmake_minimum_required(VERSION 3.20)
project(HeaderOrder LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include("''' + str(MODULE) + '''")
add_library(mpi INTERFACE IMPORTED)
set_property(TARGET mpi PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/old-uenv/include")
add_library(kokkoscore INTERFACE)
add_library(Kokkos::kokkoscore ALIAS kokkoscore)
target_include_directories(kokkoscore INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/fetched/core")
target_include_directories(kokkoscore SYSTEM INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fetched/desul>"
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/fetched/mdspan>")
add_library(ippl INTERFACE)
target_link_libraries(ippl INTERFACE kokkoscore mpi)
add_library(opalx_probe STATIC probe.cpp)
# Reproduce the actual direct/transitive dependency order and broad MPI prefix.
target_link_libraries(opalx_probe PUBLIC "$<BUILD_INTERFACE:mpi;ippl>")
if(APPLY_FIX)
  opalx_prefer_fetched_kokkos_headers(opalx_probe)
endif()
add_executable(consumer consumer.cpp)
target_link_libraries(consumer PRIVATE opalx_probe)
install(TARGETS opalx_probe EXPORT ProbeTargets DESTINATION lib)
install(EXPORT ProbeTargets DESTINATION lib/cmake/probe)
''')
            # Remove host compiler overrides so this fixture needs no CUDA installation.
            environment = dict(os.environ)
            environment.pop("CC", None)
            environment.pop("CXX", None)
            for enabled in (False, True):
                build = root / ("fixed" if enabled else "unfixed")
                configure = subprocess.run(
                    ["cmake", "-S", str(source), "-B", str(build), "-G", "Unix Makefiles",
                     "-DAPPLY_FIX=" + ("ON" if enabled else "OFF")],
                    env=environment, capture_output=True, text=True,
                )
                self.assertEqual(configure.returncode, 0, configure.stdout + configure.stderr)
                compiled = subprocess.run(["cmake", "--build", str(build), "--parallel", "2"],
                                          capture_output=True, text=True)
                with self.subTest(fix=enabled):
                    if not enabled:
                        self.assertNotEqual(compiled.returncode, 0)
                        self.assertIn("stale uenv header selected", compiled.stderr)
                        continue
                    self.assertEqual(compiled.returncode, 0, compiled.stdout + compiled.stderr)
                    commands = json.loads((build / "compile_commands.json").read_text())
                    self.assertEqual(len(commands), 2)
                    for entry in commands:
                        command = entry["command"]
                        for package in ("desul", "mdspan"):
                            self.assertLess(command.index("fetched/" + package),
                                            command.index("old-uenv/include"))
                    installed = subprocess.run(
                        ["cmake", "--install", str(build), "--prefix", str(root / "install")],
                        capture_output=True, text=True,
                    )
                    self.assertEqual(installed.returncode, 0, installed.stdout + installed.stderr)
                    exported = (root / "install/lib/cmake/probe/ProbeTargets.cmake").read_text()
                    self.assertNotIn(str(source), exported)
                    self.assertNotIn("fetched", exported)
                    self.assertNotIn("Kokkos::", exported)

    def test_missing_imported_or_non_system_kokkos_is_left_unchanged(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for kind in ("missing", "imported", "no-system-includes"):
                source = root / kind
                source.mkdir()
                if kind == "missing":
                    kokkos = ""
                elif kind == "imported":
                    kokkos = '''add_library(Kokkos::kokkoscore INTERFACE IMPORTED)
set_property(TARGET Kokkos::kokkoscore PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "/installed/kokkos")'''
                else:
                    kokkos = '''add_library(kokkoscore INTERFACE)
add_library(Kokkos::kokkoscore ALIAS kokkoscore)'''
                (source / "CMakeLists.txt").write_text('''cmake_minimum_required(VERSION 3.20)
project(Guard NONE)
include("''' + str(MODULE) + '''")
add_library(probe INTERFACE)
''' + kokkos + '''
opalx_prefer_fetched_kokkos_headers(probe)
get_target_property(includes probe INTERFACE_INCLUDE_DIRECTORIES)
if(includes)
  message(FATAL_ERROR "Unrelated include requirements changed: ${includes}")
endif()
''')
                result = subprocess.run(["cmake", "-S", str(source), "-B", str(source / "build")],
                                        capture_output=True, text=True)
                with self.subTest(kokkos=kind):
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
