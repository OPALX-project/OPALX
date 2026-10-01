"""Configure the actual OPALX options without compilers or external dependencies."""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


SOURCE = Path(__file__).resolve().parents[2]


@unittest.skipUnless(shutil.which("cmake"), "CMake is required")
class BuildConfigurationTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="opalx build options ")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "source"
        self.source.mkdir()
        (self.source / "CMakeLists.txt").write_text('''cmake_minimum_required(VERSION 3.24)
project(OPALX VERSION 0.9 LANGUAGES NONE)
list(PREPEND CMAKE_MODULE_PATH "''' + SOURCE.as_posix() + '''/cmake")
include(Messages)
include(OPALXOptions)
include(ProjectSetup)
file(WRITE "${CMAKE_BINARY_DIR}/build-type.txt" "${CMAKE_BUILD_TYPE}")
file(WRITE "${CMAKE_BINARY_DIR}/configurations.txt" "${CMAKE_CONFIGURATION_TYPES}")
''')
        self.environment = dict(os.environ)
        # Test default selection independently of the caller's build settings.
        for name in ("CMAKE_BUILD_TYPE", "CMAKE_GENERATOR", "CMAKE_GENERATOR_PLATFORM",
                     "CMAKE_GENERATOR_TOOLSET", "CMAKE_GENERATOR_INSTANCE"):
            self.environment.pop(name, None)
        self.generator = "Ninja" if shutil.which("ninja") else "Unix Makefiles"

    def configure(self, name, *options, generator=None):
        build = self.root / name
        result = subprocess.run(
            ["cmake", "-S", str(self.source), "-B", str(build),
             "-G", generator or self.generator, *options],
            env=self.environment, capture_output=True, text=True, timeout=20,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return (build / "build-type.txt").read_text()

    def test_initial_selection(self):
        cases = (
            ("default", "Release", ()),
            ("empty", "Release", ("-DCMAKE_BUILD_TYPE=",)),
            ("standard-debug", "Debug", ("-DCMAKE_BUILD_TYPE=Debug",)),
            ("legacy", "RelWithDebInfo", ("-DBUILD_TYPE=RelWithDebInfo",)),
            ("conflict", "Debug", ("-DCMAKE_BUILD_TYPE=Debug", "-DBUILD_TYPE=Release")),
        )
        for name, expected, options in cases:
            with self.subTest(case=name):
                self.assertEqual(self.configure(name, *options), expected)

    def test_reconfigure_preserves_standard_type_until_explicitly_changed(self):
        self.assertEqual(self.configure("reuse", "-DBUILD_TYPE=Debug"), "Debug")
        self.assertEqual(self.configure("reuse"), "Debug")
        self.assertEqual(self.configure("reuse", "-DBUILD_TYPE=Release"), "Debug")
        self.assertEqual(self.configure("reuse", "-DCMAKE_BUILD_TYPE=MinSizeRel"), "MinSizeRel")
        self.assertEqual(self.configure("reuse", "-DCMAKE_BUILD_TYPE="), "Release")

    def test_existing_release_cache_can_switch_to_debug(self):
        self.assertEqual(self.configure("dashboard"), "Release")
        self.assertEqual(self.configure("dashboard", "-DCMAKE_BUILD_TYPE=Debug"), "Debug")
        self.assertEqual(self.configure("dashboard"), "Debug")

    @unittest.skipUnless(shutil.which("ninja"), "Ninja is required for multi-config checks")
    def test_multi_config_does_not_force_single_build_type(self):
        for name, options, expected in (
            ("multi-default", (), ""),
            ("multi-legacy", ("-DBUILD_TYPE=Release",), ""),
            ("multi-standard", ("-DCMAKE_BUILD_TYPE=Debug", "-DBUILD_TYPE=Release"), "Debug"),
        ):
            with self.subTest(case=name):
                self.assertEqual(self.configure(name, *options, generator="Ninja Multi-Config"), expected)
                configurations = (self.root / name / "configurations.txt").read_text().split(";")
                self.assertIn("Debug", configurations)
                self.assertIn("Release", configurations)


if __name__ == "__main__":
    unittest.main()
