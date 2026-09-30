"""Unit and preprocessing checks; no OPALX build or execution is involved."""

import contextlib
import io
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from diagnose_desul_headers import main, preprocess_arguments, selected_headers


class DesulDiagnosticTest(unittest.TestCase):
    def test_preserves_launchers_and_include_precedence_without_object_outputs(self):
        command = shlex.join([
            "/kokkos/kokkos_launch_compiler", "/kokkos/nvcc_wrapper", "/env/g++",
            "/env/g++", "-isystem", "/old view/include", "-isystem", "/bundled/desul",
            "-DKEEP=1", "-std=c++20", "-arch=sm_90", "-MD", "-MMD", "-MP",
            "-MF", "object.d", "-MTobject.o", "-MQ", "object.o", "-o", "object.o",
            "-oattached.o", "-c", "/source/BeamlineCore/RFCavityRep.cpp",
        ])
        self.assertEqual(preprocess_arguments(command), [
            "/kokkos/kokkos_launch_compiler", "/kokkos/nvcc_wrapper", "/env/g++",
            "/env/g++", "-isystem", "/old view/include", "-isystem", "/bundled/desul",
            "-DKEEP=1", "-std=c++20", "-arch=sm_90",
            "/source/BeamlineCore/RFCavityRep.cpp", "-E",
        ])

    def test_rejects_non_compile_or_shell_commands(self):
        for command in ("ninja", "g++ -c source.cpp && touch marker", "g++ -c source.cpp -o"):
            with self.subTest(command=command), self.assertRaises(ValueError):
                preprocess_arguments(command)

    def test_header_report_uses_markers_not_wrapper_calls(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "probe.ii"
            output.write_text(
                '# 1 "/old/include/desul/atomics.hpp" 1\n'
                '#line 2 "/old/include/desul/atomics/Generic.hpp"\n'
                '# 3 "/kokkos/Kokkos_Atomics_Desul_Wrapper.hpp"\n'
                'desul::atomic_mod(ptr, val);\n'
                '# 8 "/old/include/desul/atomics.hpp" 2\n'
                '# 1 "/unrelated/header.hpp"\n'
            )
            self.assertEqual(selected_headers(output), [
                "/kokkos/Kokkos_Atomics_Desul_Wrapper.hpp",
                "/old/include/desul/atomics.hpp",
                "/old/include/desul/atomics/Generic.hpp",
            ])

    @unittest.skipUnless(shutil.which("c++"), "A C++ preprocessor is required")
    def test_real_preprocessor_identifies_shadowing_and_preserves_existing_outputs(self):
        with tempfile.TemporaryDirectory(prefix="desul diagnostic ") as temporary:
            root = Path(temporary)
            source = root / "BeamlineCore/RFCavityRep.cpp"
            source.parent.mkdir()
            source.write_text('#include <desul/atomics.hpp>\n')
            for name in ("old view", "bundled"):
                header = root / name / "desul/atomics.hpp"
                header.parent.mkdir(parents=True)
                header.write_text("// " + name + "\n")
            object_file = root / "object.o"
            dependency = root / "object.d"
            object_file.write_text("existing object")
            dependency.write_text("existing dependencies")
            marker = root / "launcher-used"
            launcher = root / "launcher"
            launcher.write_text(
                "#!" + sys.executable + "\nimport os,sys\n"
                "open(" + repr(str(marker)) + ", 'w').write('yes')\n"
                "os.execv(sys.argv[1], sys.argv[1:])\n"
            )
            launcher.chmod(0o755)
            ninja = root / "ninja"
            for first, second in (("old view", "bundled"), ("bundled", "old view")):
                command = shlex.join([
                    str(launcher), shutil.which("c++"), "-isystem", str(root / first),
                    "-isystem", str(root / second), "-MD", "-MF", str(dependency),
                    "-o", str(object_file), "-c", str(source),
                ])
                ninja.write_text("#!" + sys.executable + "\nprint(" + repr(command) + ")\n")
                ninja.chmod(0o755)
                report = io.StringIO()
                with patch.dict(os.environ, {"PATH": str(root) + os.pathsep + os.environ["PATH"]}):
                    with contextlib.redirect_stdout(report):
                        status = main([str(root)])
                with self.subTest(selected=first):
                    self.assertEqual(status, 0, report.getvalue())
                    headers = report.getvalue().split("Selected Desul/Kokkos wrapper headers:")[1]
                    self.assertIn(str(root / first / "desul/atomics.hpp"), headers)
                    self.assertNotIn(str(root / second / "desul/atomics.hpp"), headers)
                    self.assertEqual(marker.read_text(), "yes")
                    self.assertEqual(object_file.read_text(), "existing object")
                    self.assertEqual(dependency.read_text(), "existing dependencies")

    def test_timeout_stops_the_compiler_process_group(self):
        with tempfile.TemporaryDirectory() as temporary:
            command = "compiler -c /source/BeamlineCore/RFCavityRep.cpp -o object.o"
            process = Mock(pid=12345)
            process.wait.side_effect = [subprocess.TimeoutExpired(command, 90), -9]
            with patch("diagnose_desul_headers.subprocess.run", return_value=Mock(stdout=command)):
                with patch("diagnose_desul_headers.subprocess.Popen", return_value=process):
                    with patch("diagnose_desul_headers.os.killpg") as kill_group:
                        report = io.StringIO()
                        with contextlib.redirect_stdout(report):
                            self.assertEqual(main([temporary]), 1)
                        kill_group.assert_called_once_with(12345, signal.SIGKILL)
                        self.assertIn("Diagnostic incomplete:", report.getvalue())

    @unittest.skipUnless(shutil.which("cmake"), "CMake is required")
    def test_dashboard_publishes_notes_without_changing_build_status(self):
        dashboard = Path(__file__).resolve().parents[1] / "dashboard-configure-build.cmake"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cases = [("enabled", True, 0, 0), ("probe-fails", True, 0, 1),
                     ("disabled", False, 0, 0), ("configure-fails", True, 1, 0)]
            for name, enabled, configure_code, diagnostic_code in cases:
                case = root / name
                source, build = case / "source", case / "build"
                build.mkdir(parents=True)
                helper = source / "ci/cscs/diagnose_desul_headers.py"
                helper.parent.mkdir(parents=True)
                helper.write_text("import sys\nprint('probe report')\nsys.exit({})\n".format(diagnostic_code))
                calls = case / "calls.txt"
                harness = case / "harness.cmake"
                harness.write_text(f'''cmake_minimum_required(VERSION 3.20)
set(ENV{{CI_PROJECT_DIR}} "{source}")
set(BUILD_DIR "{build}")
set(OPALX_DIAG_DESUL_HEADERS {'ON' if enabled else 'OFF'})
set(DIAGNOSTIC_PYTHON "{sys.executable}")
macro(ctest_start)
endmacro()
function(ctest_configure)
  set(configure_result {configure_code} PARENT_SCOPE)
endfunction()
function(ctest_submit)
  file(APPEND "{calls}" "submit:${{ARGV}}\\n")
endfunction()
function(ctest_build)
  file(APPEND "{calls}" "build\\n")
  set(build_result 0 PARENT_SCOPE)
endfunction()
include("{dashboard}")
''')
                result = subprocess.run(["cmake", "-P", str(harness)], capture_output=True, text=True)
                report = build / "desul-header-diagnostic.log"
                expected_report = enabled and configure_code == 0
                with self.subTest(case=name):
                    self.assertEqual(report.exists(), expected_report, result.stderr)
                    self.assertEqual("submit:PARTS;Notes" in calls.read_text(), expected_report)
                    self.assertIn("build", calls.read_text())
                    self.assertEqual(result.returncode != 0, configure_code != 0, result.stderr)
                    if expected_report:
                        self.assertIn("probe report", report.read_text())
                        self.assertIn("Diagnostic exit status: " + str(diagnostic_code), report.read_text())

    def test_missing_ninja_is_reported_as_inconclusive(self):
        with tempfile.TemporaryDirectory() as temporary:
            report = io.StringIO()
            with patch.dict(os.environ, {"PATH": temporary}), contextlib.redirect_stdout(report):
                self.assertEqual(main([temporary]), 1)
            self.assertIn("Diagnostic incomplete:", report.getvalue())


if __name__ == "__main__":
    unittest.main()
