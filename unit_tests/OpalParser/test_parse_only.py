"""CLI integration checks; each invocation runs in an isolated temporary directory."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

EXECUTABLE = str(Path(sys.argv.pop(1)).resolve())


class ParseOnlyTest(unittest.TestCase):
    def run_input(self, text, *, files=None, startup='', valid=True, diagnostic=None,
                  arguments=None):
        with tempfile.TemporaryDirectory(prefix='opalx-parse-') as tmp:
            work = Path(tmp)
            home = work / 'home'
            home.mkdir()
            (home / 'init.opal').write_text(startup)
            (work / 'input.in').write_text(text)
            (work / 'errormsg.txt').write_text('preserve me')
            for name, content in (files or {}).items():
                (work / name).write_text(content)
            before = {str(p.relative_to(work)): p.read_bytes()
                      for p in work.rglob('*') if p.is_file()}
            result = subprocess.run(
                [EXECUTABLE, *(arguments or ['--parse-only', '--input', 'input.in'])],
                cwd=work, env={**os.environ, 'HOME': str(home)},
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=30)
            output = result.stdout
            self.assertEqual(result.returncode == 0, valid, output)
            self.assertEqual('Syntax validation succeeded' in output, valid, output)
            if diagnostic:
                self.assertIn(diagnostic, output)
            after = {str(p.relative_to(work)): p.read_bytes()
                     for p in work.rglob('*') if p.is_file()}
            self.assertEqual(before, after, output)
            self.assertEqual(sorted(str(p.relative_to(work)) for p in work.rglob('*')
                                    if p.is_dir()), ['home'], output)

    def test_definitions_and_tracks(self):
        # Missing runtime references and files are deliberately outside this mode.
        self.run_input('''REAL length = 1;
D: DRIFT, L=length;
L: LINE=(D);
TRACK, LINE=L, BEAM="missing";
RUN, FIELDSOLVER="missing";
ENDTRACK;
TRACK;
RUN;
ENDTRACK;
''')

    def test_external_data_and_physics_are_not_validated(self):
        self.run_input('''C: RFCAVITY, L=-1, FMAPFN="missing.map";
D: DISTRIBUTION, TYPE=FROMFILE, FNAME="missing.dat";
B: BEAM, PARTICLE=ELECTRON;
COF, LINE="missing", BEAM="missing";
DUMPEMFIELDS, FILE_NAME="forbidden.dat";
TRACK, LINE="missing", BEAM="missing"; RUN; ENDTRACK;
''')

    def test_unknown_command(self):
        self.run_input('NOTACOMMAND;', valid=False, diagnostic='unknown')

    def test_unknown_attribute(self):
        self.run_input('D: DRIFT, TYPO=1;', valid=False, diagnostic='TYPO')

    def test_expression(self):
        self.run_input('REAL x = 1 + ;', valid=False)

    def test_run_attribute(self):
        self.run_input('TRACK; RUN, TYPO=1; ENDTRACK;', valid=False, diagnostic='TYPO')

    def test_structured_statement_in_track(self):
        self.run_input('TRACK; IF (TRUE) { RUN; } ENDTRACK;')

    def test_run_outside_track(self):
        self.run_input('RUN;', valid=False)

    def test_unclosed_track(self):
        self.run_input('TRACK; RUN;', valid=False, diagnostic='ENDTRACK')

    def test_quit_does_not_close_track(self):
        self.run_input('TRACK; QUIT;', valid=False, diagnostic='ENDTRACK')

    def test_missing_semicolon(self):
        self.run_input('REAL x = 1', valid=False, diagnostic="terminated")

    def test_missing_bracket(self):
        self.run_input('L: LINE=(D;', valid=False, diagnostic='bracket')

    def test_missing_compound_brace(self):
        self.run_input('{ REAL x = 1;', valid=False, diagnostic="missing '}'")

    def test_missing_if_body(self):
        self.run_input('IF (TRUE)', valid=False, diagnostic='requires a statement')

    def test_missing_else_body(self):
        self.run_input('IF (TRUE) REAL x=1; ELSE', valid=False,
                       diagnostic='requires a statement')

    def test_missing_while_body(self):
        self.run_input('WHILE (FALSE)', valid=False, diagnostic='requires a statement')

    def test_include(self):
        self.run_input('CALL, FILE="child.in"; REAL y=x+1;',
                       files={'child.in': 'REAL x=2;'})

    def test_include_error(self):
        self.run_input('CALL, FILE="child.in";',
                       files={'child.in': 'BOGUS;'}, valid=False, diagnostic='child.in')

    def test_missing_include(self):
        self.run_input('CALL, FILE="missing.in";', valid=False)

    def test_macro_and_condition(self):
        self.run_input('''check(): MACRO { TRACK; RUN; ENDTRACK; };
IF (x == 2) check(); ELSE BOGUS;
''', startup='REAL x=2;')

    def test_macro_error(self):
        self.run_input('check(): MACRO { BOGUS; }; check();', valid=False)

    def test_startup_simulation_suppressed(self):
        self.run_input('REAL x=1;', startup='TRACK; RUN; ENDTRACK;')

    def test_option_and_loop(self):
        self.run_input('OPTION, ECHO=FALSE; REAL i=0; '
                       'WHILE (i < 2) { i=i+1; TRACK; RUN; ENDTRACK; }')

    def test_missing_input(self):
        self.run_input('', arguments=['--parse-only', 'missing.in'], valid=False)

    def test_restart_needs_no_checkpoint(self):
        self.run_input('TRACK; RUN; ENDTRACK;',
                       arguments=['--input', 'input.in', '--parse-only', '--restart', 'missing.h5'])

    def test_shell_commands_suppressed(self):
        self.run_input('SYSTEM, "touch forbidden"; PSYSTEM, "touch forbidden2";')

    def test_positional_and_trailing_flag(self):
        self.run_input('TRACK; RUN; ENDTRACK;',
                       arguments=['input.in', '--parse-only'])


if __name__ == '__main__':
    unittest.main()
