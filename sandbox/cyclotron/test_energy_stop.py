#!/usr/bin/env python3
"""Single-rank end-to-end rejection tests for TRACK EKINSTOP.

Uses a generated successful short-target input as the fixture. No production
input is overwritten; all outputs live in fresh per-case directories.
"""
import argparse
import os
from pathlib import Path
import re
import subprocess


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', type=Path, default=root/'opalx/acceleration590/opalx-stop-smoke-v2')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    template = (args.fixture/'acceleration.in').read_text()
    # Fixture-relative field/profile/distribution paths must remain valid when
    # the generated rejection deck is executed in its separate output directory.
    def absolute_input(match):
        path = (args.fixture / match.group(1)).resolve()
        return f'"{path}"' if path.is_file() else match.group(0)
    template = re.sub(r'"([^"\n]+)"', absolute_input, template)
    cases = [
        ('negative', re.sub(r'EKINSTOP=[^,;]+', 'EKINSTOP=-1', template), 'EKINSTOP must be finite and positive'),
        ('below_launch', re.sub(r'EKINSTOP=[^,;]+', 'EKINSTOP=0.071', template), 'EKINSTOP must exceed the launch'),
        ('turns_conflict', template.replace('FIELDSOLVER=FS0;', 'FIELDSOLVER=FS0, TURNS=1;'), 'no explicit TURNS'),
        ('budget', re.sub(r'MAXSTEPS=\d+', 'MAXSTEPS=10', template), 'EKINSTOP was not reached'),
        ('no_gaps', template.replace(',RF0,RF1,RF2,RF3,RF4)', ')'), 'EKINSTOP requires the single-container cyclotron gap path'),
        ('zero_target', re.sub(r'EKINSTOP=[^,;]+', 'EKINSTOP=0', template), 'EKINSTOP must be finite and positive'),
        ('multiple_segments', re.sub(r'DT=[^,;]+', 'DT={1e-11,1e-11}', template), 'one positive DT segment'),
    ]
    for name, text, message in cases:
        work = args.output.resolve()/name
        work.mkdir(parents=True, exist_ok=False)
        (work/'case.in').write_text(text)
        with (work/'run.log').open('w') as log:
            result = subprocess.run([str(root.parents[1]/'omp-build/src/opalx'), '--info', '1', 'case.in'],
                cwd=work, env=dict(os.environ, OMP_NUM_THREADS='1'), stdout=log, stderr=subprocess.STDOUT)
        assert result.returncode != 0, name
        assert message in (work/'run.log').read_text(), (name, message)
        print(f'{name}: PASS', flush=True)


if __name__ == '__main__':
    main()
