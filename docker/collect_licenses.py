#!/usr/bin/env python3
"""Collect the license texts that go with the opalx executable.

Usage: collect_licenses.py <source dir> <build dir> <output dir> <git revision>

Writes to <output dir>:
  README         where the source is, and what the other files are
  LICENSE        OPALX's license
  BSD-NOTICES    notices of the OPALX source files under a BSD license
  third-party/   the license of each library compiled into opalx
"""

import os
import re
import shutil
import sys
import textwrap

# Libraries compiled into opalx, with their license files relative to <build dir>/_deps.
# A missing file stops the image build, so a moved license is noticed.
THIRD_PARTY = [
    ('ippl', 'ippl-src/LICENSE'),
    ('kokkos', 'ippl-build/_deps/kokkos-src/LICENSE'),
    ('heffte', 'ippl-build/_deps/heffte-src/LICENSE'),
    ('hdf5', 'hdf5-src/LICENSE'),
    ('h5hut', 'h5hut-src/COPYING'),
    ('lapack', 'opalx_lapack-src/LICENSE'),
    ('lapacke', 'opalx_lapack-src/LAPACKE/LICENSE'),
]

README = """\
OPALX, CPU (OpenMP) build

Source code: https://github.com/OPALX-project/OPALX
Revision:    {revision}
The versions and download locations of the libraries in third-party/ are set
by the CMake files in cmake/ at that revision.

LICENSE        OPALX, GNU General Public License version 3 or later
BSD-NOTICES    OPALX source files that are under a BSD license instead
third-party/   libraries compiled into opalx: {names}

The licenses of the Ubuntu packages in this image are in
/usr/share/doc/<package>/copyright.
"""


def comment_text(comment):
    """The text of a /* ... */ comment, without the comment markers."""
    lines = [re.sub(r'^\s*(/\*+|\*+/|\*)', '', line) for line in comment.splitlines()]
    lines = [re.sub(r'\s*\*+/\s*$', '', line).rstrip() for line in lines]
    return textwrap.dedent('\n'.join(lines)).strip()


def bsd_notices(src):
    """The BSD license comments in OPALX's source, each with the files it covers."""
    notices = {}
    for root, dirs, files in os.walk(os.path.join(src, 'src')):
        dirs.sort()
        for name in sorted(files):
            path = os.path.join(root, name)
            with open(path, errors='replace') as f:
                comments = re.findall(r'/\*.*?\*/', f.read(), re.S)
            for comment in comments:
                if 'Redistribution and use in source' in comment:
                    notices.setdefault(comment_text(comment), []).append(os.path.relpath(path, src))
    return notices


def main(src, build, out, revision):
    os.makedirs(out, exist_ok=True)
    shutil.copy(os.path.join(src, 'LICENSE'), os.path.join(out, 'LICENSE'))

    for name, path in THIRD_PARTY:
        os.makedirs(os.path.join(out, 'third-party', name), exist_ok=True)
        shutil.copy(os.path.join(build, '_deps', path), os.path.join(out, 'third-party', name))

    with open(os.path.join(out, 'BSD-NOTICES'), 'w') as f:
        f.write('These OPALX source files are under a BSD license instead of the GPL.\n')
        for notice, files in bsd_notices(src).items():
            f.write('\n' + '=' * 78 + '\n\n' + '\n'.join(files) + '\n\n' + notice + '\n')

    with open(os.path.join(out, 'README'), 'w') as f:
        f.write(README.format(revision=revision, names=', '.join(name for name, _ in THIRD_PARTY)))


if __name__ == '__main__':
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
