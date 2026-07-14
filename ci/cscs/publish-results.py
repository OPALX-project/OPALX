#!/usr/bin/env python3
"""Publish OPALX regression results into the opal-live-doc docs tree.

This script maps the flat output produced by the NightlyBuildX regression
runner into the directory structure expected by opal-live-doc:

    docs/opalx-regression-test/
      overview/<branch>/<arch>/index.html
      overview/<branch>/<arch>/index.org
      overview/<branch>/index.html
      overview/<branch>/index.org
      regressionTests/<branch>/<arch>/
        results_YYYY-MM-DD_HH-MM.html
        plots_YYYY-MM-DD_HH-MM/
        ok.png, nok.png, accordion.js, nightlybuildx.css

Run from inside the opal-live-doc checkout.
"""

import argparse
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

DOCS = Path("docs/opalx-regression-test")
OVERVIEW_CSS = DOCS / "overview" / "nightlybuildx.css"


BRANCH_OVERVIEW_ORG_TEMPLATE = """#+TITLE: OPALX Regression Tests - Branch '{branch}'
#+AUTHOR: AMAS - PSI
#+HTML_HEAD: <style> #content{{max-width:2400px;}}</style>
#+HTML_HEAD: <style> .document-title {{ display: none; }}</style>
#

@@html:
<header class="topbar">
  <div class="brand"><strong>OPALX</strong><span>Published Regression Results</span></div>
  <div class="scope"><span class="readonly">read-only published data</span><span class="muted mono">{branch}</span></div>
</header>
<main class="dashboard">
  <div class="page-head">
    <div>
      <h1>Branch Results</h1>
      <div class="hint">Existing published architecture pages for OPALX branch <span class="mono">{branch}</span>.</div>
    </div>
  </div>
  <div class="summary">
    <div class="metric"><span>OPALX branch</span><b>{branch}</b><strong>published data</strong></div>
    <div class="metric"><span>Mode</span><b>read-only</b><strong>no compile or run controls</strong></div>
    <div class="metric"><span>Navigation</span><b>links</b><strong>architecture pages</strong></div>
    <div class="metric"><span>Artifacts</span><b>HTML</b><strong>unit, regression, plots</strong></div>
  </div>
  <section class="section">
    <div class="section-head"><h2>Architectures</h2><span class="hint">Open an architecture overview</span></div>
@@

#+BEGIN_EXPORT html
<div class="links">
{links}
</div>
#+END_EXPORT
#+BEGIN_EXPORT html
</section></main>
#+END_EXPORT
"""

BRANCH_OVERVIEW_HTML_TEMPLATE = """<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml" lang="" xml:lang="">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=yes" />
  <meta name="author" content="AMAS - PSI" />
  <title>OPALX Regression Tests</title>
  <link rel="stylesheet" href="../nightlybuildx.css" />
</head>
<body>
<header class="topbar">
  <div class="brand"><strong>OPALX</strong><span>Published Regression Results</span></div>
  <div class="scope"><span class="readonly">read-only published data</span><span class="muted mono">{branch}</span></div>
</header>
<main class="dashboard">
  <div class="page-head">
    <div>
      <h1>Branch Results</h1>
      <div class="hint">Existing published architecture pages for OPALX branch <span class="mono">{branch}</span>.</div>
    </div>
  </div>
  <div class="summary">
    <div class="metric"><span>OPALX branch</span><b>{branch}</b><strong>published data</strong></div>
    <div class="metric"><span>Mode</span><b>read-only</b><strong>no compile or run controls</strong></div>
    <div class="metric"><span>Navigation</span><b>links</b><strong>architecture pages</strong></div>
    <div class="metric"><span>Artifacts</span><b>HTML</b><strong>unit, regression, plots</strong></div>
  </div>
  <section class="section">
    <div class="section-head"><h2>Architectures</h2><span class="hint">Open an architecture overview</span></div>
    <div class="links">
{links}
    </div>
  </section>
</main>
</body>
</html>
"""

ARCH_OVERVIEW_ORG_TEMPLATE = """#+TITLE: OPALX Regression Tests - {arch}
#+AUTHOR: AMAS - PSI
#+HTML_HEAD: <style> #content{{max-width:2400px;}}</style>
#+HTML_HEAD: <style> .document-title {{ display: none; }}</style>
#

@@html:
<header class="topbar">
  <div class="brand"><strong>OPALX</strong><span>Published Regression Results</span></div>
  <div class="scope"><span class="readonly">read-only published data</span><span class="muted mono">{branch} / {arch}</span></div>
</header>
<main class="dashboard">
  <div class="page-head">
    <div>
      <h1>Architecture Results</h1>
      <div class="hint">Existing published runs for branch <span class="mono">{branch}</span> on <span class="mono">{arch}</span>.</div>
    </div>
  </div>
  <div class="summary">
    <div class="metric"><span>OPALX branch</span><b>{branch}</b><strong>published data</strong></div>
    <div class="metric"><span>Architecture</span><b>{arch}</b><strong>result history</strong></div>
    <div class="metric"><span>Mode</span><b>read-only</b><strong>no compile or run controls</strong></div>
    <div class="metric"><span>Artifacts</span><b>HTML</b><strong>unit, regression, plots</strong></div>
  </div>
  <section class="section">
    <div class="section-head"><h2>Published Runs</h2><span class="hint">Newest rows first</span></div>
@@

| Date | build | unit-tests | regression-tests |
|------------+-------+--------------+--------------------|
{rows}
|            |       |              |                    |
#+BEGIN_EXPORT html
</section></main>
#+END_EXPORT
"""

ARCH_OVERVIEW_HTML_TEMPLATE = """<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml" lang="" xml:lang="">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=yes" />
  <meta name="author" content="AMAS - PSI" />
  <title>OPALX - {arch}</title>
  <link rel="stylesheet" href="../../nightlybuildx.css" />
</head>
<body>
<header class="topbar">
  <div class="brand"><strong>OPALX</strong><span>Published Regression Results</span></div>
  <div class="scope"><span class="readonly">read-only published data</span><span class="muted mono">{branch} / {arch}</span></div>
</header>
<main class="dashboard">
  <div class="page-head">
    <div>
      <h1>Architecture Results</h1>
      <div class="hint">Existing published runs for branch <span class="mono">{branch}</span> on <span class="mono">{arch}</span>.</div>
    </div>
  </div>
  <div class="summary">
    <div class="metric"><span>OPALX branch</span><b>{branch}</b><strong>published data</strong></div>
    <div class="metric"><span>Architecture</span><b>{arch}</b><strong>result history</strong></div>
    <div class="metric"><span>Mode</span><b>read-only</b><strong>no compile or run controls</strong></div>
    <div class="metric"><span>Artifacts</span><b>HTML</b><strong>unit, regression, plots</strong></div>
  </div>
  <section class="section">
    <div class="section-head"><h2>Published Runs</h2><span class="hint">Newest rows first</span></div>
    <div class="result-table-scroll">
      <table>
        <thead>
          <tr><th>Date</th><th>Regression Results</th></tr>
        </thead>
        <tbody>
{rows}
        </tbody>
      </table>
    </div>
    <p><a href="../index.html">Back to branch overview</a></p>
  </section>
</main>
</body>
</html>
"""


def discover_publication_timestamp(source_dir: Path) -> str:
    """Find a YYYY-MM-DD_HH-MM style timestamp from the source directory.

    If the regression runner was invoked with a git short SHA as timestamp,
    the files will not carry a parseable datetime.  In that case we fall back
    to the file modification time of the newest result HTML, and finally to
    the current UTC time.
    """
    results = sorted(source_dir.glob("results_*.html"))
    if results:
        name = results[0].stem  # 'results_YYYY-MM-DD_HH-MM' or 'results_<sha>'
        m = re.search(r"(\d{4}-\d{2}-\d{2}_\d{2}-\d{2})", name)
        if m:
            return m.group(1)
        # fall back to mtime
        mtime = datetime.fromtimestamp(results[0].stat().st_mtime, tz=timezone.utc)
        return mtime.strftime("%Y-%m-%d_%H-%M")
    # No results at all: use current time.
    return datetime.now(tz=timezone.utc).strftime("%Y-%m-%d_%H-%M")


def copy_results(source_dir: Path, target_dir: Path, timestamp: str) -> None:
    """Copy regression result HTML and plots into the docs tree.

    Renames files to the canonical timestamp so that repeated runs with a
    git-SHA timestamp still produce chronologically named artifacts.
    """
    target_dir.mkdir(parents=True, exist_ok=True)

    # Helper: copy a source file to a canonical name in target_dir.
    def install(src: Path, dst_name: str) -> None:
        if src.is_file():
            shutil.copy2(src, target_dir / dst_name)

    # Copy the results HTML.  If the runner used a non-datetime timestamp,
    # rename to the canonical results_YYYY-MM-DD_HH-MM.html form.
    for result_html in sorted(source_dir.glob("results_*.html")):
        install(result_html, f"results_{timestamp}.html")
        break  # only one expected

    # Copy the plots directory.  Same renaming logic.
    plots_dirs = sorted(source_dir.glob("plots_*"))
    if plots_dirs:
        src_plots = plots_dirs[0]
        dst_plots = target_dir / f"plots_{timestamp}"
        if dst_plots.exists():
            shutil.rmtree(dst_plots)
        shutil.copytree(src_plots, dst_plots)

    # Static assets required by the results HTML.
    for asset in ("ok.png", "nok.png", "accordion.js"):
        src = source_dir / asset
        if src.is_file():
            shutil.copy2(src, target_dir / asset)

    # nightlybuildx.css lives in docs/opalx-regression-test/overview/
    if OVERVIEW_CSS.is_file():
        shutil.copy2(OVERVIEW_CSS, target_dir / "nightlybuildx.css")


def generate_arch_overview(branch: str, arch: str, regtests_dir: Path,
                           overview_dir: Path) -> None:
    """Generate the architecture-specific overview pages."""
    overview_dir.mkdir(parents=True, exist_ok=True)

    rows_org = []
    rows_html = []

    # Collect runs sorted newest first.
    result_files = sorted(regtests_dir.glob("results_*.html"), reverse=True)
    for result_file in result_files:
        m = re.search(r"results_(\d{4}-\d{2}-\d{2}_\d{2}-\d{2})",
                      result_file.name)
        if not m:
            continue
        ts = m.group(1)
        display = ts.replace("_", " ")
        rel = f"../../../regressionTests/{branch}/{arch}/{result_file.name}"
        label = f"[[{rel}][{arch} run {display}]]"
        rows_org.append(f"| {display} | - | - | {label} |")
        rows_html.append(
            f'          <tr><td>{display}</td><td><a href="{rel}">{arch} run {display}</a></td></tr>'
        )

    if not rows_org:
        rows_org.append("| - | - | - | No published runs yet |")
        rows_html.append("          <tr><td>-</td><td>No published runs yet</td></tr>")

    org = ARCH_OVERVIEW_ORG_TEMPLATE.format(
        branch=branch, arch=arch, rows="\n".join(rows_org))
    html = ARCH_OVERVIEW_HTML_TEMPLATE.format(
        branch=branch, arch=arch, rows="\n".join(rows_html))

    (overview_dir / "index.org").write_text(org)
    (overview_dir / "index.html").write_text(html)


def generate_branch_overview(branch: str, overview_branch_dir: Path) -> None:
    """Regenerate the branch overview listing all architecture subdirs."""
    overview_branch_dir.mkdir(parents=True, exist_ok=True)

    archs = []
    if overview_branch_dir.is_dir():
        for item in sorted(overview_branch_dir.iterdir()):
            if item.is_dir() and (item / "index.html").is_file():
                archs.append(item.name)

    if not archs:
        archs.append("no-architectures-yet")

    org_links = "\n".join(
        f'<a class="result-link" href="./{arch}/index.html"><span>{arch}</span><span class="mono">Open overview</span></a>'
        for arch in archs
    )
    html_links = "\n".join(
        f'      <a class="result-link" href="./{arch}/index.html"><span>{arch}</span><span class="mono">Open overview</span></a>'
        for arch in archs
    )

    (overview_branch_dir / "index.org").write_text(
        BRANCH_OVERVIEW_ORG_TEMPLATE.format(branch=branch, links=org_links))
    (overview_branch_dir / "index.html").write_text(
        BRANCH_OVERVIEW_HTML_TEMPLATE.format(branch=branch, links=html_links))


def publish(source_dir: str, branch: str, arch: str, commit_sha: str) -> None:
    src = Path(source_dir)
    if not src.is_dir():
        print(f"ERROR: source directory not found: {src}", file=sys.stderr)
        sys.exit(1)

    timestamp = discover_publication_timestamp(src)

    regtests_dir = DOCS / "regressionTests" / branch / arch
    overview_arch_dir = DOCS / "overview" / branch / arch
    overview_branch_dir = DOCS / "overview" / branch

    copy_results(src, regtests_dir, timestamp)
    generate_arch_overview(branch, arch, regtests_dir, overview_arch_dir)
    generate_branch_overview(branch, overview_branch_dir)

    print(f"Published {arch} results for branch {branch} (timestamp {timestamp}, commit {commit_sha})")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Publish OPALX regression results into opal-live-doc.")
    parser.add_argument("--source-dir", required=True,
                        help="Local regression output directory")
    parser.add_argument("--branch", required=True,
                        help="Target branch (e.g. master)")
    parser.add_argument("--arch", required=True,
                        help="Architecture identifier")
    parser.add_argument("--commit-sha", default="unknown",
                        help="Short commit SHA for the commit message")
    args = parser.parse_args()

    publish(args.source_dir, args.branch, args.arch, args.commit_sha)


if __name__ == "__main__":
    main()
