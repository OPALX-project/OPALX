#!/usr/bin/env python3
"""Generate the OPALX regression test master overview page.

Scans the opal-live-doc repository for architecture directories that contain
an index.html, then creates overview/<branch>/index.html with links to each
architecture's detailed results.

Run from inside the opal-live-doc checkout.
"""

import sys
from pathlib import Path


TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>OPALX Regression Tests - Branch '{branch}'</title>
<style>
body {{
  font-family: "Helvetica Neue", Helvetica, Arial, sans-serif;
  max-width: 1200px;
  margin: 2em auto;
  padding: 0 1em;
  line-height: 1.6;
}}
h1 {{
  border-bottom: 2px solid #ddd;
  padding-bottom: 0.3em;
}}
h2 {{
  margin-top: 2em;
  border-bottom: 1px solid #eee;
  padding-bottom: 0.3em;
}}
.branch-name {{
  font-family: "Consolas", "Monaco", "Courier New", monospace;
  font-size: 0.9em;
  background-color: #f4f4f4;
  padding: 2px 8px;
  border-radius: 3px;
  color: #c7254e;
  font-weight: 600;
  border: 1px solid #ddd;
}}
ul.architectures {{
  list-style: none;
  padding: 0;
}}
ul.architectures li {{
  margin: 0.5em 0;
  padding: 0.8em;
  background: #fafafa;
  border: 1px solid #eee;
  border-radius: 4px;
}}
ul.architectures a {{
  font-family: "Consolas", "Monaco", "Courier New", monospace;
  text-decoration: none;
  color: #337ab7;
}}
ul.architectures a:hover {{
  text-decoration: underline;
}}
</style>
</head>
<body>
<h1>Overview - All Architectures</h1>
<p>This page provides an overview of regression test results across all
architectures for branch <span class="branch-name">{branch}</span>.</p>

<h2>Architecture Quick Links</h2>
{links}

</body>
</html>
"""


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <branch>", file=sys.stderr)
        sys.exit(1)

    branch = sys.argv[1]
    repo_root = Path.cwd()

    # Find directories (other than .git and overview) that have an index.html.
    architectures = []
    for item in sorted(repo_root.iterdir()):
        if item.is_dir() and item.name not in {".git", "overview"}:
            if (item / "index.html").is_file():
                architectures.append(item.name)

    if architectures:
        link_items = "\n".join(
            f'<li><a href="../../{arch}/index.html">{arch}</a></li>'
            for arch in architectures
        )
        links = f'<ul class="architectures">\n{link_items}\n</ul>'
    else:
        links = "<p>No architectures available yet.</p>"

    overview_dir = repo_root / "overview" / branch
    overview_dir.mkdir(parents=True, exist_ok=True)
    overview_html = overview_dir / "index.html"
    overview_html.write_text(TEMPLATE.format(branch=branch, links=links))
    print(f"Generated {overview_html}")
    for arch in architectures:
        print(f"  - {arch}")


if __name__ == "__main__":
    main()
