"""
analysisDashboard.py
========================================================================
Collects every plot HTML file written by the rotordynamic analysis
(rotor_geometry.html, campbell.html, ucs_map.html, mode_shapes.html,
mode_shape_3d_*.html, whirl_vs_freq_ratio.html, stiffness_vs_critical.html
-- all saved to DOWNLOADS_DIR) into one scrollable dashboard.html page.

Call build_dashboard() at the end of analyzeRotor.py's main(), after all
the individual plots have been written.

Can also be run standalone to (re)build the dashboard from whatever plot
HTML files already exist:
    python analysisDashboard.py            # build only
    python analysisDashboard.py --serve    # build and open in a local server
"""

import sys
import time
import webbrowser
from pathlib import Path

from buildRotor import DOWNLOADS_DIR

HERE = Path(DOWNLOADS_DIR)
DASHBOARD = "dashboard.html"

# (filename-or-glob, friendly title, one-line description). Globs are
# expanded and sorted so any number of mode-shape files are picked up.
PLOT_SPEC = [
    ("rotor_geometry.html",       "Rotor Geometry",
     "Shaft sections, disks and bearing locations along the axis."),
    ("campbell.html",             "Campbell Diagram",
     "Damped natural frequencies vs. shaft speed with the 1x line and critical speeds."),
    ("ucs_map.html",              "Undamped Critical Speed Map",
     "Critical speeds as a function of support stiffness."),
    ("mode_shapes.html",          "Mode Shapes (2D)",
     "Normalized lateral deflection along the rotor for each computed mode."),
    ("mode_shape_3d_*.html",      "Mode Shape (3D)",
     "Three-dimensional whirl of an individual mode."),
    ("whirl_vs_freq_ratio.html",  "Whirl vs. Frequency Ratio",
     "Forward/backward whirl amplitude swept against the excitation frequency ratio."),
    ("stiffness_vs_critical.html", "Bearing Stiffness vs. Critical Speed",
     "How the first critical speeds move as bearing stiffness is varied."),
]


def collect_plots():
    """Return an ordered list of (path, title, description) for plots that exist."""
    plots = []
    seen = set()
    for pattern, title, desc in PLOT_SPEC:
        matches = sorted(HERE.glob(pattern)) if "*" in pattern else [HERE / pattern]
        multi = len(matches) > 1 or "*" in pattern
        for path in matches:
            if path.name == DASHBOARD or not path.exists() or path in seen:
                continue
            seen.add(path)
            label = title
            if multi:
                stem = path.stem.split("_")[-1]
                label = f"{title} — {stem}"
            plots.append((path, label, desc))

    # Sweep up any other *.html we didn't name explicitly.
    for path in sorted(HERE.glob("*.html")):
        if path.name == DASHBOARD or path in seen:
            continue
        seen.add(path)
        plots.append((path, path.stem.replace("_", " ").title(), ""))
    return plots


def build_dashboard():
    plots = collect_plots()
    if not plots:
        print("No plot HTML files found - nothing to build.")
        return None

    nav_items, sections = [], []
    for idx, (path, title, desc) in enumerate(plots):
        anchor = f"plot-{idx}"
        nav_items.append(
            f'<li><a href="#{anchor}">{title}</a></li>'
        )
        desc_html = f'<p class="desc">{desc}</p>' if desc else ""
        sections.append(f"""      <section id="{anchor}" class="plot">
        <div class="plot-head">
          <h2>{title}</h2>{desc_html}
          <a class="open" href="{path.name}" target="_blank" rel="noopener">open full &#8599;</a>
        </div>
        <div class="frame-wrap">
          <iframe src="{path.name}" loading="lazy" title="{title}"></iframe>
        </div>
      </section>""")

    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>EPump Rotordynamic Analysis - All Plots</title>
<style>
  :root {{ --bg:#0f1115; --panel:#171a21; --line:#2a2f3a; --text:#e7e9ee; --muted:#9aa3b2; --accent:#5b9dff; }}
  * {{ box-sizing:border-box; }}
  body {{ margin:0; background:var(--bg); color:var(--text);
         font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }}
  header.top {{ padding:18px 24px; border-bottom:1px solid var(--line);
                position:sticky; top:0; background:var(--bg); z-index:5;
                display:flex; align-items:baseline; gap:16px; flex-wrap:wrap; }}
  header.top h1 {{ font-size:16px; margin:0; }}
  header.top .meta {{ color:var(--muted); font-size:12px; }}
  header.top .controls {{ margin-left:auto; display:flex; gap:8px; }}
  header.top button {{ background:var(--panel); color:var(--text); border:1px solid var(--line);
                       border-radius:6px; padding:6px 10px; cursor:pointer; font-size:12px; }}
  header.top button:hover {{ border-color:var(--accent); }}
  .layout {{ display:grid; grid-template-columns:240px 1fr; }}
  nav {{ border-right:1px solid var(--line); padding:16px; position:sticky; top:59px;
         align-self:start; max-height:calc(100vh - 59px); overflow:auto; }}
  nav ol {{ list-style:none; margin:0; padding:0; counter-reset:n; }}
  nav li {{ margin:2px 0; }}
  nav a {{ color:var(--muted); text-decoration:none; display:block; padding:6px 8px;
           border-radius:6px; font-size:12.5px; }}
  nav a:hover {{ color:var(--text); background:var(--panel); }}
  main {{ padding:24px; display:grid; gap:28px; }}
  main.cols-2 {{ grid-template-columns:repeat(2, minmax(0, 1fr)); }}
  .plot {{ background:var(--panel); border:1px solid var(--line); border-radius:10px;
           overflow:hidden; scroll-margin-top:76px; }}
  .plot-head {{ padding:12px 16px; border-bottom:1px solid var(--line);
                display:flex; align-items:baseline; gap:12px; flex-wrap:wrap; }}
  .plot-head h2 {{ font-size:14px; margin:0; }}
  .plot-head .desc {{ margin:0; color:var(--muted); font-size:12px; flex:1 1 300px; }}
  .plot-head .open {{ color:var(--accent); text-decoration:none; font-size:12px; white-space:nowrap; }}
  .frame-wrap {{ position:relative; width:100%; height:70vh; resize:vertical; overflow:auto; }}
  main.cols-2 .frame-wrap {{ height:52vh; }}
  iframe {{ width:100%; height:100%; border:0; background:#fff; }}
  @media (max-width:900px) {{
    .layout {{ grid-template-columns:1fr; }}
    nav {{ position:static; max-height:none; border-right:0; border-bottom:1px solid var(--line); }}
    main.cols-2 {{ grid-template-columns:1fr; }}
  }}
</style>
</head>
<body>
<header class="top">
  <h1>EPump Rotordynamic Analysis</h1>
  <span class="meta">{len(plots)} plots &middot; generated {time.strftime('%Y-%m-%d %H:%M')}</span>
  <div class="controls">
    <button id="toggle-cols">Two columns</button>
    <button id="collapse">Collapse all</button>
  </div>
</header>
<div class="layout">
  <nav>
    <ol>
{chr(10).join('      ' + item for item in nav_items)}
    </ol>
  </nav>
  <main id="grid">
{chr(10).join(sections)}
  </main>
</div>
<script>
  const main = document.getElementById('grid');
  document.getElementById('toggle-cols').addEventListener('click', (e) => {{
    main.classList.toggle('cols-2');
    e.target.textContent = main.classList.contains('cols-2') ? 'One column' : 'Two columns';
  }});
  let collapsed = false;
  document.getElementById('collapse').addEventListener('click', (e) => {{
    collapsed = !collapsed;
    document.querySelectorAll('.frame-wrap').forEach(el => el.style.display = collapsed ? 'none' : '');
    e.target.textContent = collapsed ? 'Expand all' : 'Collapse all';
  }});
</script>
</body>
</html>
"""
    out = HERE / DASHBOARD
    out.write_text(html, encoding="utf-8")
    print(f"\nDashboard written: {out}")
    for _, title, _ in plots:
        print(f"    - {title}")
    return out


def serve_dashboard(port=8000):
    """Serve HERE over HTTP and open the dashboard in a browser (Ctrl+C to stop)."""
    import functools
    import http.server
    import socketserver

    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(HERE))
    with socketserver.ThreadingTCPServer(("127.0.0.1", port), handler) as httpd:
        url = f"http://127.0.0.1:{port}/{DASHBOARD}"
        print(f"\nServing {HERE}\n  {url}\nPress Ctrl+C to stop.")
        webbrowser.open(url)
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nStopped.")


if __name__ == "__main__":
    # (Re)build the dashboard from whatever plot HTML files already exist:
    #   python analysisDashboard.py
    # Build the dashboard and open it in a local web server:
    #   python analysisDashboard.py --serve
    out = build_dashboard()
    if out and "--serve" in sys.argv:
        serve_dashboard()
