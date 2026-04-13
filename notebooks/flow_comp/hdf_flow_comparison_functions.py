import h5py, glob, os
from datetime import date, datetime
import re
from pydsstools.heclib.dss import HecDss
# from pydsstools.heclib.dss.HecDss import Open
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import html as _html
from plotly.io import to_html
import sys
from contextlib import contextmanager

@contextmanager
def suppress_native_output(suppress_stdout=True, suppress_stderr=True):
    """Suppress output at the file descriptor level (works for native/C/Fortran libs)."""
    devnull = os.open(os.devnull, os.O_WRONLY)
    
    old_fds = {}
    if suppress_stdout:
        old_fds[1] = os.dup(1)  # save stdout fd
        os.dup2(devnull, 1)      # redirect stdout to /dev/null
    if suppress_stderr:
        old_fds[2] = os.dup(2)  # save stderr fd
        os.dup2(devnull, 2)      # redirect stderr to /dev/null
    
    try:
        yield
    finally:
        # Restore original file descriptors
        for fd, saved in old_fds.items():
            os.dup2(saved, fd)
            os.close(saved)
        os.close(devnull)


def get_hdf_files_with_ext(inputs_filepath):
    hdf_files_list = []
    hdfs = glob.glob(str(inputs_filepath)+'\*.hdf')
    
    for hdf in hdfs: 
        # search with regex pattern
        plan_extension = re.search(r'\.p(\d{2})\.hdf',hdf)
        # Skip those that don't match and exclude plan numbers less than 10
        if not plan_extension:
            continue
        plan_num = (re.search(r'\d+',plan_extension.group())).group() #other files with the same number are related. this is important to note
        if int(plan_num) < 10:  
            continue
        
        hdf_files_list.append(hdf)
    return hdf_files_list

def todays_date():
    today = date.today()
    
    # Format to "Month dd, yyyy"
    formatted_date = today.strftime("%B %d, %Y")
    return formatted_date

def list_to_txtfile(list_to_write, out_path):
    """Subfunction to ouput text file"""
    formatted_list = '\n'.join(list_to_write) # create a string of text separated by new lines
    with open(out_path, 'w') as file:
        file.write(formatted_list)
    return

def check_for_ref_lines_in_hdfs(hdf_files_list, outputs_filepath):
    print ("Checking hdf files for reference lines...")
    plans_wo_ref_list = ["Plans without reference lines",f"{todays_date()}"]
    plans_w_ref_list = ["Plans with reference lines",f"{todays_date()}"]
    ref_ln_subgrp = r"/Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Reference Lines"
    
    for hdf in hdf_files_list:
        # open plan and check if reference lines exist
        with h5py.File(str(hdf),'r') as hf_hc:
            try: 
                ref_ln = hf_hc[ref_ln_subgrp]
                plans_w_ref_list.append(hdf)

            except KeyError:
                plans_wo_ref_list.append(hdf)
    
    # write a list of plans with and without references to text files
    plans_w_filename = os.path.join(outputs_filepath, "Plans_with_reference_lines.txt")
    list_to_txtfile(plans_w_ref_list, plans_w_filename)
    print(f"\t{len(plans_w_ref_list)-2} plans with reference lines. List exported as a .txt file in {plans_w_filename}")

    plans_wo_filename = os.path.join(outputs_filepath, "Plans_without_reference_lines.txt")
    list_to_txtfile(plans_wo_ref_list, plans_wo_filename)
    print(f"\t{len(plans_wo_ref_list)-2} plans without reference lines.List exported as a .txt file in {plans_wo_filename}")

    return plans_w_ref_list[2:]

def get_ref_ln_info_from_hdf(hdf):
    """"
    Function to open specified hdf plans and extract reference line information.
    Params: 
        hdf filepaths (list)
    Returns:
        numpy array of reference lines if exists
    """

    hdf_subgrp_basepath = r"/Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/"
    hdf_ref_ln_flows = hdf_subgrp_basepath + r"Reference Lines/Flow"
    hdf_ref_ln_names = hdf_subgrp_basepath + r"Reference Lines/Name"
    hdf_ref_ln_timestamps = hdf_subgrp_basepath + r"Time Date Stamp"

    # for hdf in hdf_files_list:
    # capture plan number from filename
    plan_ext = hdf.split(".")[-2]
    plan_num = (re.search(r'\d+',plan_ext)).group() #other files with the same number are related. this is important to note

    with h5py.File(str(hdf),'r') as hf_hc:
        try: 
            # create numpy array with collection of flows for every reference line at each time time. Units is CFS
            # time stamps for the plan from start to end

            flows = hf_hc[hdf_ref_ln_flows][:]
            ref_lines = [str(a).strip("'b\'") for a in hf_hc[hdf_ref_ln_names][:]]
            time_stamps = [str(b).strip("'b\'") for b in hf_hc[hdf_ref_ln_timestamps][:]]

        except KeyError:
            print(f"Info not found for {os.path.basename(hdf)}")

    ref_df = make_ref_df(flows, ref_lines, time_stamps)

    return plan_num, ref_df

def make_ref_df(flows, ref_lines, time_stamps):
    ref_names = [r.split("|")[0] for r in ref_lines]
    dt_list = [datetime.strptime(s, "%d%b%Y %H:%M:%S") for s in time_stamps]

    ref_df = pd.DataFrame(flows)
    ref_df.columns = ref_names
    ref_df.index = dt_list
    # ref_df.to_excel(os.path.join(outputs_filepath, "test1.xlsx"), index=True)
    return ref_df

## Get DSS filename associated with a plan
def get_dss_file_from_plan(plan_num, inputs_filepath):
    plan = glob.glob(str(inputs_filepath)+f'\*.p{plan_num}')
    with open(plan[0],'r') as file:
        content = file.read()
    perchance = re.match('(?:Plan Title=)(.*output)',content)
    dss_file_name = perchance.group(1)
    dss_file = os.path.join(inputs_filepath, "Hydrology", dss_file_name)#+".dss")
    return dss_file

## Get timestamps and flow values for each junction in the DSS file
def get_junc_info(ref_df_col, dss_file):
    junction_name = ref_df_col.split("REF_")[1]  # Junctions have the same name as the reference line minus the "REF_"
    junc_finder = f'//{junction_name}/FLOW-COMBINE/[\d\w\S]+'
    r = re.compile(junc_finder)

    with suppress_native_output():
        with HecDss.Open(str(dss_file)) as fid:
            pathname_pattern ="/*/*/*/*/*/*/"
            dss_list = fid.getPathnameList(pathname_pattern,sort=1)
            dss_matches = list(filter(r.match, dss_list))

            times_comb = []
            values_comb = []
            for i in range(0,len(dss_matches)):
                ts = fid.read_ts(dss_matches[i], trim_missing=True)
                times_comb.extend(ts.pytimes)
                # print(len(ts.pytimes))
                if i == 0:
                    values_comb = [a for a in list(ts.values)]
                    # print(len(ts.values))
                else:
                    added_values = [a for a in list(ts.values)]
                    values_comb.extend(added_values)
                    # print(len(ts.values))

    return times_comb, values_comb

def create_summary_df():
    huc10_summary_table_columns = ['Plan Number','Event Name','Reference Line Name','Reference Peak Flow (CFS)','Reference Time of Peak Flow','Junction Name','Junction Peak Flow (CFS)','Junction Time of Peak Flow','Difference in Peak Flow (Reference - Junction)','Difference in Peak Flow Time (Referemce - Junction)']
    df = pd.DataFrame(columns=huc10_summary_table_columns)
    return df

def _format_td(td):
    if pd.isna(td):
        return ""
    total = int(td.total_seconds())
    sign = "-" if total < 0 else ""
    total = abs(total)
    h, r = divmod(total, 3600)
    m, s = divmod(r, 60)
    return f"{sign}{h:02d}:{m:02d}:{s:02d}"

def parent_dir(path, levels=1):
    for _ in range(levels):
        path = os.path.dirname(path)
    return path

def get_huc_name(inputs_filepath):
    huc_folder_name = os.path.basename(inputs_filepath)
    return huc_folder_name

def create_output_dir(inputs_filepath):
    dir_loc = parent_dir(inputs_filepath, 2)
    huc_dir_name = get_huc_name(inputs_filepath)
    out_path_name = dir_loc + "\\References\\" + huc_dir_name

    os.makedirs(out_path_name, exist_ok=True)
    print(f"\tOutputs filepath: {out_path_name}")

    return out_path_name

def create_plots_and_table(inputs_filepath, outputs_filepath, huc10, plans_w_ref_list):
    if len(plans_w_ref_list) == 0:
        print(f"\tNo plans with reference lines found in {inputs_filepath}")
    else: 
        huc_df = create_summary_df()
        rows = []
        huc10 = get_huc_name(inputs_filepath)
        excel_out_path = os.path.join(outputs_filepath,f"{huc10}_flow_summary.xlsx")
        
        for plan in plans_w_ref_list:
            plan_num, ref_df = get_ref_ln_info_from_hdf(plan)
            dss_file = get_dss_file_from_plan(plan_num, inputs_filepath)

            for col in ref_df:
                junc_times, junc_values = get_junc_info(col, dss_file)
                peak_ref_time = pd.Timestamp(ref_df.index[ref_df[col] == ref_df[col].max()][0])
                peak_junc_time = pd.Timestamp(junc_times[junc_values.index(max(junc_values))])
            
                rows.append({
                'Plan Number': plan_num,
                'Event Name': os.path.basename(dss_file).split(".")[0],
                'Reference Line Name': col,
                'Reference Peak Flow (CFS)': ref_df[col].max(),
                'Reference Time of Peak Flow': peak_ref_time,
                'Junction Name': col.split("REF_")[1] if "REF_" in col else col,
                'Junction Peak Flow (CFS)': max(junc_values),
                'Junction Time of Peak Flow': peak_junc_time,
                'Difference in Peak Flow (Reference - Junction)(CFS)': ref_df[col].max() - max(junc_values),
                'Difference in Peak Flow Time (Reference - Junction)(hh:mm:ss)': peak_ref_time - peak_junc_time,
                })

            
            html_out_path = os.path.join(outputs_filepath,f"{huc10}_Plan_{plan_num}_Reference_vs_Junction_Flow_Graphs.html")
            print(f"Generating reference vs junction plots for plan {plan_num}...")
            generate_multi_figure_html(ref_df, html_out_path, huc10, dss_file)    
        
        huc_df = pd.DataFrame(rows)
        huc_df['Difference in Peak Flow Time (Reference - Junction)(hh:mm:ss)'] = huc_df['Difference in Peak Flow Time (Reference - Junction)(hh:mm:ss)'].apply(_format_td)
            
        print(f"Creating summary table for model {huc10}...")
        huc_df.to_excel(excel_out_path, index=False)
        print(f"\t{huc10} summary table exported to {excel_out_path}")
    return

###############################################################
#################### PLOTTING #################################
###############################################################


def my_fetcher(col: str, dss_file):
    """
    Calls get_junc_info(col, dss_file).
    Returns (times, values) or None if no data.
    """
    try:
        j_times, j_values = get_junc_info(col, dss_file)
        if not j_times or not j_values:
            return None
        if len(j_times) != len(j_values):
            # Guard against mismatched lengths
            # You can log/raise if this happens in practice
            n = min(len(j_times), len(j_values))
            j_times = j_times[:n]
            j_values = j_values[:n]
        # Normalize to pandas Timestamps
        j_times = [pd.to_datetime(t) for t in j_times]

        # (Optional but recommended) sort by time and drop duplicates
        # while keeping the first occurrence.
        j_df = pd.DataFrame({"t": j_times, "v": j_values})
        j_df = j_df.sort_values("t").drop_duplicates(subset="t", keep="first")
        return j_df["t"].tolist(), j_df["v"].tolist()
    except Exception:
        return None

# ------------------------------------------------------
# 2) Helpers for maxima + annotations
# ------------------------------------------------------
def _series_max_info(xs, ys):
    """
    xs: list/array-like of datetimes
    ys: list/array-like of floats
    Returns (max_value, when_timestamp)
    """
    if ys is None or len(ys) == 0:
        return (np.nan, pd.NaT)
    y_arr = np.array(ys, dtype=float)
    if y_arr.size == 0 or np.all(np.isnan(y_arr)):
        return (np.nan, pd.NaT)
    idx = int(np.nanargmax(y_arr))
    return float(y_arr[idx]), pd.to_datetime(xs[idx])

def _nice_name_with_max(base_name: str, max_val: float, max_time: pd.Timestamp) -> str:
    if pd.isna(max_val) or pd.isna(max_time):
        return base_name
    ts = pd.to_datetime(max_time).strftime("%Y-%m-%d %H:%M")
    return f"{base_name} (max={max_val:,.3f} @ {ts})"

def _annotation_text(series_name: str, max_val: float, max_time: pd.Timestamp) -> str:
    if pd.isna(max_val) or pd.isna(max_time):
        return f"{series_name}: no data"
    ts = pd.to_datetime(max_time).strftime("%Y-%m-%d %H:%M")
    return f"<b>{series_name}</b><br>Max: {max_val:,.3f} CFS<br>When: {ts}"

# ------------------------------------------------------
# 3) One chart per column (Reference vs Junction)
# ------------------------------------------------------
def build_figure_for_pair(
    ref_times: pd.DatetimeIndex,
    ref_values: pd.Series,
    ref_name: str,
    junc_times, junc_values,
    junc_name: str,
) -> go.Figure:

    # Compute max info
    ref_max, ref_when = _series_max_info(ref_times, ref_values.tolist())
    junc_max, junc_when = _series_max_info(junc_times, junc_values)

    # Legend labels include maxima
    ref_trace_name = _nice_name_with_max("Reference", ref_max, ref_when)
    junc_trace_name = _nice_name_with_max("Junction", junc_max, junc_when)

    fig = go.Figure()

    # Reference (hourly)
    fig.add_trace(go.Scatter(
        x=ref_times,
        y=ref_values,
        mode="lines",
        name=ref_trace_name,
        line=dict(color="#1f77b4", width=2),
        hovertemplate="%{x|%Y-%m-%d %H:%M}<br>Flow: %{y:,.3f} CFS<extra>Reference</extra>"
    ))

    # Junction (15-min)
    fig.add_trace(go.Scatter(
        x=junc_times,
        y=junc_values,
        mode="lines",
        name=junc_trace_name,
        line=dict(color="#d62728", width=2),
        hovertemplate="%{x|%Y-%m-%d %H:%M}<br>Flow: %{y:,.3f} CFS<extra>Junction</extra>"
    ))

    # Title & layout
    title_txt = f"{ref_name} vs {junc_name}"
    fig.update_layout(
        title=title_txt,
        xaxis_title="Time",
        yaxis_title="Flow (CFS)",
        hovermode="x unified",
        legend_title_text="Series",
        template="plotly_white",
        margin=dict(l=60, r=60, t=60, b=50),
    )

    # Side note annotation with maxima
    side_note = (
        _annotation_text("Reference", ref_max, ref_when)
        + "<br>"
        + _annotation_text("Junction", junc_max, junc_when)
    )
    fig.add_annotation(
        x=1.0, y=1.0, xref="paper", yref="paper",
        xanchor="right", yanchor="top",
        showarrow=False,
        align="right",
        bordercolor="#bbb",
        borderwidth=1,
        borderpad=6,
        bgcolor="rgba(255,255,255,0.7)",
        text=side_note,
        font=dict(size=11)
    )

    # Range slider/selectors
    fig.update_xaxes(
        rangeslider=dict(visible=True),
        rangeselector=dict(
            buttons=list([
                dict(count=1, label="1h", step="hour", stepmode="backward"),
                dict(count=6, label="6h", step="hour", stepmode="backward"),
                dict(count=12, label="12h", step="hour", stepmode="backward"),
                dict(count=1, label="1d", step="day", stepmode="backward"),
                dict(step="all")
            ])
        )
    )

    return fig

# ------------------------------------------------------
# 4) Generate a self-contained HTML with ALL charts
# ------------------------------------------------------

# Helper: safe IDs for container elements
def _make_safe_id(s: str) -> str:
    return re.sub(r'[^A-Za-z0-9_-]+', '_', str(s))

def generate_multi_figure_html(
    ref_df: pd.DataFrame,
    output_html_path: str,
    huc10: str, 
    dss_file,
):
    # Ensure datetime index
    if not isinstance(ref_df.index, pd.DatetimeIndex):
        ref_df = ref_df.copy()
        ref_df.index = pd.to_datetime(ref_df.index)

    ref_names = list(ref_df.columns)
    ref_to_id = {col: f"chart-{_make_safe_id(col)}" for col in ref_names}

    html_parts = [
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        '<meta charset="utf-8"/>',
        f"<title>{_html.escape(str(huc10))} Reference vs Junction Flows</title>",
        f"<title>{_html.escape(todays_date())}</title>"
        "<style>",
        "  body { font-family: Arial, sans-serif; margin: 0; padding: 0 10px; }",
        "  .title { margin: 20px auto 10px auto; max-width: 1200px; }",
        "  .tip { margin: 0 auto 20px auto; max-width: 1200px; color: #444; background:#f6f7f9; border:1px solid #e2e5ea; padding:12px 14px; border-radius:8px; }",
        "  .controls { margin: 0 auto 12px auto; max-width: 1200px; display: grid; grid-template-columns: 1fr; gap: 10px; }",
        "  .controls-row { display: flex; gap: 10px; align-items: center; flex-wrap: wrap; }",
        "  .search-wrap { display:flex; align-items:center; gap:8px; }",
        "  #refSearch { width: 420px; max-width: 100%; padding: 6px 8px; font-size: 14px; }",
        "  #refFilter { min-width: 420px; width: 100%; max-width: 1200px; padding: 6px; font-size: 14px; }",
        "  .actions button { padding: 6px 10px; border: 1px solid #d0d4da; background: #fff; border-radius: 6px; cursor: pointer; }",
        "  .actions button:hover { background: #f2f4f7; }",
        "  .chart { margin: 20px auto; max-width: 1200px; }",
        "</style>",
        "</head>",
        "<body>",
        f"<div class='title'><h1>{_html.escape(str(huc10))} Reference vs Junction Flows</h1></div>",

        # --- Tip at the top ---
        ("<div class='tip'>"
         "<b>How to use:</b> Use the search box to find references, then select one or more in the dropdown "
         "(hold <i>Ctrl/Cmd</i> or <i>Shift</i> to multi‑select). "
         "Click the <b>legend</b> inside each chart to toggle series on/off "
         "(double‑click a legend item to isolate it). Use the range slider and buttons to zoom."
         "</div>"),

        # --- Controls: search + multi-select + actions ---
        "<div class='controls'>",

        "<div class='controls-row search-wrap'>",
        "<label for='refSearch'><b>Search references:</b></label>",
        "<input type='text' id='refSearch' placeholder='Type to filter options (case-insensitive)…' spellcheck='false' />",
        "</div>",

        "<div class='controls-row'>",
        "<label for='refFilter'><b>Select references:</b></label>",
        # Multi-select list; show up to 12 rows (auto shrinks if fewer)
        f"<select id='refFilter' multiple size='{min(12, max(6, len(ref_names)))}'>",
        # options are appended next
    ]

    # Inject options INSIDE the <select>
    for col in ref_names:
        html_parts.append(
            f"<option value=\"{_html.escape(col)}\">{_html.escape(col)}</option>"
        )

    html_parts += [
        "</select>",
        "</div>",  # controls-row

        "<div class='controls-row actions'>",
        "<button id='selectAllVisible' type='button'>Select All Visible</button>",
        "<button id='clearSelection' type='button'>Clear Selection</button>",
        "<button id='showAllCharts' type='button'>Show All Charts</button>",
        "</div>",

        "</div>",  # .controls
    ]

    include_plotlyjs = "cdn"  # include Plotly only once

    # Build each chart container
    for col in ref_names:
        fetched = my_fetcher(col, dss_file)  # calls get_junc_info(col, dss_file) via adapter
        if fetched is None:
            html_parts.append(
                f"<div class='chart' id='{ref_to_id[col]}' data-ref='{_html.escape(col)}'>"
                f"<h3>{_html.escape(col)}: no matching junction data found</h3>"
                "</div>"
            )
            continue

        j_times, j_values = fetched

        # Junction name
        junc_name = col.split("REF_")[1] if "REF_" in col else col

        fig = build_figure_for_pair(
            ref_times=ref_df.index,
            ref_values=ref_df[col].astype(float),
            ref_name=col,
            junc_times=j_times,
            junc_values=j_values,
            junc_name=junc_name
        )

        frag = to_html(
            fig,
            full_html=False,
            include_plotlyjs=include_plotlyjs,
            default_height="600px"
        )

        html_parts.append(
            f"<div class='chart' id='{ref_to_id[col]}' data-ref='{_html.escape(col)}'>{frag}</div>"
        )
        include_plotlyjs = False

    # --- Client-side JS: search + multi-select filter ---
    html_parts += [
        "<script>",
        "  const searchEl  = document.getElementById('refSearch');",
        "  const selectEl  = document.getElementById('refFilter');",
        "  const btnAllVis = document.getElementById('selectAllVisible');",
        "  const btnClear  = document.getElementById('clearSelection');",
        "  const btnShowAll= document.getElementById('showAllCharts');",
        "  const charts    = Array.from(document.querySelectorAll('.chart'));",

        "  function visibleOption(opt) { return !opt.hidden; }",
        "  function getSelectedValues() {",
        "    return Array.from(selectEl.selectedOptions).map(o => o.value);",
        "  }",
        "  function showAllCharts() {",
        "    charts.forEach(div => {",
        "      const ref = div.getAttribute('data-ref');",
        "      if (!ref) return; // skip non-chart elements",
        "      div.style.display = '';",
        "    });",
        "  }",
        "  function applyFilter() {",
        "    const selected = getSelectedValues();",
        "    if (selected.length === 0) {",
        "      // Nothing selected => show all charts",
        "      showAllCharts();",
        "      return;",
        "    }",
        "    charts.forEach(div => {",
        "      const ref = div.getAttribute('data-ref');",
        "      if (!ref) return;",
        "      div.style.display = selected.includes(ref) ? '' : 'none';",
        "    });",
        "  }",
        "  function filterOptions() {",
        "    const q = (searchEl.value || '').toLowerCase().trim();",
        "    Array.from(selectEl.options).forEach(opt => {",
        "      const text = opt.textContent.toLowerCase();",
        "      const match = q === '' || text.includes(q);",
        "      // Use the 'hidden' property to hide options without removing them",
        "      opt.hidden = !match;",
        "    });",
        "  }",
        "  function selectAllVisible() {",
        "    Array.from(selectEl.options).forEach(opt => {",
        "      if (!opt.hidden) opt.selected = true;",
        "    });",
        "    applyFilter();",
        "  }",
        "  function clearSelection() {",
        "    Array.from(selectEl.options).forEach(opt => opt.selected = false);",
        "    applyFilter();",
        "  }",

        "  // Wire up events",
        "  searchEl.addEventListener('input', filterOptions);",
        "  selectEl.addEventListener('change', applyFilter);",
        "  btnAllVis.addEventListener('click', selectAllVisible);",
        "  btnClear.addEventListener('click', clearSelection);",
        "  btnShowAll.addEventListener('click', () => {",
        "    clearSelection(); // clears selection and shows all charts",
        "    searchEl.value = '';",
        "    filterOptions();",
        "  });",

        "  // Initial state: show all charts, no selection, no search filter",
        "  filterOptions();",
        "  applyFilter();",
        "</script>",
        "</body>",
        "</html>",
    ]

    with open(output_html_path, "w", encoding="utf-8") as f:
        f.write('\n'.join(html_parts))

    print(f"\tWrote HTML with {len(ref_names)} chart(s) to: {output_html_path}")