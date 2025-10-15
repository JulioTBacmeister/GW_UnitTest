from __future__ import annotations
from pathlib import Path
import re
from typing import Dict, Any, Tuple

_BOOL_MAP = {'.true.': True, '.false.': False}

def _strip_comment(line: str) -> str:
    """Remove ! comments unless inside quotes."""
    out, in_quote = [], False
    for ch in line:
        if ch in ("'", '"'):
            in_quote = not in_quote
            out.append(ch)
        elif ch == '!' and not in_quote:
            break
        else:
            out.append(ch)
    return ''.join(out).strip()

def _parse_value(raw: str) -> Any:
    s = raw.strip()
    # quoted string
    if (len(s) >= 2) and ((s[0] == s[-1] == "'") or (s[0] == s[-1] == '"')):
        return s[1:-1]
    # logicals
    low = s.lower()
    if low in _BOOL_MAP:
        return _BOOL_MAP[low]
    # Fortran D exponent → E exponent
    s_num = re.sub(r'([0-9])d([+-]?[0-9]+)', r'\1e\2', s, flags=re.IGNORECASE)
    # Try int, then float
    try:
        if re.fullmatch(r'[+-]?\d+', s_num):
            return int(s_num)
        return float(s_num)
    except ValueError:
        # Fallback: raw string (e.g., unquoted paths with spaces—rare)
        return s

def read_namelist(path: str | Path) -> Dict[str, Dict[str, Any]]:
    """
    Minimal Fortran namelist reader for the common case:
    - Groups start with &name and end with /
    - key = value pairs per line
    - Comments start with ! (outside quotes)
    """
    path = Path(path)
    txt = path.read_text()

    groups: Dict[str, Dict[str, Any]] = {}
    current: Dict[str, Any] | None = None
    current_name: str | None = None

    for raw in txt.splitlines():
        line = _strip_comment(raw)
        if not line:
            continue
        if line.lstrip().startswith('&'):
            # Start of a group
            current_name = line.lstrip()[1:].strip()
            current = {}
            groups[current_name] = current
            continue
        if line.strip() == '/':
            current = None
            current_name = None
            continue
        if current is None:
            # ignore stray lines outside groups
            continue
        # key = value (only first '=' counts)
        if '=' in line:
            key, val = line.split('=', 1)
            key = key.strip()
            val = val.strip()
            current[key] = _parse_value(val)
        else:
            # ignore lines without '=' inside a group
            pass
    return groups

def choose_active_group(nml: Dict[str, Dict[str, Any]]) -> Tuple[str, Dict[str, Any]]:
    """
    Use top_ctl_nl%calculation_type to pick the initfiles group.
    We match against each group's 'ncdata_type' (case-insensitive).
    """
    if 'top_ctl_nl' not in nml:
        raise KeyError("Missing group 'top_ctl_nl' in namelist.")
    top = nml['top_ctl_nl']
    if 'calculation_type' not in top:
        raise KeyError("Missing 'calculation_type' in group 'top_ctl_nl'.")

    calc = str(top['calculation_type']).strip().lower()

    # Search for a group with matching ncdata_type (case-insensitive)
    for gname, gvals in nml.items():
        if gname == 'top_ctl_nl':
            continue
        ncdata_type = gvals.get('ncdata_type', None)
        if ncdata_type is None:
            continue
        if str(ncdata_type).strip().lower() == calc:
            return gname, gvals

    # Fallback: try suffix match on group name like cam_initfiles_nl_<calc>
    for gname, gvals in nml.items():
        if gname == 'top_ctl_nl':
            continue
        if gname.lower().endswith('_' + calc):
            return gname, gvals

    raise ValueError(
        f"No group found with ncdata_type matching calculation_type='{calc}'. "
        "Check your namelist."
    )
