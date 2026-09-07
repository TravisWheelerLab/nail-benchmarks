import math
from collections.abc import Iterable
import re

import matplotlib as mpl
import matplotlib.pyplot as plt

TOL_ORANGE = "#EE7733"
TOL_BLUE = "#0077BB"
TOL_CYAN = "#33BBEE"
TOL_TEAL = "#009988"
TOL_RED = "#CC3311"
TOL_MAGENTA = "#EE3377"
TOL_GREY = "#BBBBBB"


TOL_VIBRANT = [
    TOL_ORANGE,
    TOL_BLUE,
    TOL_CYAN,
    TOL_TEAL,
    TOL_RED,
    TOL_MAGENTA,
    TOL_GREY,
]

COLORS = TOL_VIBRANT

MARKERS = [
    'o', 's', '^', 'v', '<', '>',
    'D', 'H', '*', 'p', 'X', 'd', 'h', '8'
]

TOOL_COLORS = {
    "hmmer": TOL_ORANGE,
    "nail": TOL_MAGENTA,
    "mmseqs": TOL_CYAN,
    "blast": TOL_TEAL,
    "diamond": TOL_BLUE,
    "last": TOL_GREY,
}

FLOAT_RE = r'\s*(-?\d+(?:\.\d+)?)\s*'


def theta(ax, x1: float, y1: float, x2: float, y2: float) -> float:
    x1, y1 = ax.transData.transform((x1, y1))
    x2, y2 = ax.transData.transform((x2, y2))

    dx = x2 - x1
    dy = y2 - y1

    theta = math.degrees(math.atan2(dy, dx))

    return theta


def annotate(
    ax,
    text: str,
    point: tuple[float, float],
    offset_px: tuple[float, float],
    color: str,
    fontsize: float = 16,
    rotation: float = 0,
    ha: str = None,
    va: str = None,
    linestyle: str = None,
    arrowstyle: str = None,
):

    # try to auto-set horizontal alignment
    if ha is None:
        if offset_px[0] == 0:
            ha = "center"
        elif offset_px[0] > 0:
            ha = "left"
        else:
            ha = "right"

    # try to auto-set vertical alignment
    if va is None:
        if offset_px[1] == 0:
            va = "center"
        elif offset_px[1] > 0:
            va = "bottom"
        else:
            va = "top"

    # compute rotation stuff
    theta = math.radians(rotation)
    ox, oy = offset_px

    roffset_px = (
        ox * math.cos(theta) - oy * math.sin(theta),
        ox * math.sin(theta) + oy * math.cos(theta),
    )

    # draw the label
    ax.annotate(
        text,
        point,
        xytext=roffset_px,
        textcoords="offset pixels",
        va=va,
        ha=ha,
        color=color,
        fontsize=fontsize,
        fontweight="bold",
        rotation=rotation,
        rotation_mode="anchor",
    )

    # (optionally) draw the line
    if linestyle is not None:
        ax.annotate(
            "",
            point,
            xytext=roffset_px,
            textcoords="offset pixels",
            arrowprops=dict(
                arrowstyle="-",
                linestyle=linestyle,
                color=color,
            ),
        )

    # (optionally) draw the arrowtip
    if arrowstyle is not None:
        head_frac = 0.12
        head_offset_px = (
            roffset_px[0] * head_frac,
            roffset_px[1] * head_frac,
        )

        ax.annotate(
            "",
            point,
            xytext=head_offset_px,
            textcoords="offset pixels",
            arrowprops=dict(
                arrowstyle=arrowstyle,
                linewidth=0,
                color=color,
            ),
        )


def prefix_label(
        prefix: str,
        exclude=[],
        reformat=[]
) -> str:
    rest, search_type = prefix.rsplit(".", 1)
    tokens = rest.split("-")

    tool = tokens[0]
    args = tokens[1:]

    args = []
    for arg in tokens[1:]:
        if any([e in arg for e in exclude]):
            continue

        args.append(arg)

    nail_params = {
        "s4.0": "--mmseqs-s 4.0",
        "s5.7": "--mmseqs-s 5.7",
        "s7.5": "--mmseqs-s 7.5",
        "s8.0": "--mmseqs-s 8.0",
        "s9.0": "--mmseqs-s 9.0",
        "s10.0": "--mmseqs-s 10.0 (default)",
        "s11.0": "--mmseqs-s 11.0",
        "s12.0": "--mmseqs-s 12.0",
        "s14.0": "--mmseqs-s 14.0",
        "ms2000": "--mmseqs-max-seqs 2000",
        "prog": "--prog-seed",
    }

    mmseqs_params = {
        "s4.0": "-s 4.0 (fast)",
        "s5.7": "-s 5.7 (default)",
        "s7.5": "-s 7.5 (sensitive)",
        "s8.0": "-s 8.0",
        "s9.0": "-s 9.0",
        "s10.0": "-s 10.0",
        "s11.0": "-s 11.0",
        "s12.0": "-s 12.0",
        "s14.0": "-s 14.0",
        "ms300": "--max-seqs 300 (default)",
        "ms2000": "--max-seqs 2000",
    }

    ret_args = []
    if tool == "hmmer":
        if search_type == "prf":
            tool = "hmmsearch (profile)"
        elif search_type == "seq":
            tool = "phmmer (sequence)"
    elif tool == "nail":
        ret_args = [nail_params[p] for p in args]
        if search_type == "prf":
            tool = "nail (profile)"
        elif search_type == "seq":
            tool = "nail (sequence)"
    elif tool == "mmseqs":
        ret_args = [mmseqs_params[p] for p in args]
        if search_type == "prf":
            tool = "mmseqs (profile)"
        elif search_type == "seq":
            tool = "mmseqs (sequence)"
    elif tool == "blast":
        if search_type == "prf":
            tool = "psiblast (profile)"
        elif search_type == "seq":
            tool = "blastp (sequence)"
    elif tool == "diamond":
        if search_type == "seq":
            tool = "diamond (sequence)"

    for (a, b) in reformat:
        ret_args = [arg.replace(a, b) for arg in ret_args]

    return tool, ret_args


def axes():
    scale = 1.75
    mpl.rcParams.update({
        "font.size": mpl.rcParams["font.size"] * scale,
    })

    fig, ax = plt.subplots(
        figsize=(16, 9),
        # constrained_layout=True
    )

    P = 0.075
    fig.subplots_adjust(
        left=P,
        right=1.0 - P,
        bottom=P,
        top=1.0 - P
    )

    return fig, ax


class Scatter:
    x: [float]
    y: [float]

    def __init__(self, lines: Iterable[str]) -> None:
        point_re = re.compile(rf"{FLOAT_RE},{FLOAT_RE}")
        self.x = []
        self.y = []
        for line in lines:
            x, y = re.search(point_re, line).groups()
            self.x.append(float(x))
            self.y.append(float(y))


class CurveND:
    data: [[float]]
    prefix: str
    tool: str
    params: dict[str, float]
    extra: [str]

    def __init__(self, line: str, n: int):
        info, _ = line.split(",(", 1)

        info_tokens = info.split(',')
        self.prefix = info_tokens[0]

        self.params = {}
        self.extra = info_tokens[1:]

        prefix_tokens = self.prefix.split('-')
        self.tool = prefix_tokens[0]

        for tok in prefix_tokens[1:]:
            m = re.match(r'([a-zA-Z]+)([0-9.]+)', tok)
            if m:
                p, v = m.groups()
                self.params[p] = float(v)
            else:
                self.extra.append(tok)

        mutli_point_re = rf"\({','.join(FLOAT_RE for _ in range(n))}\)"

        self.data = [[] for _ in range(n)]

        for pt in re.findall(mutli_point_re, line):
            for (i, v) in enumerate(pt):
                self.data[i].append(float(v))


class Curve:
    x: [float]
    y: [float]
    prefix: str
    tool: str
    search_type: str
    params: dict[str, float]
    extra: [str]

    def __init__(self, line):
        self.prefix, points = line.split(",", 1)

        params = self.prefix
        if self.prefix.endswith(".prf"):
            params = params.removesuffix(".prf")
            self.search_type = "prf"
        elif self.prefix.endswith(".seq"):
            params = params.removesuffix(".seq")
            self.search_type = "seq"

        param_tokens = params.split("-")

        self.params = {}
        self.extra = []
        self.tool = param_tokens[0]

        for tok in param_tokens[1:]:
            m = re.match(r'([a-zA-Z]+)([0-9.]+)', tok)
            if m:
                p, v = m.groups()
                self.params[p] = float(v)
            else:
                self.extra.append(tok)

        points = [
            (float(a), float(b))
            for a, b in re.findall(
                r'\(\s*(-?\d+(?:\.\d+)?)\s*,\s*(-?\d+(?:\.\d+)?)\s*\)',
                points
            )
        ]

        self.x, self.y = map(list, zip(*points))

    def approx_y(self, x: float) -> float:
        for (x1, x2), (y1, y2) in zip(
            zip(self.x, self.x[1:]),
            zip(self.y, self.y[1:])
        ):
            if x1 <= x <= x2:
                return y1 + (x - x1) * (y2 - y1) / (x2 - x1)
        raise ValueError("x out of range")


class Point:
    x: float
    y: float
    prefix: str
    tool: str
    params: dict[str, float]
    extra: [str]

    def __init__(self, line: str):
        self.prefix, point = line.split(",", 1)

        x, y = map(float, re.match(
            r'\(\s*([-+]?\d*\.?\d+)\s*,\s*([-+]?\d*\.?\d+)\s*\)', point)
            .groups()
        )
        self.x = float(x)
        self.y = float(y)

        prefix_tokens = self.prefix.rsplit(".", 1)[0].split('-')

        self.params = {}
        self.extra = []
        self.tool = prefix_tokens[0]

        for tok in prefix_tokens[1:]:
            m = re.match(r'([a-zA-Z]+)([0-9.]+)', tok)
            if m:
                p, v = m.groups()
                self.params[p] = float(v)
            else:
                self.extra.append(tok)


def partition(l: [], pred):
    a, b = [], []
    for x in l:
        (a if pred(x) else b).append(x)

    return (a, b)


def link(points: [Point], key: str):
    aa = list(set([p.params[key] for p in points]))
    aaa = [
        [
            [p.x for p in points if p.params[key] == v],
            [p.y for p in points if p.params[key] == v],
            f"{key}-{v}"
        ]
        for v in aa
    ]

    return aaa
