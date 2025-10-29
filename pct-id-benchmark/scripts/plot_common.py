from dataclasses import dataclass

TOL_VIBRANT = [
    "#EE7733",  # orange
    "#0077BB",  # blue
    "#33BBEE",  # cyan
    "#009988",  # teal
    "#CC3311",  # red
    "#EE3377",  # magenta
    "#BBBBBB",  # grey
]

COLORS = TOL_VIBRANT

MARKERS = ['o', 's', '^', 'v', '<', '>', 'D', 'H', '*', 'p', 'X', 'd', 'h', '8']

TOOLS = [
    'diamond',
    'last',
    'blast',
    'hmmer',
    'mmseqs',
    'nail',
]


@dataclass
class Style:
    label: str
    marker: str
    color: str
    filled: bool = True


styles = {
    #
    "diamond seq": Style("diamond (seq)", "o", COLORS[0]),
    #
    "last seq": Style("last (seq)", "o", COLORS[1]),
    #
    "blast prf": Style("psiblast", "D", COLORS[2]),
    "blast seq": Style("blastp", "o"  , COLORS[2]),
    #
    "hmmer prf": Style("hmmsearch", "D", COLORS[3]),
    "hmmer seq": Style("phmmer", "o", COLORS[3]),
    #
    "nail-s12.0 prf": Style("nail (prf)", "D", COLORS[4]),
    "nail-s12.0 seq": Style("nail (seq)", "o", COLORS[4]),
    #
    "mmseqs-default prf": Style("mmseqs2 (prf) | default", "D", COLORS[5], filled=False),
    "mmseqs-sens prf"   : Style("mmseqs2 (prf) | -s 7.5", "D", COLORS[5]),
    "mmseqs-nail prf"   : Style("mmseqs2 (prf) | -s 12.0 --max-seqs 2000", "<", COLORS[5], filled=False),
    #
    "mmseqs-default seq": Style("mmseqs2 (seq) | default", "o", COLORS[1], filled=False),
    "mmseqs-s7.5 seq"   : Style("mmseqs2 (seq) | -s 7.5", "o", COLORS[1]),
    "mmseqs-s12.0 seq"  : Style("mmseqs2 (seq) | -s 12.0 --max-seqs 2000", "<", COLORS[1], filled=False),
}

_marker_seen = {}
_marker_idx = 0


def marker_for(s):
    global _marker_idx
    if s not in _marker_seen:
        _marker_seen[s] = MARKERS[_marker_idx % len(MARKERS)]
        _marker_idx += 1
    return _marker_seen[s]


_color_seen = {}
_color_idx = 0


def color_for(s):
    global _color_idx
    if s not in _color_seen:
        _color_seen[s] = COLORS[_color_idx % len(MARKERS)]
        _color_idx += 1
    return _color_seen[s]


def _style(label):
    if label not in styles:
        return ({'label': label}, None)

    prefix, search_type = label.split(" ")
    style = styles[label]

    d = {
        'label': style.label,
        'marker': style.marker,
        'color': style.color,
    }

    return (d, style)


def plot_style(label):
    (d, style) = _style(label)

    if style is not None:
        if not style.filled:
            d['mfc'] = 'white'

    return d


def scatter_style(label):
    (d, style) = _style(label)

    if style is not None:
        if not style.filled:
            d['facecolors'] = 'white'

    return d


def parse_point(point):
    label, rest = point.split(",", 1)
    p = rest.strip().lstrip("(").rstrip(")")
    x, y = p.split(",")
    x, y = float(x), float(y)

    return label.strip(), x, y


def parse_curve(line):
    label, rest = line.split(",", 1)
    pairs = rest.strip().split("),")
    xs = []
    ys = []
    for p in pairs:
        p = p.strip().lstrip("(").rstrip(")")
        if p:
            x, y = p.split(",")
            x, y = float(x), float(y)
            xs.append(x)
            ys.append(y)

    return label.strip(), xs, ys
