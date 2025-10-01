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

TOOLS = [
    'diamond',
    'blast',
    'hmmer',
    'mmseqs',
    'nail',
]


@dataclass
class Style:
    label: str
    marker: str
    filled: bool = True


styles = {
    #
    "diamond seq": Style("diamond (seq)", "o"),
    #
    "blast prf": Style("psiblast", "D"),
    "blast seq": Style("blastp", "o"),
    #
    "hmmer prf": Style("hmmsearch", "D"),
    "hmmer seq": Style("phmmer", "o"),
    #
    "mmseqs-p7 prf": Style("mmseqs2 (p7 prf) | -s 7.5", "*"),
    "mmseqs prf": Style("mmseqs2 (prf) | -s 7.5", "D"),
    "mmseqs seq": Style("mmseqs2 (seq) | -s 7.5", "o"),
    #
    "nail-dbl prf": Style("nail (profile; double seed)", "D"),
    "nail-dbl seq": Style("nail (sequence; double seed)", "o"),
    "nail prf": Style("nail (prf)", "D"),
    "nail seq": Style("nail (seq)", "o"),
    #
    "nail-k80 prf": Style("nail (prf) | default", "D"),
    "nail-k60 prf": Style("nail (prf) | --k-score 60", "v"),
    "nail-dbl-k80 prf": Style("nail (prf) | --double-seed", "D", filled=False),
    "nail-dbl-k60 prf": Style("nail (prf) | --k-score 60 --double-seed", "v", filled=False),
    "nail-k80 seq": Style("nail (seq) | default", "o"),
    "nail-k60 seq": Style("nail (seq) | --k-score 60", "v"),
    "nail-dbl-k80 seq": Style("nail (seq) | --double-seed", "o", filled=False),
    "nail-dbl-k60 seq": Style("nail (seq) | --k-score 60 --double-seed", "v", filled=False),
}


def _style(label):
    prefix, search_type = label.split(" ")
    tool = prefix.split('-')[0]
    style = styles[label]
    color = COLORS[TOOLS.index(tool)]

    return {
        'label': style.label,
        'marker': style.marker,
        'color': color,
    }


def plot_style(label):
    style = styles[label]
    d = _style(label)

    if not style.filled:
        d['mfc'] = 'white'

    return d


def scatter_style(label):
    style = styles[label]
    d = _style(label)

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
