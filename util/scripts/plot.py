import re


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

TOOL_COLORS = {
    "diamond": COLORS[0],
    "last": COLORS[5],
    "blast": COLORS[2],
    "mmseqs": COLORS[1],
    "hmmer": COLORS[3],
    "nail": COLORS[4],
}

FLOAT_RE = r'\s*(-?\d+(?:\.\d+)?)\s*'


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

        point_re = r"\({}\)".format(",".join(FLOAT_RE for _ in range(n)))

        self.data = [[] for _ in range(n)]

        for pt in re.findall(point_re, line):
            for (i, v) in enumerate(pt):
                self.data[i].append(float(v))


class Curve:
    x: [float]
    y: [float]
    prefix: str
    tool: str
    params: dict[str, float]
    extra: [str]

    def __init__(self, line):
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

        points = [
            (float(a), float(b)) for a, b in re.findall(
                r'\(\s*(-?\d+(?:\.\d+)?)\s*,\s*(-?\d+(?:\.\d+)?)\s*\)', line
            )]

        self.x, self.y = map(list, zip(*points))


class Point:
    x: float
    y: float
    prefix: str
    tool: str
    params: dict[str, float]
    extra: [str]

    def __init__(self, line: str):
        self.prefix, search_type, point = line.split(",", 2)

        x, y = map(float, re.match(r'\(\s*([-+]?\d*\.?\d+)\s*,\s*([-+]?\d*\.?\d+)\s*\)', point).groups())
        self.x = float(x)
        self.y = float(y)

        prefix_tokens = self.prefix.split('-')

        self.params = {}
        self.extra = [search_type]
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
