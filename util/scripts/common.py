def index_lines(path):
    offs, pos = [0], 0
    with open(path, 'rb') as f:
        for line in f:
            pos += len(line)
            offs.append(pos)
    return offs


class FastaIndex:
    def __init__(self, path):
        self.path = path
        self.line_offsets = index_lines(path)
        self.ranges = []
        start = 0
        length = 0

        with open(path) as f:
            for (i, line) in enumerate(f):
                if i == 0:
                    continue

                if line.startswith(">"):
                    self.ranges.append((start, i - 1, length))
                    length = 0
                    start = i
                else:
                    length += len(line) - 1

        self.ranges.append((start, i, length))

    def len(self):
        return len(self.ranges)

    def read_lines(self, start, end):
        with open(self.path, 'rb') as f:
            f.seek(self.line_offsets[start])
            return "".join([f.readline().decode() for _ in range(end - start + 1)])

    def split(self, n):
        splits = [[] for _ in range(n)]
        for i in range(self.len()):
            splits[i % n].append(self.ranges.pop())

        return splits


class HmmIndex:
    def __init__(self, path):
        self.path = path
        self.line_offsets = index_lines(path)
        self.ranges = []
        name = None
        start = 0
        length = 0

        with open(path) as f:
            for (i, line) in enumerate(f):
                if "//" in line:
                    self.ranges.append((start, i, length, name))
                    name = None
                    length = 0
                    start = i + 1

                elif line.startswith("NAME"):
                    name = line.split()[1]
                elif "LENG" in line:
                    length = int(line.split()[1])

    def len(self):
        return len(self.ranges)

    def read_lines(self, start, end):
        with open(self.path, 'rb') as f:
            f.seek(self.line_offsets[start])
            return "".join([f.readline().decode() for _ in range(end - start + 1)])

    def split(self, n):
        splits = [[] for _ in range(n)]
        for i in range(self.len()):
            splits[i % n].append(self.ranges.pop())

        return splits
