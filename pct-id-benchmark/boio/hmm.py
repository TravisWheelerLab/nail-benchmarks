from dataclasses import dataclass, field, fields
from enum import Enum
from typing import Dict, List, Optional
import math


class Alphabet(Enum):
    AMINO = "amino"
    DNA = "dna"


@dataclass
class HmmHeader:
    name: Optional[str] = None
    length: Optional[int] = None
    n_seq: Optional[int] = None
    alphabet: Optional[Alphabet] = None

    def complete(self) -> bool:
        return all(getattr(self, f.name) is not None for f in fields(self))


@dataclass
class HmmRecord:
    header: HmmHeader = field(default_factory=HmmHeader)
    mat_prob: List[List[float]] = field(default_factory=list)
    ins_prob: List[List[float]] = field(default_factory=list)
    trans_prob: list[float] = field(default_factory=list)
    msa_indices: list[int] = field(default_factory=list)
    consensus: str = ""

    @staticmethod
    def parse(buf: List[str]) -> 'HmmRecord':
        assert (buf[0].startswith("HMMER3/f"))

        record = HmmRecord()

        for (i, line) in enumerate(buf[1:]):
            tag, value = line.split(maxsplit=1)
            if tag == "NAME":
                record.header.name = value
            elif tag == "LENG":
                record.header.length = int(value)
            elif tag == "ALPH":
                record.header.alphabet = Alphabet(value)
            elif tag == "NSEQ":
                record.header.n_seq = int(value)
            elif tag == "HMM":
                break

        assert (record.header.complete())

        start = i + 6
        for i in range(start, len(buf), 3):
            mat, ins, trans = [line.split() for line in buf[i:i + 3]]
            record.mat_prob.append([math.exp(-float(v)) for v in mat[1:21]])
            record.ins_prob.append([math.exp(-float(v)) for v in ins])
            record.trans_prob.append([0.0 if v == "*" else math.exp(-float(v)) for v in trans])
            record.consensus += mat[22]
            record.msa_indices.append(int(mat[21]))

        return record


@dataclass
class Hmm:
    records: Dict[str, HmmRecord] = field(default_factory=dict)

    @staticmethod
    def parse(path: str) -> 'Hmm':
        records = {}

        buf = []
        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if line == "//":
                    rec = HmmRecord.parse(buf)
                    records[rec.header.name] = rec
                    buf = []
                else:
                    buf.append(line)

        return Hmm(records)

    def __iter__(self):
        return iter(self.records.values())
