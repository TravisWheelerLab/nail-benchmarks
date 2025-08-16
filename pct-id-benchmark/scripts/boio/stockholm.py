from .fasta import FastaRecord

from dataclasses import dataclass, field, fields
from typing import Dict


@dataclass
class StockholmMeta:
    id: str = ""
    ac: str = ""
    de: str = ""
    au: str = ""
    ga: str = ""
    nc: str = ""
    tc: str = ""


@dataclass
class SequenceMeta:
    name: str
    ac: str


@dataclass
class SequenceEntry:
    name: str
    sequence: str

    def merge(self, other: "SequenceEntry"):
        self.sequence += other.sequence

    def flat_seq(self) -> str:
        return self.sequence.replace(".", "").replace("-", "")

    def fasta_record(self) -> FastaRecord:
        return FastaRecord(self.name, "", self.flat_seq())


@dataclass
class StockholmRecord:
    meta: StockholmMeta = field(default_factory=StockholmMeta)
    sequence_meta: Dict[str, SequenceMeta] = field(default_factory=dict)
    sequences: Dict[str, SequenceEntry] = field(default_factory=dict)

    def __str__(self) -> str:
        lines = ["# STOCKHOLM 1.0"]

        for key in [f.name for f in fields(StockholmMeta)]:
            val = getattr(self.meta, key)
            if val:
                lines.append(f"#=GF {key.upper()} {val}")
        lines.append("")

        for meta in self.sequence_meta.values():
            lines.append(f"#=GS {meta.name} AC {meta.ac}")
        lines.append("")

        for seq in self.sequences.values():
            lines.append(f"{seq.name:<25} {seq.sequence}")

        lines.append("//\n")
        return "\n".join(lines)

    def valid(self) -> bool:
        it = iter(self.sequences.values())
        first_len = len(next(it).sequence)
        return all(len(v.sequence) == first_len for v in it)

    def msa_len(self) -> int:
        it = iter(self.sequences.values())
        return len(next(it).sequence)

    def col(self, idx) -> [str]:
        return [seq.sequence[idx - 1] for seq in self.sequences.values()]


@dataclass
class Stockholm:
    records: Dict[str, StockholmRecord]

    @staticmethod
    def from_path(path: str) -> "Stockholm":
        records = {}
        record = StockholmRecord()

        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("# STOCKHOLM"):
                    continue

                if line.startswith("#=GF"):
                    _, tag, *rest = line.split()
                    value = " ".join(rest)
                    setattr(record.meta, tag.lower(), value)

                elif line.startswith("#=GS"):
                    _, name, tag, *rest = line.split()
                    if tag == "AC":
                        ac = " ".join(rest)
                        record.sequence_meta[name] = SequenceMeta(name, ac)

                elif line == "//":
                    records[record.meta.id] = record
                    record = StockholmRecord()

                elif not line.startswith("#"):
                    name, seq = line.split(maxsplit=1)
                    seq = SequenceEntry(name, seq.strip())
                    if name not in record.sequences:
                        record.sequences[name] = seq
                    else:
                        record.sequences[name].merge(seq)

        return Stockholm(records)

    def __iter__(self):
        return iter(self.records.values())

    def write(self, file):
        for rec in self.records.values():
            file.write(str(rec))
