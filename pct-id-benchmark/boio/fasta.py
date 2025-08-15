from dataclasses import dataclass, field
from typing import Dict


@dataclass
class FastaRecord:
    name: str
    extra: str
    sequence: str

    def __str__(self) -> str:
        lines = [f">{self.name} {self.extra}"]
        lines += [self.sequence[i:i + 80] for i in range(0, len(self.sequence), 80)]
        lines += [""]
        return "\n".join(lines)


@dataclass
class Fasta:
    records: Dict[str, FastaRecord] = field(default_factory=dict)

    @staticmethod
    def from_path(path: str) -> "Fasta":
        fasta = Fasta()
        with open(path, "r") as f:
            header = None
            seq_chunks = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header:
                        parts = header.split(maxsplit=1)
                        name = parts[0]
                        extra = parts[1] if len(parts) > 1 else ""
                        fasta.records[name] = FastaRecord(name, extra, ''.join(seq_chunks))
                    header = line[1:]
                    seq_chunks = []
                else:
                    seq_chunks.append(line)

            if header:
                parts = header.split(maxsplit=1)
                name = parts[0]
                extra = parts[1] if len(parts) > 1 else ""
                fasta.records[name] = FastaRecord(name, extra, ''.join(seq_chunks))

        return fasta

    def write(self, file):
        for rec in self.records.values():
            file.write(str(rec))

    def __iter__(self):
        return iter(self.records.values())
