#!/usr/bin/env bash
set -euo pipefail

A="${1:?usage: $0 fileA fileB fileC}"
B="${2:?usage: $0 fileA fileB fileC}"
C="${3:?usage: $0 fileA fileB fileC}"

awk -v B="$B" -v C="$C" '
BEGIN { OFS="," }
{
  if ((getline b < B) <= 0) exit
  if ((getline c < C) <= 0) exit

  split(b, fb)
  split(c, fc)

  print fb[3] * fc[3], $10
}
END { close(B); close(C) }
' "$A"
