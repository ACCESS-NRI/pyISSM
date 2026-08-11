"""One-off: generate src/pyissm/assets/love_numbers.csv from the upstream
ISSM Love-number table.

Run once, commit the resulting CSV. Not part of the pyissm package.

Safety:
- Fetches from a pinned commit SHA, never a mutable branch, so the content
  cannot change between review and execution.
- Verifies the fetched file's SHA-256 against a recorded-at-review-time
  value before parsing anything.
- Extracts only the numeric array literal and parses it with
  `ast.literal_eval` (numeric/list literals only) instead of `eval()` -
  arbitrary code in the fetched text cannot execute.
"""
import ast
import hashlib
import re
import urllib.request

import numpy as np
import pandas as pd

# Pinned to the only commit that has ever touched this file (verified via
# GitHub API history for src/m/boundaryconditions/getlovenumbers.py on
# access-development). Re-pin deliberately if upstream ever changes it.
COMMIT = '6e16985dce80db4c564a9951f6e8d0a51308db0e'
URL = (f'https://raw.githubusercontent.com/ACCESS-NRI/ISSM/{COMMIT}/'
       'src/m/boundaryconditions/getlovenumbers.py')
EXPECTED_SHA256 = (
    '5923c3321e7f42d84892d5ceb169b4165a9171704aeb8ee29de909d83602155d')

COLUMNS = ['h', 'k', 'l', 'th', 'tk', 'tl']
OUT_PATH = 'src/pyissm/assets/love_numbers.csv'


def fetch_verified_source():
    with urllib.request.urlopen(URL) as resp:
        raw = resp.read()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != EXPECTED_SHA256:
        raise RuntimeError(
            f'Checksum mismatch for pinned commit {COMMIT}: '
            f'expected {EXPECTED_SHA256}, got {digest}. Upstream content '
            'at this commit should never change - do not proceed without '
            're-reviewing the diff.')
    return raw.decode()


def extract_love_numbers_array(source):
    m = re.search(r'love_numbers = np\.array\((\[.*?\])\)', source, re.S)
    if m is None:
        raise RuntimeError('Could not locate the love_numbers array literal '
                            'in the fetched source.')
    # ast.literal_eval only accepts literal numbers/strings/lists/tuples/
    # dicts/sets/booleans/None - no names, no calls, no attribute access -
    # so it cannot execute arbitrary code, unlike eval().
    nested = ast.literal_eval(m.group(1))
    return np.array(nested, dtype=float)


def main():
    source = fetch_verified_source()
    arr = extract_love_numbers_array(source)

    assert arr.shape == (10001, 6), f'unexpected shape {arr.shape}'
    assert np.allclose(arr[0], 0.0), 'row 0 (degree 0) should be all zeros'
    assert np.isclose(arr[1, 0], -1.28740059, atol=1e-8), arr[1, 0]
    assert np.isclose(arr[10000, 0], -6.27342778, atol=1e-8), arr[10000, 0]

    pd.DataFrame(arr, columns=COLUMNS).to_csv(OUT_PATH, index=False)
    print(f'Wrote {OUT_PATH}: shape={arr.shape}')


if __name__ == '__main__':
    main()
