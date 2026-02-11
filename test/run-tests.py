import argparse
import glob
import numpy as np
from sys import float_info
import os
import re
import pyissm

# Define ISSM archive directory & local test directory
issm_dir = pyissm.tools.config.get_issm_dir()
archive_path = issm_dir + '/test/Archives/'
test_path = './'

# -------------------------
# NEW: selector parsing
# -------------------------
def parse_test_selector(selector: str) -> list[int]:
    """
    Parse a comma-separated selector of:
      - single ids: "123"
      - ranges: "100:120" (inclusive)
    Example: "101,105:110,200" -> [101,105,106,107,108,109,110,200]
    """
    if selector is None or str(selector).strip() == "":
        return []

    out: list[int] = []
    parts = [p.strip() for p in selector.split(",") if p.strip()]
    for p in parts:
        if ":" in p:
            a, b = p.split(":", 1)
            a = a.strip()
            b = b.strip()
            if not a.isdigit() or not b.isdigit():
                raise ValueError(f"Bad range '{p}'. Expected like 100:120")
            start = int(a)
            end = int(b)
            if end < start:
                raise ValueError(f"Bad range '{p}'. End < start.")
            out.extend(range(start, end + 1))
        else:
            if not p.isdigit():
                raise ValueError(f"Bad test id '{p}'. Expected an integer.")
            out.append(int(p))

    # dedupe while preserving order
    seen = set()
    deduped = []
    for t in out:
        if t not in seen:
            deduped.append(t)
            seen.add(t)
    return deduped

# Define functions
def id_to_name(test_id):
    """
    Return the name of the test by reading the first line of the test file.
    """
    test_file = os.path.join(test_path, f'test{test_id}.py')
    with open(test_file, 'r') as f:
        first_line = f.readline().strip()

    prefix = '#Test Name:'
    if first_line.startswith(prefix):
        return first_line[len(prefix):].strip()
    else:
        return f'Test_{test_id}'

def _extract_test_id_from_path(test_file: str) -> int:
    """
    Extract test ID from a filename like test123.py (any path).
    """
    base = os.path.basename(test_file)
    m = re.fullmatch(r"test(\d+)\.py", base)
    if not m:
        raise ValueError(f"Expected filename like 'test123.py', got: {base}")
    return int(m.group(1))

def run_test(test_id,
             benchmark='nightly',
             procedure='check',
             output=None):
    """
    Run a single test and return a list of errors
    """
    if benchmark != 'nightly':
        raise NotImplementedError(f"'{benchmark}' benchmark is not implementd yet.")
    if procedure != 'check':
        raise NotImplementedError(f"'{procedure}' procedure is not implementd yet.")
    if output is not None:
        raise NotImplementedError(f"'{output}' output is not implementd yet.")

    print(f"\n{'-'*100}\n{'RUNNING TEST ' + str(test_id):^{100}}\n{'-'*100}\n", flush=True)
    errors = []
    archive_name = f'Archive{test_id}'
    test_name = id_to_name(test_id)

    test_file = os.path.join(test_path, f'test{test_id}.py')
    archive_file = os.path.join(archive_path, f'{archive_name}.arch')

    if not os.path.isfile(test_file):
        print(f"ERROR: Test file does not exist: {test_file}", flush=True)
        errors.append({
            'field': 'TEST_FILE_MISSING',
            'error_diff': None,
            'tolerance': None,
            'test_name': test_name,
            'test_id': test_id
        })
        return errors

    if not os.path.isfile(archive_file):
        print(f"ERROR: Archive file does not exist: {archive_file}", flush=True)
        errors.append({
            'field': 'ARCHIVE_FILE_MISSING',
            'error_diff': None,
            'tolerance': None,
            'test_name': test_name,
            'test_id': test_id
        })
        return errors

    try:
        local_globals = {}
        exec(compile(open(test_file, 'rb').read(), test_file, 'exec'), local_globals, local_globals)
        field_names = local_globals.get('field_names', [])
        field_values = local_globals.get('field_values', [])
        field_tolerances = local_globals.get('field_tolerances', [])
    except Exception as e:
        print(f'ERROR executing test {test_id}: {e}', flush=True)
        errors.append({'field': 'EXECUTION_FAILED', 'error_diff': None, 'tolerance': None,
                       'test_name': test_name, 'test_id': test_id})
        return errors

    for k, fieldname in enumerate(field_names):
        try:
            field = np.array(field_values[k])
            if field.ndim == 1:
                field = field.reshape((-1, 1)) if field.size else field.reshape((0, 0))
            tolerance = field_tolerances[k]

            archive = np.array(pyissm.tools.archive.arch_read(archive_file, f'{archive_name}_field{k + 1}'))
            if archive is None:
                raise ValueError(f"Field '{archive_name}_field{k + 1}' does not exist in archive file.")

            if field.shape != archive.shape and field.shape not in [(1, 1), (0, 0), (1, 0), (0, 1)]:
                field = field.T
                if field.shape != archive.shape:
                    raise RuntimeError(
                        f"Field '{fieldname}' from test {archive_name[7:]} malformed; "
                        f"shape is {field.shape}, expected {archive.shape} or {archive.T.shape}"
                    )

            error_diff = np.amax(np.abs(archive - field), axis=0) / (np.amax(np.abs(archive), axis=0) + float_info.epsilon)
            if not np.isscalar(error_diff):
                error_diff = error_diff[0]

            if np.any(error_diff > tolerance) or np.isnan(error_diff):
                print(f"ERROR   difference: {error_diff:7.2g} > {tolerance:7.2g} "
                      f"test id: {test_id} test name: {test_name} field: {fieldname}", flush=True)
                errors.append({
                    'field': fieldname,
                    'error_diff': error_diff,
                    'tolerance': tolerance,
                    'test_name': test_name,
                    'test_id': test_id
                })
            else:
                print(f"SUCCESS difference: {error_diff:7.2g} < {tolerance:7.2g} "
                      f"test id: {test_id} test name: {test_name} field: {fieldname}", flush=True)

        except Exception as e:
            print(f"Error while checking field '{fieldname}' from test {archive_name[7:]}: {e}", flush=True)
            errors.append({
                'field': fieldname,
                'error_diff': None,
                'tolerance': None,
                'test_name': test_name,
                'test_id': test_id
            })

    return errors

def run_test_file(test_file: str):
    """
    Run exactly one test file (path), expecting it to be named test<ID>.py.
    """
    test_id = _extract_test_id_from_path(test_file)
    global test_path
    test_path = os.path.dirname(os.path.abspath(test_file)) or "."
    return run_test(test_id)

# -------------------------
# UPDATED: run_tests signature + logic
# -------------------------
def run_tests(tests: str = None, exclude: str = None):
    """Run multiple tests and return non-zero count of failed fields."""

    available_tests = [int(f.split('test')[1].split('.py')[0])
                       for f in glob.glob(os.path.join(test_path, 'test*.py'))]

    if tests:
        requested_tests = parse_test_selector(tests)
        tests_to_run = [t for t in requested_tests if t in available_tests]
    else:
        tests_to_run = available_tests

    exclude_list = parse_test_selector(exclude) if exclude else []
    tests_to_run = [t for t in tests_to_run if t not in exclude_list]
    tests_to_run.sort()

    all_errors = []
    for t in tests_to_run:
        all_errors.extend(run_test(t))

    if all_errors:
        print("\n=== SUMMARY OF FAILED TESTS ===")
        for e in all_errors:
            print(f"Test ID: {e.get('test_id', '?')}, Test Name: {e.get('test_name')}, "
                  f"Field: {e.get('field')}, "
                  f"Error Diff: {e.get('error_diff')}, Tolerance: {e.get('tolerance')}", flush=True)
        print(f"\nTOTAL FAILED FIELDS: {len(all_errors)}")
    else:
        print("\nALL TESTS PASSED!")

    return len(all_errors)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ISSM tests.")
    parser.add_argument("--tests", type=str,
                        help="Comma-separated list of tests and/or ranges, e.g. 101,105:110,200. "
                             "If omitted, runs all test*.py in the directory.")
    parser.add_argument("--exclude", type=str,
                        help="Comma-separated list of tests and/or ranges to exclude, e.g. 107,200:210")
    parser.add_argument("--testfile", type=str,
                        help="Run exactly one test file, e.g., ./test123.py or path/to/test123.py")
    args = parser.parse_args()

    if args.testfile:
        errors = run_test_file(args.testfile)
        num_errors = len(errors)
        if num_errors:
            print("\n=== SUMMARY OF FAILED TEST ===")
            for e in errors:
                print(f"Test ID: {e.get('test_id', '?')}, Test Name: {e.get('test_name')}, "
                      f"Field: {e.get('field')}, "
                      f"Error Diff: {e.get('error_diff')}, Tolerance: {e.get('tolerance')}", flush=True)
            print(f"\nTOTAL FAILED FIELDS: {num_errors}")
        else:
            print("\nTEST PASSED!")
    else:
        num_errors = run_tests(tests=args.tests, exclude=args.exclude)

    exit(num_errors)
