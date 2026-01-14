import argparse
import glob
import numpy as np
import sys
import os
import pyissm

issm_dir = pyissm.tools.config.get_issm_dir()
archive_path = issm_dir + '/test/Archives/'
test_path = issm_dir + './test/'


def id_to_name(test_id):
    """
    Return the name of the test by reading the first line of the test file.
    """

    # Open file and read first line
    test_file = os.path.join(test_path, f'{test_id}.py')
    with open(test_file, 'r') as f:
        first_line = f.readline().strip()
    
    # Remove "#Test Name" prefix if it exists and return test name
    prefix = '#Test Name:'
    if first_line.startswith(prefix):
        return first_line[len(prefix):].strip()
    # Otherwise, just return test number
    else:
        return f'Test_{test_id}'


def run_test(test_id,
             benchmark = 'nightly',
             procedure = 'check',
             output = None):
    """
    Run a single test and return a list of errors
    """

    # Check input validity
    if benchmark != 'nightly':
        raise NotImplementedError(f"'{benchmark}' benchmark is not implementd yet.")
    
    if procedure != 'check':
        raise NotImplementedError(f"'{procedure}' procedure is not implementd yet.")
    
    if output is not None:
        raise NotImplementedError(f"'{output}' output is not implementd yet.")
    
    # Initialise fields
    print(f'--- RUNNING TEST {test_id} ---\n')
    errors = []
    archive_name = f'Archive{test_id}'
    test_name = id_to_name(test_id)

    # Execute test
    test_file = os.path.join(test_path, f'test{test_id}.py')
    try:
        exec(compile(open(test_file, 'rb').read(), test_file, 'exec'), globals())
    except Exception as e:
        print(f'ERROR executing test {test_id}: {e}')
        errors.append({'field': 'EXECUTION_FAILED', 'error_diff': None, 'tolerance': None, 'test_name': test_name})
        return errors
    
    ## Define archive file
    archive_file = os.path.join(archive_path, f'{archive_name}.arch')

    ## Loop over all fields and compute error
    for k, fieldname in enumerate(field_names):
        try:
            field = np.array(field_values[k])
            if field.ndim == 1:
                field = field.reshape((-1, 1)) if field.size else field.reshape((0, 0))
            tolerance = field_tolerances[k]

            archive = np.array(pyissm.tools.archive.arch_read(archive_file, f'{archive_name}_field{k + 1}'))
            if archive is None:
                raise ValueError(f"Field '{archive_name}_field{k + 1}' does not exist in archive file.")

            ### reshape field if necessary
            if field.shape != archive.shape and field.shape not in [(1, 1), (0, 0), (1, 0), (0, 1)]:
                field = field.T
                if field.shape != archive.shape:
                    raise RuntimeError(
                        f"Field '{fieldname}' from test {archive_name[7:]} malformed; "
                        f"shape is {field.shape}, expected {archive.shape} or {archive.T.shape}"
                    )

            ### Compute error
            error_diff = np.amax(np.abs(archive - field), axis=0) / (np.amax(np.abs(archive), axis=0) + float_info.epsilon)
            if not np.isscalar(error_diff):
                error_diff = error_diff[0]

            ### If error, append to errors list
            if np.any(error_diff > tolerance) or np.isnan(error_diff):
                print(f"ERROR   difference: {error_diff:7.2g} > {tolerance:7.2g} "
                      f"test id: {test_id} test name: {test_name} field: {fieldname}")
                errors.append({
                    'field': fieldname,
                    'error_diff': error_diff,
                    'tolerance': tolerance,
                    'test_name': test_name,
                    'test_id': test_id
                })
            ### If no error, print success message
            else:
                print(f"SUCCESS difference: {error_diff:7.2g} < {tolerance:7.2g} "
                      f"test id: {test_id} test name: {test_name} field: {fieldname}")
        
        ### If error, append to errors list
        except Exception as e:
            print(f"Error while checking field '{fieldname}' from test {archive_name[7:]}: {e}")
            errors.append({
                'field': fieldname,
                'error_diff': None,
                'tolerance': None,
                'test_name': test_name,
                'test_id': test_id
            })

    return errors

def run_tests(test_range: str = None, exclude: str = None):
    """Run multiple tests and return a list of all errors."""
    
    # Parse range of tests
    if test_range:
        start, end = map(int, test_range.split(':'))
        tests = list(range(start, end + 1))
    else:
        # If no range is specified, auto-detect all test files
        tests = [int(f.split('test')[1].split('.py')[0]) for f in glob.glob(os.path.join(test_path, 'test*.py'))]

    # Parse exclusions
    exclude_list = [int(x) for x in exclude.split(',')] if exclude else []
    tests_to_run = [t for t in tests if t not in exclude_list]

    # Run all tests and compile list of all errors
    all_errors = []
    for t in tests_to_run:
        all_errors.extend(run_test(t))

    # If any errors exist, print summary of errors
    if all_errors:
        print("\n=== SUMMARY OF FAILED TESTS ===")
        for e in all_errors:
            print(f"Test ID: {e.get('test_id', '?')}, Test Name: {e.get('test_name')}, "
                  f"Field: {e.get('field')}, "
                  f"Error Diff: {e.get('error_diff')}, Tolerance: {e.get('tolerance')}")
        print(f"\nTOTAL FAILED FIELDS: {len(all_errors)}")
    else:
        print("\nALL TESTS PASSED!")

    return len(all_errors) # If any errors occurred, return non-zero exit

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ISSM tests.")
    parser.add_argument("--range", type=str, help="Range of tests to run, e.g., 100:201")
    parser.add_argument("--exclude", type=str, help="Comma-separated list of tests to exclude, e.g., 107,201")
    args = parser.parse_args()

    num_errors = run_tests(test_range=args.range, exclude=args.exclude)
    exit(num_errors)  # non-zero exit code if any errors for CI