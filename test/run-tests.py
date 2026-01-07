# import sys
# sys.path.append('/Users/lawrence.bird/pyISSM/src/')
import pyissm
import numpy as np
from sys import float_info

# import os
# os.chdir('/Users/lawrence.bird/pyISSM/')

# archive_path = '/Users/lawrence.bird/ISSM/test/Archives/'
archive_path = '../../../issm/test/Archives/'

def IdToName(test_id):
    """IdToName - return name of test

    Usage:
        name = IdToName(test_id)
    """
    infile = open('./test' + str(test_id) + '.py', 'r')
    file_text = infile.readline()

    string = '#Test Name:'
    name = file_text[len(string) + 1:-1]
    return name

def run_tests(id = None,
              exclude = None,
              benchmark = 'nightly',
              procedure = 'check',
              output = None):
    
    print(f'-- RUNNING TEST {id} --')
    
    archive_name = 'Archive' + str(id)
    id_string = IdToName(id)
    
    # Execute test to generate results
    exec(compile(open('./test{}.py'.format(id), "rb").read(), './test{}.py'.format(id), 'exec'), globals())

    # Load archive file
    archive_file = archive_path + 'Archive{}.arch'.format(id)

    for k, fieldname in enumerate(field_names):
        try:
            # Get field and tolerance
            field = np.array(field_values[k])
            if len(field.shape) == 1:
                if np.size(field):
                    field = field.reshape((np.size(field), 1))
                else:
                    field = field.reshape((0, 0))
            tolerance = field_tolerances[k]

            # Read data from archive
            archive = np.array(pyissm.tools.archive.arch_read(archive_file, archive_name + '_field' + str(k + 1)))
            if str(archive) == 'None':
                raise NameError("Field name '" + archive_name + '_field' + str(k + 1) + "' does not exist in archive file.")
            if np.shape(field) != np.shape(archive) and not np.shape(field) in [(1, 1), (0, 0), (1, 0), (0, 1)]:
                field = field.T
                if np.shape(field) != np.shape(archive):
                    raise RuntimeError("Field '{}' from test {} is malformed; shape is {}, should be {} or {}".format(fieldname, archive_name[7:], np.shape(field.T), np.shape(archive), np.shape(archive.T)))

            error_diff = np.amax(np.abs(archive - field), axis=0) / (np.amax(np.abs(archive), axis=0) + float_info.epsilon)
            if not np.isscalar(error_diff):
                error_diff = error_diff[0]

                                # Display test result
            if (np.any(error_diff > tolerance) or np.isnan(error_diff)):
                print(('ERROR   difference: {:7.2g} > {:7.2g} test id: {} test name: {} field: {}'.format(error_diff, tolerance, id, id_string, fieldname)))
                error_count += 1
            else:
                print(('SUCCESS difference: {:7.2g} < {:7.2g} test id: {} test name: {} field: {}'.format(error_diff, tolerance, id, id_string, fieldname)))

        except Exception as e:
            print("Error while checking field '{}' from test {}: {}".format(fieldname, archive_name[7:], str(e)))
            continue

run_tests(id = 201)
# run_tests(id = 208)
# run_tests(id = 212)
# run_tests(id = 222)
# run_tests(id = 223)