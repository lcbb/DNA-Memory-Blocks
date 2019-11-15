import csv

import numpy as np

def parse_csv_table(path, ignore_blanks=True, delimiter=','):
  # Looks for the header line beginning with "Well," to determine column headers
  # Subsequent lines are added to the data dict.
  inpath = open(path)
  csvreader = csv.reader(inpath, delimiter=delimiter)

  # skip any extraneous lines at top of file
  # waits for a line starting with the first CSV column equal to "Well" (case-sensitive)
  for line in csvreader:
    if len(line) == 0:  continue

    if line[0] == 'Well':
      column_headers = [v.upper() for v in line[1:]]
      break

  # read each line and add data to the growing dictionary
  data_dict = {}
  for line in csvreader:
    well = line[0]
    vals = line[1:]
    if ignore_blanks:
      well_data = {header: val for header,val in zip(column_headers, vals) if val!=''}
      if len(well_data.keys()) > 1: 
        data_dict[well] = well_data
    else:
      data_dict[well] = {header: val for header,val in zip(column_headers, vals)}

  inpath.close()

  return data_dict
     
def get_wells_by_setup_info(setup_info, well_position = None, sample = None, target = None):
  wells = setup_info.keys()
  if well_position is not None:
    wells = filter(lambda w: setup_info[w]['WELL POSITION'] == well_position, wells) 
  if sample is not None:
    wells = filter(lambda w: setup_info[w]['SAMPLE NAME'] == sample, wells)
  if target is not None:
    wells = filter(lambda w: setup_info[w]['TARGET NAME'] == target, wells)

  return list(wells)


