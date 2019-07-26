## NOTE: 
##   Edit BLANK_NAMES if you want to change the sample names corresponding to blank wells 
##        (default is 'BLANK')
##   Edit PULL_NAMES if you want to set the names of the pulls used in this plate
##        (default is to automatically determine the names of the pulls)

import sys

import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

import qpcr_utils

def preprocess_data(raw_data):
  for well in raw_data:
    for field in raw_data[well]:
      if field in ['CT', 'Tm1', 'Tm2', 'Tm3']:
        val = raw_data[well][field]
        if val == 'Undetermined':
          print "WARNING: Invalid {} value '{}' encountered in well {} ({})".format(field, val, well, raw_data[well]['WELL POSITION'])
          raw_data[well][field] = np.nan
        else:
          raw_data[well][field] = float(raw_data[well][field])

def get_standard_curve(target):
  ## Returns the slope (m) and y-intercept (b) of the fit line:
  ##   CT = m*log_2(amt) + b
  ## Three vals are returned:
  ##   a function mapping CT to amt
  ##   m
  ##   b
  wells_100fg = qpcr_utils.get_wells_by_setup_info(raw_data, sample = '100fg', target = target)
  wells_1pg = qpcr_utils.get_wells_by_setup_info(raw_data, sample = '1pg', target = target)
  wells_10pg = qpcr_utils.get_wells_by_setup_info(raw_data, sample = '10pg', target = target)

  # check and remove wells if they don't show up in raw_data. This shouldn't really happen
  wells_100fg = list(filter(lambda w: w in raw_data, wells_100fg))
  wells_1pg = list(filter(lambda w: w in raw_data, wells_1pg))
  wells_10pg = list(filter(lambda w: w in raw_data, wells_10pg))

  if len(wells_100fg) + len(wells_1pg) + len(wells_10pg) == 0:
    print "WARNING: No qPCR data found for target", target
    return lambda ct: np.nan

  amts = np.array([1e-13]*len(wells_100fg) + [1e-12]*len(wells_1pg) + [1e-11]*len(wells_10pg))
  CT = np.array([raw_data[w]['CT'] for w in wells_100fg+wells_1pg+wells_10pg])

  A = np.concatenate((np.log2(amts).reshape((len(amts),1)), np.ones((len(amts),1))), axis=1)
  B = CT
  res, _, _, _ = np.linalg.lstsq(A, B)

  if VERBOSITY >= 1:
    plt.figure()
    plt.plot(amts, CT, 'ko')
    plt.plot(amts, np.log2(amts)*res[0] + res[1], 'k--')
    plt.xscale('log')
    plt.xlabel('Amount (g)')
    plt.ylabel('CT')
    plt.legend([target, 'fit: CT = {:.2f}*$log_2$(amt) + {:.2f}'.format(res[0],res[1])])
    plt.title('standard curve ({})'.format(target))

    plt.savefig('standard-curve_{}.pdf'.format(target), bbox_inches='tight')
    plt.close()

  return lambda ct: 2**((ct-res[1])/res[0])

def analyze_pull(pull_name, blank_names, targets, correct_targets=None, standard_curves = {}):
  pull_amounts = {}
  pull_amounts_norm = {}
  for target in targets:
    blank_wells = [
        w
        for b in blank_names 
        for w in qpcr_utils.get_wells_by_setup_info(raw_data, sample=b, target=target)
        if w in raw_data
    ]

    if len(blank_wells) == 0:
      print "Warning: No blank wells found for target {}. Is BLANK_NAMES set correctly?".format(target)

    pull_wells = list(filter(
        lambda w: w in raw_data,
        qpcr_utils.get_wells_by_setup_info(raw_data, sample=pull_name, target=target)
    ))

    if target not in standard_curves:
      standard_curves[target] = get_standard_curve(target)
    standard_curve = standard_curves[target]

    blank_amount = np.mean([standard_curve(raw_data[w]['CT']) for w in blank_wells])
    pull_amount = np.mean([standard_curve(raw_data[w]['CT']) for w in pull_wells])

    pull_amounts[target] = pull_amount
    pull_amounts_norm[target] = pull_amount / blank_amount

  return pull_amounts, pull_amounts_norm

VERBOSITY = 1
BLANK_NAMES = ['Blank2']
PULL_NAMES = None # None = automatically determine pull names; replace with list of strings to explicitly list out the pulls

#if len(sys.argv) != 3:
#  print "Usage: python analyze_pulldown_qpcr.py <results_path.csv> <setup_path.csv>"
#  sys.exit(0)
#
#data_path, setup_path = sys.argv[1:3]

if len(sys.argv) != 2:
  print "Usage: python analyze_pulldown_qpcr.py <results_path.csv>"
  sys.exit(0)

data_path = sys.argv[1]


#setup_info = qpcr_utils.parse_csv_table(setup_path)
raw_data = qpcr_utils.parse_csv_table(data_path)
preprocess_data(raw_data)

targets = sorted(set(raw_data[w]['TARGET NAME'] for w in raw_data if 'TARGET NAME' in raw_data[w]))

if PULL_NAMES is None:
  pulls = sorted(set(raw_data[w]['SAMPLE NAME'] for w in raw_data if 'SAMPLE NAME' in raw_data[w]) - set(BLANK_NAMES + ['100fg', '1pg', '10pg']))
  print "Automatically determined pulls on this plate:", pulls
else:
  pulls = PULL_NAMES

pull_data = {}
pull_data_norm = {}
for pull in pulls:
  print "Processing pull: {}...".format(pull),
  pull_data[pull], pull_data_norm[pull] = analyze_pull(pull_name = pull, blank_names = BLANK_NAMES, targets=targets)
  print "Done!"

pull_matrix = np.array([[pull_data[pull][target] for pull in pulls] for target in targets]) * 10**12
pull_matrix_norm = np.array([[pull_data_norm[pull][target] for pull in pulls] for target in targets])
percentage_matrix = 100*np.array([[pull_data[pull][target] / np.nansum(pull_data[pull].values()) for pull in pulls] for target in targets])

plt.figure(figsize=(12,4.8))

plt.subplot(1,3,1)
plt.imshow(pull_matrix, interpolation='nearest', norm=matplotlib.colors.LogNorm(vmin=np.nanmin(pull_matrix), vmax=np.nanmax(pull_matrix)))
plt.tick_params(bottom=False,labeltop=True,labelbottom=False,left=False)
plt.xticks(np.arange(len(pulls)), pulls, rotation='vertical')
plt.yticks(np.arange(len(targets)), targets)
cbar = plt.colorbar()
plt.xlabel('Amount (raw), pg')

plt.subplot(1,3,2)
plt.imshow(pull_matrix_norm, interpolation='nearest', norm=matplotlib.colors.LogNorm(vmin=1.0, vmax=np.nanmax(pull_matrix_norm)))
plt.tick_params(bottom=False,labeltop=True,labelbottom=False,left=False)
plt.xticks(np.arange(len(pulls)), pulls, rotation='vertical')
plt.yticks(np.arange(len(targets)), targets)
plt.colorbar()
plt.xlabel('Enrichment Over Blank')

plt.subplot(1,3,3)
plt.imshow(percentage_matrix, interpolation='nearest', vmin=0.0, vmax=100.0)
plt.tick_params(bottom=False,labeltop=True,labelbottom=False,left=False)
plt.xticks(np.arange(len(pulls)), pulls, rotation='vertical')
plt.yticks(np.arange(len(targets)), targets)
plt.colorbar()
plt.xlabel('Percentage of Sample')

plt.savefig('pull_analysis.pdf', bbox_inches='tight')
