import cPickle as cp
import scipy.io

# Fields inside structure
#       xlabel: [502x1 double]
#            X: [502x16 double]
#            C: 'cpncnnccnncpnnpc'
#   cat_labels: {[]  []  []  {10x1 cell}  []  [] {4x1 cell}  ... []}
#       ylabel: {'stage'  'rx'  'dtime' 'status'  'age'  'wt'  'pf'... 'bm'}
#  ylabel_long: {1x16 cell}
input_file = '../../datasets/mat/prostate_v3.mat'
tmp = scipy.io.loadmat(input_file)
tmp = tmp['data'][0,0]

data = dict() # data is a dictionary with the following keys
data['X'] = tmp[1]
data['C'] = tmp[2]
data['C'] = str(data['C'][0])
data['ylabel'] = tmp[5]
data['cat_labels'] = tmp[3]

idx_toKeep = [0, 1, 3, 12, 14]
data['X'] = data['X'][:,idx_toKeep]
data['C'] = ''.join([data['C'][i] for i in idx_toKeep])
data['ylabel'] = ['Stage','Drug level', 'Prognostic Status',\
        'Size of Primary Tumor (cm^2)', 'Serum Prostatic Acid Phosphatase']
data['cat_labels'] = [['3','4'],['0','1','5'],['alive','vascular',\
        'prostatic ca','others'],None,None]

with open('../../datasets/py/prostate.pk','wb') as f:
    cp.dump(data,f)
