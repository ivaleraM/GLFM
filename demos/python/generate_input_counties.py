import cPickle as cp
import scipy.io

import pdb

# Fields inside structure
#       xlabel: [502x1 double]
#            X: [502x16 double]
#            C: 'cpncnnccnncpnnpc'
#   cat_labels: {[]  []  []  {10x1 cell}  []  [] {4x1 cell}  ... []}
#       ylabel: {'stage'  'rx'  'dtime' 'status'  'age'  'wt'  'pf'... 'bm'}
#  ylabel_long: {1x16 cell}
input_file = '../../datasets/mat/counties.mat'
tmp = scipy.io.loadmat(input_file)
tmp = tmp['data'][0,0]

data = dict() # data is a dictionary with the following keys
data['X'] = tmp[0]
D = data['X'].shape[1]

data['C'] = tmp[1]
data['C'] = str(data['C'][0])
data['ylabel'] = map(lambda x: str(x.tolist()[0]), tmp[3][0].tolist())

#data['cat_labels'] = map(lambda x: x.tolist(), tmp[2][0].tolist()
data['cat_labels'] = [None] * D
for d in xrange(D):
    if len(tmp[2][0].tolist()[d].tolist()[0]) > 0:
        data['cat_labels'][d] = map(lambda x: str(x[0][0].tolist()), tmp[2][0].tolist()[d].tolist())

pdb.set_trace()

## ADAPT INPUT DATA --> put bias
data['X'][:,8] = data['X'][:,7] + data['X'][:,8]
data['ylabel'][8] = 'age >= 65'

# idx_to_remove = [1,3,4, 6, 7,8, 10,19]; # [2, 5, 9, 11, 12, 13, 14, 15, 16, 17, 18]
idx_toKeep = [1, 4, 8, 10, 11, 12, 13, 14, 15, 16, 17]
data['X'] = data['X'][:,idx_toKeep]
data['C'] = ''.join([data['C'][i] for i in idx_toKeep])
data['cat_labels'] = [data['cat_labels'][i] for i in idx_toKeep]
data['ylabel'] = [data['ylabel'][i] for i in idx_toKeep]

with open('../../datasets/py/counties.pk','wb') as f:
    cp.dump(data,f)
