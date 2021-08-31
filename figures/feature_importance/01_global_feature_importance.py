import sys
import numpy as np
import pandas as pd 
import pytorch_tabnet
import pickle
sys.path.append('/bioinformatics/Users/mclaurt/projects/deep_learning_tumor-only_variant_calling/')

from pytorch_tabnet.metrics import Metric
from sklearn.metrics import average_precision_score

sys.path.append('/bioinformatics/Users/mclaurt/projects/deep_learning_tumor-only_variant_calling/test/no_ps_fit_with_snp-vaf-bin/tcga_pancancer')
from train_features import vaf_bin_mode_features

import matplotlib.pyplot as plt


       
class APS(Metric):
    """          
    Average Precision Score.  A more reliable 'area under the precision recall curve'.
    """          
    def __init__(self):
        self._name = "aps" 
        self._maximize = True
                 
    def __call__(self, y_true, y_score):
        aps = average_precision_score(y_true, y_score[:, 1])
        return aps 


pkl_classifier = '/bioinformatics/Users/mclaurt/projects/deep_learning_tumor-only_variant_calling/test/no_ps_fit_with_snp-vaf-bin/tcga_pancancer/tabnet_trained_tcga_v3.pkl'

features = vaf_bin_mode_features               

clf = pickle.load(open(pkl_classifier,'rb'))

imp = clf.feature_importances_

print(imp)

indices = np.argsort(imp)
print(indices)
print(np.array(features)[indices])

for x in indices:
    print(np.array(features)[x], imp[x])

df = pd.DataFrame({'feature' : np.array(features), 'importance' : imp})
df.to_csv('feature_importance.csv')

# Plot the feature importances of the forest
plt.figure(figsize = (11,12))
plt.title("TabNet global feature importances")
plt.barh(range(len(features)), imp[indices],
       color="r", align="center")

# change indices to a list of labels on the following line.
plt.yticks(range(len(features)), np.array(features)[indices])
plt.ylim([-1, len(features)])
plt.show()
plt.savefig('tabnet_global_feature_importances.png')
