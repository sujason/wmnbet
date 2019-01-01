"""
Paths, names, optimized hyper-parameters, and other constants for THOMAS.
"""
import os
import json

image_name = 'WMnMPRAGE_bias_corr.nii.gz'
# Find path for priors
this_path = os.path.dirname(os.path.realpath(__file__))
config_json = os.path.join(this_path, 'brain.json')
assert os.path.exists(config_json)
with open(config_json, 'r') as f:
    config = json.load(f)
prior_path = os.path.join(this_path, config['path'])
assert os.path.exists(prior_path)
template = os.path.join(prior_path, config['template'])
assert os.path.exists(template)
subjects = [el for el in os.listdir(prior_path) if os.path.isdir(os.path.join(prior_path, el)) and not el.startswith('.')]
assert len(subjects) > 0


if __name__ == '__main__':
    print subjects
