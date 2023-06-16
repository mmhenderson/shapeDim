import sys, os
import numpy as np
import time
from sklearn import decomposition
import pandas as pd
import torch
import PIL.Image

print(sys.path)

# root directory is 2 dirs up from this file
path = os.path.realpath(__file__).split('/')
root = '/'+os.path.join(*path[0:-3])

import extract_resnet_features
import image_utils

def get_resnet_features(debug=False, training_type='simclr', \
                        image_set_name = 'Images_grid3'):
    
    # image_set_name = 'images_expt1'
    n_comp_keep = 500;
    
    feat_dir = os.path.join(root, 'Analysis', 'image_similarity', 'features')
    image_dir = os.path.join(root, 'Stimuli','AmpGrid3_adj_full_grey_small')

    if not os.path.exists(feat_dir):
        os.makedirs(feat_dir)


    if training_type=='imgnet':
        feat_path = os.path.join(feat_dir, 'resnet')
    elif training_type=='simclr':
        feat_path = os.path.join(feat_dir, 'simclr')
        
    if debug:
        feat_path = os.path.join(feat_path,'DEBUG')

    if not os.path.exists(feat_path):
        os.makedirs(feat_path)

    shape_df = image_utils.make_shape_labels(image_dir)
    fn2save = os.path.join(feat_dir, 'Image_labels_grid3.csv')
    shape_df.to_csv(fn2save)

    image_data = image_utils.load_images(shape_df, debug=debug)

    blocks_to_do = [2,6,12,15]
    
    # loop over resnet blocks
    for ll in blocks_to_do:

        block_inds = [ll]

        # if ll<=2:
        # reduce size of larger feature maps
        pooling_op = torch.nn.MaxPool2d(kernel_size=4, stride=4, padding=0)
        # else:
        #     pooling_op = None

        if debug:
            image_data_use = image_data[0:2,:,:,:]
        else:
            image_data_use = image_data
            
        # first extract features for all pixels/feature channels 
        features_raw = extract_resnet_features.extract_features(image_data_use,\
                                        block_inds,\
                                        pooling_op = pooling_op, 
                                        save_dtype=np.float32,\
                                        training_type=training_type, \
                                        debug=debug)

        # reduce the dimensionality of the activs here
        scores, wts, pre_mean, ev = compute_pca(features_raw)

        n_keep = np.min([scores.shape[1], n_comp_keep])
        
        scores = scores[:,0:n_keep]

        feat_file_name = os.path.join(feat_path, \
                                      '%s_%s_block%d_pca.npy'%(image_set_name,\
                                                               training_type, \
                                                               ll))
        print('size of scores is:')
        print(scores.shape)
        print('saving to %s'%feat_file_name)
        np.save(feat_file_name, scores)

        
def compute_pca(values, max_pc_to_retain=None, copy_data=False):
    """
    Apply PCA to the data, return reduced dim data as well as weights, var explained.
    """
    n_features_actual = values.shape[1]
    n_trials = values.shape[0]
    
    if max_pc_to_retain is not None:        
        n_comp = np.min([np.min([max_pc_to_retain, n_features_actual]), n_trials])
    else:
        n_comp = np.min([n_features_actual, n_trials])
         
    print('Running PCA: original size of array is [%d x %d], dtype=%s'%\
          (n_trials, n_features_actual, values.dtype))
    sys.stdout.flush()
    t = time.time()
    pca = decomposition.PCA(n_components = n_comp, copy=copy_data)
    scores = pca.fit_transform(values)           
    elapsed = time.time() - t
    print('Time elapsed: %.5f'%elapsed)
    values = None            
    wts = pca.components_
    ev = pca.explained_variance_
    ev = ev/np.sum(ev)*100
    pre_mean = pca.mean_
  
    return scores, wts, pre_mean, ev
