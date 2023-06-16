import numpy as np
import sys, os
import argparse
import gc
import torch
import time
import h5py
import copy
from collections import OrderedDict
import torch.nn as nn
import torchvision.models as models

# install from https://pypi.org/project/simclr/
# https://github.com/Spijkervet/SimCLR
import simclr
from simclr import SimCLR
# from simclr.modules import get_resnet

# clip implemented in this package, from:
# https://github.com/openai/CLIP
import clip

# Define stuff about the resnet layers here
n_features_each_resnet_block = [256,256,256, 512,512,512,512, 1024,1024,1024,1024,1024,1024, 2048,2048,2048]
resnet_block_names = ['block%d'%nn for nn in range(len(n_features_each_resnet_block))]


# only need this because it specifies where simclr pretrained weights are stored
path = os.path.realpath(__file__).split('/')
path = '/'+os.path.join(*path[0:-1])
print(path)
simclr_weights_path = os.path.join(path, 'features','simclr')

def init_cuda():
    
    # Set up CUDA stuff
    
    print ('#device:', torch.cuda.device_count())
    print ('device#:', torch.cuda.current_device())
    print ('device name:', torch.cuda.get_device_name(torch.cuda.current_device()))

    torch.manual_seed(time.time())
    device = torch.device("cuda:0") #cuda
    torch.backends.cudnn.enabled=True

    print ('\ntorch:', torch.__version__)
    print ('cuda: ', torch.version.cuda)
    print ('cudnn:', torch.backends.cudnn.version())
    print ('dtype:', torch.get_default_dtype())

    return device

try:
    device = init_cuda()
except:
    device = 'cpu:0'
print('device = %s'%device)


def extract_features(image_data, \
                    block_inds, \
                    pooling_op = None,
                    save_dtype=np.float32,\
                    training_type='imgnet', \
                    batch_size=100, \
                    debug=False):
    """ 
    Extract features for resnet model.
    Image data is [n_images x 3 x n_pix x n_pix]
    Return features [n_images x n_channels]
    
    """

    debug = debug==1
    print('debug=%s'%debug)
    
    assert(len(block_inds)==1)
    ll = block_inds[0]
    
    n_images = image_data.shape[0]
    
    model_architecture='RN50'
    
    do_tile=True
    model_filename=None
   
    n_batches = int(np.ceil(n_images/batch_size))

    # figure out how big features will be, by passing a test image through
    if do_tile:
        image_template = np.tile(image_data[0:1,:,:,:], [1,3,1,1])
    else:
        image_template = image_data[0:1,:,:,:]
    activ_template = get_resnet_activations_batch(image_template, block_inds, \
                                                 model_architecture, training_type, \
                                                  model_filename=model_filename, device=device)
    if pooling_op is not None:
        activ_template = [pooling_op(activ_template[0])]
    n_features_total = np.prod(activ_template[0].shape[1:])  
    print('number of features total: %d'%n_features_total)
    
    features = np.zeros((n_images, n_features_total),dtype=save_dtype)

    with torch.no_grad():

        
        for bb in range(n_batches):

            if debug and bb>1:
                continue
            print('Processing images for batch %d of %d'%(bb, n_batches))
            sys.stdout.flush()
            
            batch_inds = np.arange(batch_size * bb, np.min([batch_size * (bb+1), n_images]))

            # using grayscale images for better comparison w my other models.
            # need to tile to 3 so model weights will be right size
            if do_tile:
                image_batch = np.tile(image_data[batch_inds,:,:,:], [1,3,1,1])
            else:
                image_batch = image_data[batch_inds,:,:,:]

            gc.collect()
            torch.cuda.empty_cache()

            activ_batch = get_resnet_activations_batch(image_batch, block_inds, \
                                                 model_architecture, training_type, \
                                                 model_filename=model_filename, device=device)
            if bb==0:
                print('size of activ this batch raw:')
                print(activ_batch[0].shape)
              
            if pooling_op is not None:
                
                activ_batch = [pooling_op(activ_batch[0])]
                
                if bb==0:
                    print('size of activ pooled:')
                    print(activ_batch[0].shape)
              
            activ_batch_reshaped = torch.reshape(activ_batch[0], [len(batch_inds), -1])
           
            if bb==0:
                print('size of activ reshaped:')
                print(activ_batch_reshaped.shape)
                
            features[batch_inds,:] = activ_batch_reshaped.detach().cpu().numpy()
        
    return features

        
def get_resnet_activations_batch(image_batch, \
                               block_inds, \
                               model_architecture, \
                               training_type, \
                               model_filename = None, \
                               device=None):

    """
    Get activations for images passed through pretrained resnet model.
    Specify which which layers to return.
    """

    if device is None:
        device = torch.device('cpu:0')
       
    if training_type=='clip':        
        
        print('Using CLIP model')
        model, preprocess = clip.load(model_architecture, device=device)
        model = model.visual
        
    elif training_type=='imgnet':
        
        # normal pre-trained model from pytorch
        print('Using pretrained Resnet50 model')
        model = models.resnet50(weights=models.ResNet50_Weights.IMAGENET1K_V2).float().to(device)
        
    elif training_type=='simclr':
        
        print('using SimCLR model')
        encoder = simclr.modules.get_resnet('resnet50',pretrained=False)
        projection_dim = 64
        n_features = encoder.fc.in_features  # get dimensions of last fully-connected layer
        model = SimCLR(encoder, projection_dim, n_features)
        model_fp = os.path.join(simclr_weights_path, 'checkpoint_100.tar')
        print('checkpoint path: %s'%model_fp)
        if device=='cpu:0':
            device_use = torch.device('cpu')
        else:
            device_use = device
        print(device_use)
        model.load_state_dict(torch.load(model_fp, map_location=device_use))
        model = model.encoder.to(device_use)
        
    else:
        raise ValueError('training type %s not recognized'%training_type)
        
    model.eval()
    
    # The 16 residual blocks are segmented into 4 groups here, which have different numbers of features.
    blocks_each= [len(model.layer1), len(model.layer2), len(model.layer3),len(model.layer4)]
    which_group = np.repeat(np.arange(4), blocks_each)

    activ = [[] for ll in block_inds]
    hooks = [[] for ll in block_inds]
    
    # first making this subfunction that is needed to get the activation on a forward pass
    def get_activ_fwd_hook(ii,ll):
        def hook(self, input, output):
            # the relu operation is used multiple times per block, but we only 
            # want to save its output when it has this specific size.
            if output.shape[1]==n_features_each_resnet_block[ll]:
                print('executing hook for %s'%resnet_block_names[ll])  
                activ[ii] = output
                print(output.shape)
        return hook

    image_tensors = torch.Tensor(image_batch).to(device)
    with torch.no_grad():

        # adding a "hook" to the module corresponding to each layer, so we'll save activations at each layer.
        # For resnet, going to save output of each residual block following last relu operation.
        for ii, ll in enumerate(block_inds):
            if which_group[ll]==0:            
                h = model.layer1[ll].relu.register_forward_hook(get_activ_fwd_hook(ii,ll))
            elif which_group[ll]==1:            
                h = model.layer2[ll-blocks_each[0]].relu.register_forward_hook(get_activ_fwd_hook(ii,ll))
            elif which_group[ll]==2:            
                h = model.layer3[ll-sum(blocks_each[0:2])].relu.register_forward_hook(get_activ_fwd_hook(ii,ll))
            elif which_group[ll]==3:            
                h = model.layer4[ll-sum(blocks_each[0:3])].relu.register_forward_hook(get_activ_fwd_hook(ii,ll))
            else:
                h=None
            hooks[ii] = h

        # Pass images though the model (hooks get run now)
        image_features = model(image_tensors)

        # Now remove all the hooks
        for ii, ll in enumerate(block_inds):
            hooks[ii].remove

    # Sanity check that we grabbed the right activations - check their sizes against expected
    # output size of each block
    exp_size = np.array(n_features_each_resnet_block)[block_inds]
    actual_size = [activ[bb].shape[1] for bb in range(len(activ))]
    assert(np.all(np.array(actual_size)==np.array(exp_size)))

    return activ
