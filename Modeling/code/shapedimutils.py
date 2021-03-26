# -*- coding: utf-8 -*-
"""
Spyder Editor

General use functions related to shapedim neural network modeling.

"""

import torch
import numpy as np
import scipy
import torchvision.models as models
from torchvision import transforms
from torch.utils.data import Dataset, DataLoader
from PIL import Image
from skimage import io
import matplotlib.pyplot as plt
from collections import OrderedDict
from sklearn import decomposition
from sklearn import svm
from matplotlib import cm
import torch.nn.functional as F
from torch import nn
import math
import pandas as pd
import os


  
class shape_dataset(Dataset):
  
  ## Creating a custom dataset class that can be used for my shape images

  def __init__(self, csv_file, image_dir, transform=None):

      self.shape_labels = pd.read_csv(csv_file)
      self.image_dir = image_dir
      self.transform = transform

  def __len__(self):

      return len(self.shape_labels)

  def __getitem__(self, idx):

      if torch.is_tensor(idx):
        idx = idx.tolist()

      imfn = os.path.join(self.image_dir, 'Shape_%.2f_%.2f.png'%(self.shape_labels.coord_axis1[idx], self.shape_labels.coord_axis2[idx]))
      im = io.imread(imfn)        
      im = np.tile(np.expand_dims(im,2), [1,1,3]).astype('double')/255
      # tile to add RGB as third dim 
      # note that once the transformation ToTensor is applied, the dimensions 
      # get automatically switched so RGB is first [batch_size x 3 x H x W]
      # So, the code will break if you don't include ToTensor in the transform
      # im = np.moveaxis(np.tile(np.expand_dims(im,2), [1,1,3]),[0,1,2],[1,2,0]).astype('double')
      task_labels = np.array([self.shape_labels['labels_task1'][idx], 
                              self.shape_labels['labels_task2'][idx], 
                              self.shape_labels['labels_task3'][idx]]).astype('int')
     
      if self.transform:
        im = self.transform(im)

      item = {'image':im, 'task_labels':task_labels}

      return item
    
    

def make_shape_labels(image_dir):
  
  ## Create a csv file that labels all the images in dataset
  # only need to run this once to make the file (saving to same folder where images are)
  
  # create image labels
  start=0; stop=5; step=0.1
  grid_x = np.round(np.arange(start,stop+step,step),1)
  grid_y = np.round(np.arange(start,stop+step,step),1)
  x, y = np.meshgrid(grid_x, grid_y)
  all_grid_points = np.column_stack((x.ravel(),y.ravel()))
  
  center = 2.5  # center of shape space is the "boundary"
  all_quadrant = np.zeros([np.shape(all_grid_points)[0],1]);
  all_quadrant[np.logical_and(all_grid_points[:,0]>center, all_grid_points[:,1]>center)] = 1;
  all_quadrant[np.logical_and(all_grid_points[:,0]<center, all_grid_points[:,1]>center)] = 2;
  all_quadrant[np.logical_and(all_grid_points[:,0]<center, all_grid_points[:,1]<center)] = 3;
  all_quadrant[np.logical_and(all_grid_points[:,0]>center, all_grid_points[:,1]<center)] = 4;
  
  labels_task1 = np.zeros(np.shape(all_quadrant))
  labels_task1[np.isin(all_quadrant,[1,4])] = 0
  labels_task1[np.isin(all_quadrant,[2,3])] = 1
  
  labels_task2 = np.zeros(np.shape(all_quadrant))
  labels_task2[np.isin(all_quadrant,[1,2])] = 0
  labels_task2[np.isin(all_quadrant,[3,4])] = 1
  
  labels_task3 = np.zeros(np.shape(all_quadrant))
  labels_task3[np.isin(all_quadrant,[1,3])] = 0
  labels_task3[np.isin(all_quadrant,[2,4])] = 1
  
  shape_labels = pd.DataFrame.from_dict({'coord_axis1':all_grid_points[:,0], 'coord_axis2':all_grid_points[:,1],
     'quadrant':np.squeeze(all_quadrant), 'labels_task1':np.squeeze(labels_task1),
     'labels_task2':np.squeeze(labels_task2), 'labels_task3':np.squeeze(labels_task3)})
  
  savefn = os.path.join(image_dir, 'shape_labels_all.csv')
  print('saving shape labels dataframe to %s'%savefn)
  shape_labels.to_csv(savefn)
  print('saved')
  # shape_labels = pd.read_csv(savefn)
  


def get_alexnet_activations(image_dir, transform=None):
    
    # Get activations for all images in the dataset in batches
    # will return a nested ordered dict containing all layer activations

    # first making this subfunction that is needed to get the activation on a forward pass
    def get_activ_fwd_hook(ii):

      def hook(self, input, output):

          activ[ii] = output.detach().numpy() # detach just means make a copy whose gradients we stop measuring

      return hook

    # first loading pre-trained model from torch model zoo
    model = models.alexnet(pretrained=True).double() 
    model.eval()
   
    # will loop over all the layers in the "features" module (which is a Sequential container of modules)
    nfeaturelayers = 13
    # layer_names = ['1_Conv2d','1_ReLU','1_MaxPool2d',
    #                '2_Conv2d','2_ReLU','2_MaxPool2d',
    #                '3_Conv2d','3_ReLU',
    #                '4_Conv2d','4_ReLU',
    #                '5_Conv2d','5_ReLU','5_MaxPool2d']
    
    # initialize this - will be big nested ordereddict to store all activs
    all_activ = OrderedDict()

    csv_file = os.path.join(image_dir, 'shape_labels_all.csv') 
    sds = shape_dataset(csv_file, image_dir)
    
    if transform==None:
        # create a transformation that will normalize the image intensity (z-score)
        # these values were copied from https://pytorch.org/tutorials/beginner/finetuning_torchvision_models_tutorial.html
        # supposed to match pre-training data for these models
        transform = transforms.Compose([transforms.ToTensor(), 
                                        transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])])
            
    sds = shape_dataset(csv_file, image_dir, transform = transform)

    # create a dataloader object that breaks whole set into smaller batches
    batch_size=100
    # note this is NOT shuffled, so we will know exactly which images are in which batch
    dloader = DataLoader(sds, batch_size = batch_size, shuffle=False)
    dgenerator = iter(dloader)
    nbatches = len(dloader)

    for bb in range(nbatches):
      
        print('getting activs for batch %d of %d'%(bb, nbatches))
        curr_batch = next(dgenerator)

        # get image and labels for this batch
        # image_tensors is [batch_size x 3 x 224 x 224]
        image_tensors =  curr_batch['image']
        hooks = OrderedDict()
        activ = OrderedDict()

        model.eval()
        # adding this "hook" to the module corresponding to each layer, so we'll save activations at each layer
        # this only modifies the "graph" e.g. what the model code does when run, but doesn't actually run it yet.
        for ii in range(nfeaturelayers):  
            hooks[ii] = model.features[ii].register_forward_hook(get_activ_fwd_hook(ii))

        # do the forward pass of model, which now includes the forward hooks
        # now the "activ" variable will get modified, because it gets altered during the hook function
        model(image_tensors)

        all_activ[bb] = activ

        # removing the hooks now, because otherwise they get added multiple times
        for ii in range(nfeaturelayers):
            hooks[ii].remove()
    
    
    return all_activ


def save_pca_activ(all_activ, save_activ_dir, model_name='AlexNet', save_big_activ=False):

    ## Take the output of "get_alexnet_activations", restructure it into big arrays, do PCA and save.
    
    # figure out number of images/number of batches based on dim of the input
    nbatches = len(all_activ)
    batch_size = np.shape(all_activ[0][0])[0]
    nims = (nbatches-1)*batch_size + np.shape(all_activ[nbatches-1][0])[0]
    nfeaturelayers = len(all_activ[0])
    
    if model_name=='AlexNet':
        layer_names = ['1_Conv2d','1_ReLU','1_MaxPool2d',
                       '2_Conv2d','2_ReLU','2_MaxPool2d',
                       '3_Conv2d','3_ReLU',
                       '4_Conv2d','4_ReLU',
                       '5_Conv2d','5_ReLU','5_MaxPool2d']
    else:
        raise ValueError('the model has to be alexnet for now')    
        
    # looping over layers, will save separate files for each layer
    for ll in range(nfeaturelayers):

        # going to loop over batches and put the whole big matrix back together 
        activ= np.zeros([nims, np.shape(all_activ[0][ll])[1], 
                             np.shape(all_activ[0][ll])[2],np.shape(all_activ[0][ll])[3]])
        for bb in range(nbatches):
            activ[bb*batch_size:(bb+1)*batch_size,:,:,:] = all_activ[bb][ll]

      
        if save_big_activ:
            # note these are very big files, only save if you have a reason for wanting them later
            print('saving activs for layer %s, size is:\n'%layer_names[ll])
            print(np.shape(activ))
            fn2save = os.path.join(save_activ_dir, '%s_%s_fullactiv.npy'%(model_name,layer_names[ll]))
            np.save(fn2save, activ)
            print('saved')

        # Now PCA
        # first, reshape to [nIms x nUnits]
        # disregarding difference between channels/spatial dims for now
        nUnitsTotal = np.prod(np.shape(activ)[1:])
        activ_full = np.reshape(activ,[nims, nUnitsTotal])
        np.shape(activ_full)

        # reduce dimensionality of whole matrix w PCA
        pca = decomposition.PCA()
        print('running pca (original size %d by %d)'%(np.shape(activ_full)[0], np.shape(activ_full)[1]))
        pca.fit(activ_full)
        scores = pca.transform(activ_full)
        pct_expl_var = pca.explained_variance_
        pct_expl_var = pct_expl_var/sum(pct_expl_var)*100
        print('finished')   
        fn2save = os.path.join(save_activ_dir, '%s_%s_pca.npy'%(model_name, layer_names[ll]))
        print('saving to %s'%fn2save)
        torch.save({'scores': scores,'pct_expl_var': pct_expl_var },
                    fn2save)
        print('saved')
        
        
def get_dataset_splits(image_dir, transform=None, rndseed=None):
   
    ## Make dataset and split into training/validation/testing sets
    
    csv_file = os.path.join(image_dir, 'shape_labels_all.csv')
    
    if transform==None:
        # create a transformation that will normalize the image intensity (z-score)
        # these values were copied from https://pytorch.org/tutorials/beginner/finetuning_torchvision_models_tutorial.html
        # supposed to match pre-training data for these models
        transform = transforms.Compose([transforms.ToTensor(), 
                                        transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])])
        
    sds = shape_dataset(csv_file, image_dir, transform = transform)
    
    nTotal = len(sds)
    nTst = int(0.10*nTotal)  # this part is left out throughout all of training, used to test at very end
    nVal = int(0.10*(nTotal-nTst)) # this part is left out during training but also used to evaluate loss throughout training
    nTrn = nTotal - nTst - nVal  # this part is for training
    
    if rndseed==None:
        rndseed = 904589  # make sure we do this the same way every time  

    # random_split will create 3 new datasets, each of which is a subset of shape_dataset - have same fields etc but only subset of items.
    trn, val, tst = torch.utils.data.random_split(sds, [nTrn, nVal, nTst], generator=torch.Generator().manual_seed(rndseed))
    
    return trn, val, tst
  
  
class shape_classifier(nn.Module):
    
    ## Create a module to do binary shape classification - will use for fine-tuning
    # for now this is basically same as torch.nn.Linear, but might want to change it later?
    
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.weights = nn.Parameter(torch.randn(in_channels, out_channels) / math.sqrt(in_channels))
        self.bias = nn.Parameter(torch.zeros(out_channels))

    def forward(self, xb):

        return xb.double() @ self.weights.double() + self.bias.double()
      
      
def fine_tune_alexnet_binary(image_dir, save_model_dir, task, 
                             maxsteps=201, val_every_nsteps=10, lr=0.0001, rndseed=None):

    # Using pre-trained AlexNet, will replace last layer w binary classifier and 
    # fine tune for a particular binary task (specified by task=[0,1,2])

    if not np.isin(task, [0,1,2]):
        raise ValueError('task has to be 0, 1 or 2 to indicate the categorization scheme')

    # first get the dataset splits
    trn, val, tst  = get_dataset_splits(image_dir)

    # define the pre-trained model
    model = models.alexnet(pretrained=True).double()  # has to be double otherwise throws error
    model.train()  # use train mode to be able to get gradients

    # create the linear classifier that will read out a shape category
    if rndseed==None:
        rndseed = 993494       
    torch.manual_seed(rndseed) # using a manual seed to initializing weights of last layer, for reproducibility
    clf = shape_classifier(4096,2)
    clf.weights.requires_grad_() # make sure these are trainable
    clf.bias.requires_grad_()

    # replace last layer of original model with this new classifier
    model.classifier[6] = clf

    # create a dataloader object that breaks training set into shuffled batches
    batch_size=100
    trnloader = DataLoader(trn, batch_size = batch_size, shuffle=True)
    
    # same for validation set (doing this in one batch for now, but could split into multiple)
    val_batch_size = len(val) 
    valloader = DataLoader(val, batch_size=val_batch_size, shuffle=True)

    # counting total passes over the entire training set (nepochs = nsteps/nbatches)
    nepochs_done = 0 
    # create an iterator object that loops over training set batches (this gets re-randomized every epoch)
    trngenerator = iter(trnloader)

    # loop over training steps 
    # one step=one weight update, based on one batch of data
    for step in range(maxsteps):

      # first get the next batch of data from trn set
      try:
          # try taking next sample from current iterator object
          curr_batch_trn = next(trngenerator)
      except StopIteration:
          # if the previous iterator object is exhausted, then the epoch is finished, 
          # need to make a new iterator for next epoch.
          trngenerator = iter(trnloader)
          curr_batch_trn = next(trngenerator)
          nepochs_done = nepochs_done+1;
          print('starting epoch %d'%nepochs_done)

      # get image and labels for this batch
      # image_tensors is [batch_size x 3 x 224 x 224]
      image_tensors =  curr_batch_trn['image']
      # task_labels is [batch_size x 3] where columns are diff tasks  
      task_labels = curr_batch_trn['task_labels'][:,task] 

      # do the forward pass of model 
      model.train() 
      out = model(image_tensors)
      shape_output = out

      # computing some basic performance metrics
      pred_labels = torch.argmax(shape_output, 1)
      acc = torch.sum(torch.eq(pred_labels,task_labels))/len(task_labels)
      loss = F.cross_entropy(shape_output, task_labels)

      # print out the performance of the model, just before this weight updating step was performed
      print('step %d, trn loss = %.5f, trn accuracy = %.5f'%(step,loss.item(),acc))

      # print out an example gradient before loss.backward is run (all grads should be zero here)
      # if step>0:
        # print('before update: grad is %.6f, weight is %.6f'%(clf.weights.grad[0,0].item(),clf.weights[0,0].item()))
        # print('before update: grad is %.6f, weight is %.6f'%(model.features[0].weight.grad[0,0,0,0].item(),model.features[0].weight[0,0,0,0].item()))
        # # for any layers after the one we attached the classifier to, this should fail since the gradient for those layers is not defined during the backward sweep
        # print('before update: grad is %.6f, weight is %.6f'%(model.classifier[1].weight.grad[0,0].item(),model.classifier[1].weight[0,0].item()))

      # run the graph backward, to get gradient values. 
      loss.backward()

      with torch.no_grad():
        # here is where we update the weights for the final layer, based on grad values
        clf.weights -= clf.weights.grad * lr
        clf.bias -= clf.bias.grad * lr

      # print out the same example weight just after updating (this is to double check that only desired layer weights change)
      # if step>0:
        # print('after update: grad is %.6f, weight is %.6f'%(clf.weights.grad[0,0].item(),clf.weights[0,0].item()))
        # print('after update: grad is %.6f, weight is %.6f'%(model.features[0].weight.grad[0,0,0,0].item(),model.features[0].weight[0,0,0,0].item()))
        # print('after update: grad is %.6f, weight is %.6f'%(model.classifier[1].weight.grad[0,0].item(),model.classifier[1].weight[0,0].item()))

      # zero out the gradients before moving to the next step
      with torch.no_grad():
        model.zero_grad()
        clf.weights.grad.zero_()
        clf.bias.grad.zero_()

      # decide whether to validate this step (not every step, since it takes a little while)
      if step>0 and (np.mod(step, val_every_nsteps)==0 or step==maxsteps-1):

        # do validation    
        curr_batch_val = next(iter(valloader))
        val_image_tensors =  curr_batch_val['image']
        val_task_labels = curr_batch_val['task_labels'][:,task] 

        model.eval() # set to eval mode before we run val set through, so we don't mess up gradients
        with torch.no_grad():
          val_out = model(val_image_tensors)
          val_pred_labels = torch.argmax(val_out, 1)
          val_acc = torch.sum(torch.eq(val_pred_labels,val_task_labels))/len(val_task_labels)
          val_loss = F.cross_entropy(val_out, val_task_labels)

        print('step %d, val loss = %.5f, val accuracy = %.2f'%(step,val_loss.item(),val_acc))
        model.train() # set back to trn mode
    
    # save model after the finetuning
    fn2save = os.path.join(save_model_dir, 'AlexNet_finetune_task%d.pt'%task)
    print('saving to %s ...'%fn2save)
    torch.save({'val_acc': val_acc, 
                'val_loss': val_loss.item(),
                'trn_loss': loss.item(),
                'trn_acc': acc,
                'step': step,
                'epoch': nepochs_done,            
                'model_state_dict': model.state_dict()}, fn2save)
    print('saved')
    
    
    
def eval_alexnet_fine_tune(image_dir, save_model_dir, task):

    ## Load saved model after fine-tuning last layer, and evaluate on the test set. 
    # Print accuracy/loss.
    
    # first get the dataset splits 
    # this function has a fixed default random seed, if you don't over-ride that seed then this splitting 
    # will be identical to how it was during training.
    trn, val, tst  = get_dataset_splits(image_dir)
    
    model = models.alexnet(pretrained=True).double()
    # have to make the model exactly same structure as how it was initially created
    # so put in a shape classifier module like we did before (if skip this the state dict won't load)
    clf = shape_classifier(4096,2)
    clf.weights.requires_grad_()
    clf.bias.requires_grad_()
    model.classifier[6] = clf

    fn2load = os.path.join(save_model_dir, 'AlexNet_finetune_task%d.pt'%task)
    print('loading from %s\n'%fn2load)
    ckpt = torch.load(fn2load)
    model.load_state_dict(ckpt['model_state_dict'])

    model.eval()

    tst_batch_size = len(tst)
    tstloader = DataLoader(tst, batch_size=tst_batch_size, shuffle=False)
    curr_batch_tst = next(iter(tstloader))
    tst_image_tensors =  curr_batch_tst['image']
    tst_task_labels = curr_batch_tst['task_labels'][:,task] 
#     print(tst_task_labels[0:10])
    model.eval()
    with torch.no_grad():
      tst_out = model(tst_image_tensors)
      tst_pred_labels = torch.argmax(tst_out, 1)
      tst_acc = torch.sum(torch.eq(tst_pred_labels,tst_task_labels))/len(tst_task_labels)
      tst_loss = F.cross_entropy(tst_out, tst_task_labels).item()

    print('step %d, val loss = %.5f, val accuracy = %.2f'%(ckpt['step'],ckpt['val_loss'],ckpt['val_acc']))
    print('step %d, test loss = %.5f, test accuracy = %.2f'%(ckpt['step'],tst_loss, tst_acc))