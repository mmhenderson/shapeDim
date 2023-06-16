import sys, os
import numpy as np
import time
from sklearn import decomposition
import pandas as pd
import torch
import PIL.Image
import h5py

# root directory is 2 dirs up from this file
path = os.path.realpath(__file__).split('/')
root = '/'+os.path.join(*path[0:-3])


def prep_brick(debug=False):
    
    print('debug=%s'%debug)
    # make image brick that is used by get_gist.m
    feat_dir = os.path.join(root, 'Analysis', 'image_similarity', 'features')
    image_dir = os.path.join(root, 'Stimuli','AmpGrid3_adj_full_grey_small')
    
    shape_df = make_shape_labels(image_dir)
    image_data = load_images(shape_df, debug=debug)

    fn2save = os.path.join(feat_dir,'Images_grid3_all.h5py')
    print('saving to %s'%fn2save)

    with h5py.File(fn2save, 'w') as f:

        f.create_dataset("stimuli", data=image_data, dtype='i')
        f.close()

def load_images(shape_df, debug=False):
    
    n_images = shape_df.shape[0]
    if debug:
        n_images_load = 10
    else:
        n_images_load = n_images
        
    first_image = PIL.Image.open(shape_df['filename_full'][0])
    n_pix = first_image.size[0]
    
    # [images x color_channels x height x width]
    image_array = np.zeros((n_images,1,n_pix,n_pix),dtype=np.float32)
    
    for ii in range(n_images_load):
        
        im = PIL.Image.open(shape_df['filename_full'][ii])
        imdat = np.reshape(np.array(im.getdata()), im.size)
        
        image_array[ii,0,:,:] = imdat
        
    return image_array



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
    labels_task1[np.isin(all_quadrant,[1,4])] = 2
    labels_task1[np.isin(all_quadrant,[2,3])] = 1

    labels_task2 = np.zeros(np.shape(all_quadrant))
    labels_task2[np.isin(all_quadrant,[1,2])] = 2
    labels_task2[np.isin(all_quadrant,[3,4])] = 1

    labels_task3 = np.zeros(np.shape(all_quadrant))
    labels_task3[np.isin(all_quadrant,[1,3])] = 2
    labels_task3[np.isin(all_quadrant,[2,4])] = 1

    filenames_full = [os.path.join(image_dir, 'Shape_%.2f_%.2f.png'%(x,y)) \
                      for x,y in zip(all_grid_points[:,0], all_grid_points[:,1])]
    
    shape_df = pd.DataFrame.from_dict({'coord_axis1':all_grid_points[:,0], 'coord_axis2':all_grid_points[:,1],
     'quadrant':np.squeeze(all_quadrant), 'labels_task1':np.squeeze(labels_task1),
     'labels_task2':np.squeeze(labels_task2), 'labels_task3':np.squeeze(labels_task3),
     'filename_full': filenames_full})

    return shape_df
  