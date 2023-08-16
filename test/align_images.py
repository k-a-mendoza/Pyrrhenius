import numpy as np
from scipy.ndimage import correlate
from PIL import Image


def _get_image_array(location):
    img = 255-np.asarray(Image.open(location).convert('L')).astype(float)
    img/=255
    return img

def get_corner_weight_masks():
    corner_1 = np.zeros((200,200))# upper left
    corner_1[:,:2]=1
    corner_1[:2,:]=1

    corner_2 = np.zeros((200,200))# lower right
    corner_2[:,-2:]=1
    corner_2[-2:,:]=1
    corner_3 = np.zeros((200,200)) # upper right
    corner_3[:,-2:]=1
    corner_3[:3,:]=1
    corner_4 = np.zeros((200,200)) # lower left
    corner_4[:,:3]=1
    corner_4[-2:,:]=1
    return corner_4, corner_2, corner_1, corner_3

def get_corner_indices(img_location):
    image_array = _get_image_array(img_location)
    shape = image_array.shape
    corner_4, corner_2, corner_1, corner_3 =  get_corner_weight_masks()
    lower_left  = correlate(image_array,corner_4, mode='constant')
    lower_left = np.unravel_index(lower_left.argmax(), shape)

    lower_right = correlate(image_array,corner_2,  mode='constant')
    lower_right = np.unravel_index(lower_right.argmax(), shape)
    
    upper_left = correlate(image_array,corner_1,  mode='constant')
    upper_left = np.unravel_index(upper_left.argmax(), shape)

    upper_right = correlate(image_array,corner_3,  mode='constant')
    upper_right = np.unravel_index(upper_right.argmax(), shape)
    return lower_left, lower_right, upper_left, upper_right, image_array

# image to real coordinates:

def get_y_extent(img,y0,y1,yprime0,yprime1):
    ymin = y0 + (y1-y0)/(yprime1 - yprime0)*-yprime0
    ymax = y0 + (y1-y0)/(yprime1 - yprime0)*(img.shape[0]-yprime0)
    return ymin, ymax

def get_x_extent(img,x0,x1,xprime0,xprime1):
    xmin = x0 + (x1-x0)/(xprime1 - xprime0)*-xprime0
    xmax = x0 + (x1-x0)/(xprime1 - xprime0)*(img.shape[1]-xprime0)
    return xmin, xmax

def get_extent(x0,x1,y0,y1,img_location):
    lower_left, lower_right, upper_left, upper_right, image_array = get_corner_indices(img_location)
    xprime0 = lower_left[1]-100
    xprime1 = lower_right[1]+100

    yprime0 = lower_left[0]+100
    yprime1 = upper_right[0]-100

    ymin, ymax = get_y_extent(image_array,y0,y1,yprime0,yprime1)
    xmin, xmax = get_x_extent(image_array,x0,x1,xprime0,xprime1)
    return xmin, xmax, ymax, ymin