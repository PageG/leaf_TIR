#!/usr/bin/python


# USAGE:
#  ./run_me.py COMPOSITE_IMAGE     
#      -- or --
#  ./run_me.py THERMAL_IMAGE COLOR_IMAGE
# 
# for example,
#  ./run_me.py FLIR1544.jpg


import cv2
import numpy as np
import os, sys

TMP_OUTPUT = True # Should we store temporary images?

if (len(sys.argv) == 3):
    ## supply the color and thermal images separately
    therm_img_path = sys.argv[1]
    color_img_path = sys.argv[2]
else:
    ## extract the color image from the thermal image
    therm_img_path = sys.argv[1]
    color_img_path = os.path.splitext(therm_img_path)[0]+'_color.jpg'
    os.system("exiftool -b -EmbeddedImage "+therm_img_path+" > "+color_img_path)
    print("Exported color data to "+color_img_path)

# Read the color image
im1 =  cv2.imread(color_img_path);

# This processes the thermal image to output it in relative scale:
# 1) Extract the thermal image with exiftool
# 2) correct the quirks in the FLIR data file (convert+awk)
# 3) scale it in a relative scale, output to png, and read it back
tmp_therm_img_path = os.path.splitext(therm_img_path)[0]+'_thermal.pgm'
new_therm_img_path = os.path.splitext(therm_img_path)[0]+'_thermal.png'
#os.system("exiftool -b -RawThermalImage "+therm_img_path+" | convert - -compress none "+tmp_therm_img_path)
#os.system("awk -f png.txt "+tmp_therm_img_path+" | convert -  -auto-level "+new_therm_img_path)
os.system("exiftool -b -RawThermalImage "+therm_img_path+" | convert - -compress none "+tmp_therm_img_path)
os.system("awk -f png.txt "+tmp_therm_img_path+" | convert -  -auto-level "+new_therm_img_path)
os.system("rm "+tmp_therm_img_path)
# os.system("exiftool -b -RawThermalImage "+therm_img_path+" > "+new_therm_img_path)
print("Exported thermal data to "+new_therm_img_path+ " (dont forget: the thermal scale is relative)")
im2 =  cv2.imread(new_therm_img_path)

# now some manual transformations to *roughly* align the color image with the thermal data
fx = fy = 0.55 
cx = 550.; cy = 540. 
im1 = cv2.resize(im1,None,fx=fx, fy=fy, interpolation = cv2.INTER_CUBIC)
M = np.float32([[1,0,-cx/(1./fx)],[0,1,-cy/(1./fy)]])
rows,cols,_ = im1.shape
im1 = cv2.warpAffine(im1,M,(cols,rows))
im1 = im1[0:480, 0:640]

# Convert images to grayscale
im1_gray = cv2.cvtColor(im1,cv2.COLOR_BGR2GRAY)
im2_gray = cv2.cvtColor(im2,cv2.COLOR_BGR2GRAY)
if TMP_OUTPUT:
    cv2.imwrite('tmp1a.jpg', im1_gray)
    cv2.imwrite('tmp2a.jpg', im2_gray)

# Generate Laplacian filtered images (edge detection)
im1_gray = cv2.convertScaleAbs(cv2.Laplacian(im1_gray,3),10)
im2_gray = cv2.convertScaleAbs(cv2.Laplacian(im2_gray,3),10)
if TMP_OUTPUT:
    cv2.imwrite('tmp1b.jpg', im1_gray)
    cv2.imwrite('tmp2b.jpg', im2_gray)
 
# Get the size of image (should be 640x480 !)
sz = im1.shape
 
# Define the motion model
warp_mode = cv2.MOTION_HOMOGRAPHY # homography gives best results
 
# Define 2x3 or 3x3 matrices and initialize the matrix to identity
if warp_mode == cv2.MOTION_HOMOGRAPHY :
    warp_matrix = np.eye(3, 3, dtype=np.float32)
else :
    warp_matrix = np.eye(2, 3, dtype=np.float32)
 
# Specify the number of iterations.
number_of_iterations = 1000; # maybe could decrease (long for now)
 
# Specify the threshold of the increment in the correlation coefficient between two iterations
termination_eps = 1e-10;
 
# Define termination criteria
criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, number_of_iterations,  termination_eps)
 
# Run the ECC algorithm. The results are stored in warp_matrix.
(cc, warp_matrix) = cv2.findTransformECC(im2_gray,im1_gray,warp_matrix, warp_mode, criteria)
 
if warp_mode == cv2.MOTION_HOMOGRAPHY :
    # Use warpPerspective for Homography 
    im1_aligned = cv2.warpPerspective(im1, warp_matrix, (sz[1],sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP)
else :
    # Use warpAffine for Translation, Euclidean and Affine
    im1_aligned = cv2.warpAffine(im1, warp_matrix, (sz[1],sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP);

if TMP_OUTPUT:
    cv2.imwrite('tmp1c.jpg', cv2.warpAffine(im1_gray, np.resize(warp_matrix, (2,3)), (sz[1],sz[0]), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP))
    print('Homography matrix:\n')
    print(str(warp_matrix))

new_color_img_path = os.path.splitext(color_img_path)[0]+'_aligned.png'
cv2.imwrite(new_color_img_path, im1_aligned)
print("Saved aligned color image to "+new_color_img_path)
 
# Show intermediate and final results
#cv2.imshow("Color Image", im1)
#cv2.imshow("Thermal Image", im2)
#cv2.imshow("Edges in Color Image", im1_gray)
#cv2.imshow("Edges if Thermal Image", im2_gray)
#cv2.imshow("Aligned Color Image", im1_aligned)
#cv2.waitKey(0)
