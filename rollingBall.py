import numpy as np
import math
from scipy.ndimage import zoom
from scipy.ndimage import uniform_filter,minimum_filter
# from tqdm import tqdm
import scipy.ndimage as ndi
from skimage.measure import block_reduce
from tqdm.notebook import tqdm_notebook as tqdm
def normalize_local_contrast(image, blockRadiusX, blockRadiusY, meanFactor, center, stretch):
    image_original = np.float32(image)
    if image.ndim == 2:
        block_size = (2*blockRadiusY+1, 2*blockRadiusX+1)
    else:
        block_size = (1,2*blockRadiusY+1, 2*blockRadiusX+1)
    mean = uniform_filter(image_original, block_size)
    squared_mean = uniform_filter(image_original**2, block_size)
    std = np.sqrt(squared_mean - mean**2)
    d = meanFactor * std

    if center and stretch:
        min_ = mean - d
        return ((image_original - min_) / (2 * d) * image_original.ptp() + image_original.min())
    elif center:
        return image_original - mean + image_original.mean()
    elif stretch:
        return (image_original - mean) / (2 * d) * image_original.ptp() + mean


def shrink_image(ip, shrink_factor):
    # Use skimage's block_reduce function to find the minimum of each block
    small_image = block_reduce(ip, block_size=(shrink_factor, shrink_factor), func=np.min, cval=np.inf)
    return small_image

def enlarge_image(small_image, ip, shrink_factor):
    height, width = ip.shape
    small_height, small_width = small_image.shape

    x_small_indices, x_weights = make_interpolation_arrays(width, small_width, shrink_factor)
    y_small_indices, y_weights = make_interpolation_arrays(height, small_height, shrink_factor)

    # Precalculate these values to avoid repeated computation
    inv_x_weights = 1.0 - x_weights

    line0 = np.zeros(width, dtype=np.float32)
    line1 = np.zeros(width, dtype=np.float32)
    # Vectorized x-interpolation of the first smallImage line
    line1 = small_image[0, x_small_indices] * x_weights + small_image[0, x_small_indices + 1] * inv_x_weights

    y_small_line0 = -1  # line0 corresponds to this y of smallImage
    for y in range(height):
        if y_small_line0 < y_small_indices[y]:
            # swap lines
            line0, line1 = line1, line0
            y_small_line0 += 1
            s_y_pointer = y_small_indices[y] + 1  # points to line0 + 1 in smallImage
            # Vectorized x-interpolation of the new smallImage line -> line1
            line1 = small_image[s_y_pointer, x_small_indices] * x_weights + small_image[s_y_pointer, x_small_indices + 1] * inv_x_weights
        weight = y_weights[y]
        inv_weight = 1.0 - weight
        ip[y] = line0 * weight + line1 * inv_weight  # Vectorized y-interpolation
    return ip


def make_interpolation_arrays(length, small_length, shrink_factor):
    small_indices = np.zeros(length, dtype=np.int32)
    weights = np.zeros(length, dtype=np.float32)
    for i in range(length):
        small_index = (i - shrink_factor // 2) // shrink_factor
        if small_index >= small_length - 1:
            small_index = small_length - 2
        small_indices[i] = small_index
        distance = (i + 0.5) / shrink_factor - (small_index + 0.5)  # distance of pixel centers (in smallImage pixels)
        weights[i] = 1.0 - distance
    return small_indices, weights

class RollingBall:
    def __init__(self, radius):
        self.radius = radius
        self.ball = self.buildRollingBall(radius)


    @staticmethod
    def buildRollingBall(radius):
        y, x = np.ogrid[-radius: radius+1, -radius: radius+1]
        ball = x**2 + y**2 <= radius**2
        return ball.astype(float)


def getShrinkFactor(radius):
    if (radius<=10):
        return 1
    if (radius<=30):
        return 2
    elif (radius<=100):
        return 4
    return 8
    

def rolling_ball_float_background(fp, radius, invert=False,ball=None,shrink_factor=None):
    # print("hi1")
    if ball == None:
        ball = RollingBall(radius)
    if shrink_factor ==None:
        shrink_factor = getShrinkFactor(radius)
    if len(fp.shape) ==3:
        # pixels = [print(x) for x in tqdm(fp)]
        # print("hi2")
        pixels = [rolling_ball_float_background(x, radius, invert,ball) for x in tqdm(fp)]
        return pixels
    # print("hi3")
    pixels = np.array(fp).astype(np.float32).copy()  # make a copy
    if invert:
        pixels = -pixels

    shrunk = shrink_image(pixels, shrink_factor) if shrink_factor > 1 else pixels
    shrunk = rollBall(ball, shrunk)
    pixels = enlarge_image(shrunk, pixels, shrink_factor) if shrink_factor > 1 else pixels

    if invert:
        pixels = -pixels
    return fp-pixels


def rollBall(ball, image):
    ball = ball.ball
    height, width = image.shape
    ball_width = ball.shape[0]
    radius = ball_width // 2

    # Subtracts the minimum value within the sliding window defined by the ball from each pixel in the image.
    final_image = minimum_filter(image, footprint=ball, mode='constant', cval=np.inf) - ball.min()

    return final_image

