import os
import glob

import matplotlib.pyplot as plt
import imageio

def pngs_to_mp4(path_to_pngs, video_filename, fps=5):
    """
    Converts Portable Network Graphic (PNG) files into an MP4 video.
    The order of the files depends on the file names, in alphabetical order.

    Parameters
    ---------
    path_to_pngs (string)
        A string representing the relative or absolute path to 
        a directory containing PNG files.

    video_filename (string)
        A string representing the name of the mp4 file to be saved.

    fps (int)
        An integer number of frames per second.

    """
    # Check that `path_to_pngs` is indeed a directory, and not, say, a file.
    assert(os.path.isdir(path_to_pngs))

    # Get a list of all the png files.
    png_files = sorted(glob.glob(os.path.join(path_to_pngs, '*.png')))

    # Ensure video_filename ends with .mp4
    if not video_filename.endswith(".mp4"): video_filename += ".mp4"

    # Make directory if needed
    video_directory = os.path.dirname(video_filename)
    if not os.path.isdir(video_directory): os.makedirs(video_directory)

    # Loop over files to make video.
    with imageio.get_writer(video_filename, fps=fps) as writer:
        for file in png_files:
            image = imageio.imread(file)
            writer.append_data(image)

