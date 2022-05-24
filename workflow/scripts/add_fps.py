# Send stdout and stderr to log file
import sys,os
import logging, traceback
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(exc_type, exc_value, exc_traceback)
                         ])
                 )
# Install exception handler
sys.excepthook = handle_exception

# Import libraries

import cv2 as cv
import pandas as pd

# Get variables

## Debug
#IN_VIDS = ["/hps/nobackup/birney/users/ian/MIKK_F2_tracking/raw_videos/20211117_1326_R.avi",
#      "/hps/nobackup/birney/users/ian/MIKK_F2_tracking/raw_videos/20220405_1005_R.avi",
#      "/hps/nobackup/birney/users/ian/MIKK_F2_tracking/raw_videos/20220323_1421_L.avi"]
#SAMPLES_CONV = "config/samples_converted.csv"

## True
IN_VIDS = snakemake.input.videos
SAMPLES_CONV = snakemake.input.samples_file[0]
OUT = snakemake.output[0]

# Sort file names

IN_VIDS = sorted(IN_VIDS)

# Create DF with FPS and sample name

FPS = []
VIDEOS = []
for VID in IN_VIDS:
    cap = cv.VideoCapture(VID)
    fps = str(int(cap.get(cv.CAP_PROP_FPS)))
    video = os.path.basename(VID)
    # add to lists
    FPS.append(fps)
    VIDEOS.append(video)

fps_df = pd.DataFrame({
    'fps' : FPS,
    'movie_name' : VIDEOS
})

# Read in converted metadata
meta_df = pd.read_csv(SAMPLES_CONV)

out = meta_df.set_index('movie_name').join(fps_df.set_index('movie_name'))

# Write to file

out.to_csv(OUT)



