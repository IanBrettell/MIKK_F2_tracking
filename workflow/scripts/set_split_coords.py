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
import pandas as pd
import numpy as np
import cv2 as cv
import os

# Get variables

## Debugging
IN_FILE = "/hps/nobackup/birney/users/ian/MIKK_F2_tracking/recoded/20211117_1326_R.avi"
SAMPLES_FILE = "config/samples_long.csv"
ASSAY = "open_field"
VIDEO = "20211117_1326_R"

## True
IN_FILE = snakemake.input.video[0]
SAMPLES_FILE = snakemake.params.samples_file
ASSAY = snakemake.params.assay
VIDEO = snakemake.params.video
OUT_FILE = snakemake.output.fig

# Read samples_file
samples_df = pd.read_csv(SAMPLES_FILE)

# Get movie name

MOVIE_NAME = VIDEO + ".avi"

# Get target rows

rows = samples_df[
    (samples_df['movie_name'] == MOVIE_NAME) & \
    (samples_df['assay'] == ASSAY)
]

# Get start frame and crop adjustment values
## note: Negative values for top/bottom shift boundary up
## note: Negative values for left/right shift boundary left
start = int(rows['frame_start'].max())
adj_top = int(rows['adj_top'].max())
adj_right = int(rows['adj_right'].max())

# Read video from file
cap = cv.VideoCapture(IN_FILE)

# Frame width and height
wid = int(cap.get(cv.CAP_PROP_FRAME_WIDTH))
hei = int(cap.get(cv.CAP_PROP_FRAME_HEIGHT))

# Set adjusted midpoints
mid_x = round(((wid - 1) / 2) + adj_right)
mid_y = round(((hei - 1) / 2) + adj_top)

# Capture start frame
cap.set(cv.CAP_PROP_POS_FRAMES, start)

# Read frame
ret, frame = cap.read()

# Add vertical line 
start_point = (mid_x, 0)
end_point = (mid_x, hei)
color = (62,90,248)
thickness = 1
frame = cv.line(frame, start_point, end_point, color, thickness)

# Add horizontal line
start_point = (0, mid_y)
end_point = (wid, mid_y)
color = (62,90,248)
thickness = 1
frame = cv.line(frame, start_point, end_point, color, thickness)

# Get width of boundaries in pixels
bleft = int(rows["bound_left"].max())
bright = int(rows["bound_right"].max())
btop = int(rows["bound_top"].max())
bbottom = int(rows["bound_bottom"].max())

# Get line coordinates
left_start = (bleft, 0)
left_end = (bleft, hei)
right_start = (wid - bright, 0)
right_end = (wid - bright, hei)
top_start = (0, btop)
top_end = (wid, btop)
bottom_start = (0, hei - bbottom)
bottom_end = (wid, hei - bbottom)

# Vertical lines
frame = cv.line(frame, left_start, left_end, color, thickness)
frame = cv.line(frame, right_start, right_end, color, thickness)

# Horizontal lines
frame = cv.line(frame, top_start, top_end, color, thickness)
frame = cv.line(frame, bottom_start, bottom_end, color, thickness)

# Write frame
cv.imwrite(OUT_FILE, frame)

cap.release()
