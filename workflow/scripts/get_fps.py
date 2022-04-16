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

# Get variables

## Debug
#IN = ["/hps/nobackup/birney/users/ian/MIKK_F2_tracking/raw_videos/20211117_1326_R.avi",
#      "/hps/nobackup/birney/users/ian/MIKK_F2_tracking/raw_videos/20220405_1005_R.avi",
#      "/hps/nobackup/birney/users/ian/MIKK_F2_tracking/raw_videos/20220323_1421_L.avi"]

## True
IN = snakemake.input
OUT = snakemake.output[0]

# Sort file names

IN = sorted(IN)

# Write lines to file and close

file = open(OUT, 'w')

for VID in IN:
    cap = cv.VideoCapture(VID)
    fps = str(int(cap.get(cv.CAP_PROP_FPS)))
    sample = os.path.basename(VID).replace('.avi', '')
    file.writelines(sample + ',' + fps + '\n')

file.close()



