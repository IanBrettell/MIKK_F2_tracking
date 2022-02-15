# To run on Codon:
# ssh codon
# module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
# bsub -M 20000 -q datamover -Is bash
# cd /hps/software/users/birney/ian/repos/MIKK_F0_tracking
# conda activate opencv_4.5.5.yaml
# python workflow/scripts/join_20191120_1217_106-2_R_C.py  

# Add packages

import cv2 as cv

# Set variables

IN_1 = "/nfs/ftp/private/indigene_ftp/upload/behaviour/transfer/20191111_panel_1/20191120/20191120_1217_106-2_R_C.avi"
IN_2 = "/nfs/ftp/private/indigene_ftp/upload/behaviour/transfer/20191111_panel_1/20191120/20191120_1217_106-2_R_C_part2.avi"
OUT_FILE = "/nfs/ftp/private/indigene_ftp/upload/behaviour/transfer/20191111_panel_1/all_to_analyse/20191120_1217_106-2_R_C.avi"
FPS = 30

# Write first part

cap = cv.VideoCapture(IN_1)

start = 0
end = int(cap.get(cv.CAP_PROP_FRAME_COUNT))
wid = int(cap.get(cv.CAP_PROP_FRAME_WIDTH))
hei = int(cap.get(cv.CAP_PROP_FRAME_HEIGHT))
size = (wid, hei)

# Define the codec and create VideoWriter object

fourcc = cv.VideoWriter_fourcc('h', '2', '6', '4')
out = cv.VideoWriter(OUT_FILE, fourcc, FPS, size, isColor=True)


i = start
while i in range(start,end):
    cap.set(cv.CAP_PROP_POS_FRAMES, i)
    # Capture frame-by-frame
    ret, frame = cap.read()
    # if frame is read correctly ret is True
    if not ret:
        print("Can't receive frame (stream end?). Exiting ...")
        break
    # Write frame
    out.write(frame)
    # Add to counter
    i += 1

cap.release()

# Write second part

cap = cv.VideoCapture(IN_2)

start = 0
end = int(cap.get(cv.CAP_PROP_FRAME_COUNT))

i = start
while i in range(start,end):
    cap.set(cv.CAP_PROP_POS_FRAMES, i)
    # Capture frame-by-frame
    ret, frame = cap.read()
    # if frame is read correctly ret is True
    if not ret:
        print("Can't receive frame (stream end?). Exiting ...")
        break
    # Write frame
    out.write(frame)
    # Add to counter
    i += 1

cap.release()
out.release()
