## 14 April 2022

Alicia sent through updated metadata and new videos:

`config/New\ Fish\ Movies.xlsx`

## 20220222

Previous tracking run had 121/274 (44%) perfect tracking, and 182/274 (66%) above 99%.

After adjusting the parameters for each of the remaining 92 videos, I am uncertain we can obtain further improvements in the tracking after this second round.

Videos that will likely continue to track poorly are listed in `config/poor_segmentation.csv`, with notes for each.

## 20220221

Background subtraction usually works better with the novel object. If not used, the fish object can disappear when it moves close to the object.

Best parameters:

OF:
* bgsub: False
* intensity: [0, 155]
* area: [20, 500]

NO:
* bgsub: True
* intensity: [0, 130]
* area: [20, 500]

The videos likely to segment poorly are where:
* One fish remains still, which creates a shadow when bgsub = True; and
* The other fish moves close to or behind the novel object which makes it disappear.