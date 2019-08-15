# Imaging_tania_project

-We record all the tracks for each channels without allowing gaps (trackmate)

-We stitch trajectories within a channel using the following information:
1) We privilege the tracks where the distances across channels are small
2) We privilege to stitch together tracks  that are close (distance btw centre of mass) within a given channels
3) We privilege to stitch tracks that are close in time

- run on all files within a directory or directories in a directory
- save stitched tracks and id of tracks that are used for each cell
- filter out timepoints with high movement speed (derivative of space
with respect to time)
- save quality control for each cell
