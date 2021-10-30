# Single-pulse candidate viewer #

This directory contains a simple viewer for single-pulse candidate plots as generated by the MeerTRAP detection pipeline. It runs on the command line and uses `ImageMagick`'s `display` tool to show the candidate plots.

## Requirements ##

* `bash` shell
* [ImageMagick](https://imagemagick.org/)

## Usage ##

Download the script, mark it as executable, and copy it to some directory that is in your $PATH. Restart your terminal/shell. If you are vetting MeerTRAP candidates, navigate to the directory that contains the candidate plots to be reviewed. Execute the script in the candidate directory in your terminal. The first candidate plot will open in a new window, and you will see the candidate number, filename, a progress indicator, and other information in your terminal. Select the candidate plot window and hit the `q` key on your keyboard to advance to the next candidate plot. The candidate window will automatically close once you have finished vetting all candidates. If you want to check a previous candidate plot, simply copy its filename from the status list in your terminal, and open it in a separate terminal using `display $candidatefile`.

The script shows you the candidate plots in decreasing DM order by default, but that can be changed easily in the code. It appends the filenames of each reviewed candidate plot to a state file that is stored in your home directory, `~/fj_viewer.state`, so that you can easily check what files you have vetted already.

For the candidate plots, the scripts expects the following file naming convention that is standard for MeerTRAP:

`<detection MJD>_DM_<detection DM>_beam_<beam ID>.jpg`, e.g. `59515.5468834399_DM_196.17_beam_51C.jpg`.

The tool assumes JPEG candidate plots by default, but `ImageMagick` supports a variety of other file formats too (PNG, BMP, etc.).