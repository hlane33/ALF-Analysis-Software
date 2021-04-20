# ALF-Analysis-Software
A software package, written by Zihao Liu, using MATLAB to view and analyze ALF data from Rutherford Appleton Labs.

How to use this GUI can be found in our paper, this README file is to provide some technical details regarding the GUI.

1.  When using the GUI, please put the 'grain(1).mat' and 'grain(2).mat' files in the same folder of the neutron counts '.nxs' files. The GUI can be put in any folder, but you need to add the folder that contains the '.nxs' files to the path. This can be done using the 'Folder' or 'Load' button in each tab to choose add the folder to path.

2.  After running the second tab, files called 'Range_record.mat' and 'totalcounts.mat' will be generated.

    'Range_record.mat' has 5 columns:

    First column: index of files. Second column: maximum of phi. Third column: minimum of phi Forth column: maximum of theta Fifth column: minimum of theta

    and this file records the range of theta and phi of the vector G plotted in each file.

    'totalcounts.mat' has the same length of the grain, it has only one column that contains the neutron counts each point of the grain contains.

    When analysing the plotted figures, the GUI will read this two files, so please do not remove them. However, you can copy and store them seperately for use later as they     contain all the information we can get from the experiments data.

3. In this cases we have some corrupted files, you can comment them out in line 554 of the GUI.

4. To understand the comments regarding how the GUI plot the mosaic angle, please have a look first at 'grainmaker.m' regarding how we created the grains.
