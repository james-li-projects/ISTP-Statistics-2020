# ISTP-Statistics-2020

README for Running the Code
TL;DR: Download all the files on the GitHub repository and use the “Instructions for running code using GitHub files” to run Step 1 and Step 2 in the code (pre-processing steps) rather than directly running the code in the master_script.R file. Make sure to run these specific instructions before running any of the subsequent analysis in Steps 3 to 6!

Description and background for our files (as well as links to the original larger files on Google Drive):
Step 1 and Step 2: these include the initial processing and annotation steps of our data. These steps are really computationally intensive and our GitHub repository only contains data files that are already processed in order to save run time. Here are Google Drive links to the raw data files in case you would like to run them (EXTREMELY time intense in Step 2): 
https://drive.google.com/file/d/1q08eXj0qvNMuH-LuVP6Fu6mSO-W01B-j/view?usp=sharing 
https://drive.google.com/file/d/1fBdi0DzRu-jTROAQ5_80nuUH-x3NM9Ot/view?usp=sharing
https://drive.google.com/file/d/1fBdi0DzRu-jTROAQ5_80nuUH-x3NM9Ot/view?usp=sharing
The instructions below are essential to run code using only the GitHub files and you shouldn’t need to run any other commands within Steps 1 and 2 besides the ones listed below.

Instructions for running code using GitHub files (run only these lines in order to bypass the extremely time intensive Steps 1 and 2): 
1) Install and run all necessary packages using lines 11-36
2) Load our central gene expression data table for all the patients using lines 66-68
3) Load pairwise cancer cohort comparison data for every single gene using line 102
4) Congrats on finishing Step 1 and 2! Now you can run the code normally and however you want from Step 3 and after!
