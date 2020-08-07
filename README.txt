# COMP9517_Project
Check the requirements.txt 
To run open jupyter file and in the main function enter the (dataset you want to see, sequence from that dataset)
Note from dataset 1, sequence 1 and 2 were used for developing the neural net hence, test with sequence 3 and 4 

Our program shows the segmentation, tracking, mitosis and analysis of cell motion for different image sequence.
We developed a program that detects cells and shows a bounding box. Cells are tracked from frame to frame and mitosis events are identified. Each cells track is shown as well as the amount of cells found in the image. 
A predicted cell mitosis is shown by changing the colour of the cell bounding box from green to white. We also see this as an output in the terminal. Information about mother and daughter cells are also shown. When a mitosis occurs, each daughter cell is a new cell with a unique id. 
Tracking can follow each cells. If you press "p" on your keyboard, you can choose a cell to analyse.
Analysis of this cell will present information in the terminal as per task 3 requirements, ie speed of cell, total distance travelled, net distance travelled, and the confinement ratio of the cells motion. 

#NOTE: images are stored in cell_Images but as they are large we did not include them in the submission. On the forums it said that this was okay. 
