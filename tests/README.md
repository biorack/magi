# Some notes for tests
* The full workflow test is designed to assess the validity of the entire workflow, as many interconnected pieces rely on each other. It can catch some corner cases, but has not been catching some of the hdf5-related errors when saving an hdf5 file that came up with Manuel's data, mainly NaN-related. Need to add in some ringers for this, particularly making columns with all sorts of different datatypes. 
* full workflow tester script with parameters suffix is the preferred testing script because it also changes every parameter available to magi
* Obviously, need to make unit tests for all functions.
