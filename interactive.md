Lipydomics Interactive Interface
========

Interactive Interface of Lipydomics package offers easy to follow and user-friendly way of using Lipydomics. 
	

Getting Started
--------

- Clone Lipydomics from: https://github.com/dylanhross/lipydomics

- Or `git clone https://github.com/dylanhross/lipydomics.git` from commandline

- Navigate to lipydomics folder and type `python interactive.py` from commandline to start the interface

Functionalities
------------

- **The Dataset**:

	- The user can either make a new lipydomics dataset or load previous one in the beginning of the interface by specifying a path of a file (either .csv or .pickle). 
	
	- `! INFO: Loaded a new Dataset` message is shown upon successful load.
	
	- Example of making a new Dataset
	```
	> python interactive.py
	What would you like to do?
		1. Make a new Dataset
		2. Load a previous Dataset
	> 1
	Please enter the path to the csv file you want to work with.
	> /path/to/the/dataset.csv
	What ESI mode was used for this data? (pos/neg)
	> pos
	! INFO: Loaded a new Dataset from .csv file: "dataset.csv"
	Would you like to automatically assign groups from headers? (y/N)
	> y
	```


- **Managing groups**:

	- Lipydomics allows user to assign groups to specific columns of the dataset. Groups need to be assigned so that the user can view, filter and compute statistics by group.
	
	- `! INFO: Assigned indices: (specified indices) to group: "specified group name"` message is shown upon successful assigning.
	
	- Example of making new groups
	```
	What would you like to do with this Dataset? 
		1.  Manage Groups
		2.  Filter Data
		3.  Manage Statistics
		4.  Make Plots
		5.  Lipid Identification
		6.  Normalize Intensities
		7.  Calibrate Retention Time
		8.  Overview of Dataset
		9.  Export Current Dataset to Spreadsheet
		10. Save Current Dataset to File
		"exit" to quit the interface
	> 1
	Managing groups... What would you like to do?
		1. Assign group
		2. View assigned groups
		3. Get data by group(s)
		"back" to go back
	> 1
	Please provide a name for a group and its indices in order of name > starting index > ending index.
		* group name should not contain spaces
		* indices start at 0
		* example: 'A 1 3'
	> A 0 3
	! INFO: Assigned indices: [0, 1, 2, 3] to group: "A"
	Managing groups... What would you like to do?
		1. Assign group
		2. View assigned groups
		3. Get data by group(s)
		"back" to go back
	```
	
- **Filtering data**:

	- The interactive interface allows user to filter data on given range of M/Z, retention time and ccs or filter on PLS-DA loadings and group correlation values (S-plot). 

	- The user can also make batch-queries with a csv file with a correct format. 

	- The filtered data is available for download.
	
	- `! INFO: Successfully filtered data. Would you like to download the result as csv? (y/N)` Input "y" to download the result as a .csv file
	
	- `! INFO: Successfully downloaded the filtered data` message is shown upon successful download
	
	- [Sample .csv for making batch queries](./sample_batch_query.csv)


- **Computing Statistics**:

	- The user can choose specific groups assigned previously to compute statistic of choosing on those groups.
	
	- Groups need to be assigned before computing any statistic. 
	
	- Available statistics are:

		- Anova-P
		
		- PCA3
		
		- PLS-DA
		
		- Two Group Correlation
		
		- PLS-RA (using an external continuous variable)


- **Plotting**:

	- After computing the necessary statitics, the user can make plots out of the data. 
	
	- The plot will be automatically generated to the user specified path.

	- Available plots are:
	
		- Bar plot feature(s) by group
		
		- Scatter PCA3 Projections by group
		
		- Scatter PLS-DA Projections by group
		
		- S-Plot PLSA-DA and Pearson correlation by group
		
		- Scatter PLS-RA projections by group


- **Lipid Identification**:

	- Allows the user to perform lipid identification at a variety of levels of confidence
	
	- Levels of Confidence:
	
		- `theo_mz` - match on theoretical m/z
		
		- `theo_mz` - match on theoretical m/z and retention time
		
		- `theo_mz_ccs` - match on theoretical m/z and CCS
		
		- `theo_mz_rt_ccs` - match on theoretical m/z, retention time, and CCS
		
		- `meas_mz_ccs` - match on measured m/z and CCS
		
		- `meas_mz_rt_ccs` - match on measured m/z, retention time, and CCS
		
		- `any` - try all criteria (highest confidence first)


- **Normalization**

	- The interactive module allows user to do two types of normalizations:
	
		- Internal Normalization
			
			- User is to specify a feature range. The module then will find maximum value of intensity that matches the specified feature range. All intensities is to be divided with the maximum intensity and those become normalization weights.
			
		- External Normalization
		
			- User is to specify normalization weights in the .txt file each weight in a new line. Make sure that number of weights match number of intensity values in the dataset. 


- **Calibrating Retention Time**

	- The interactive module allows user to calibrate retention time in the given data set. 
	
	- The user can:
		
		- Build new retention time calibration
		
		- View retention time calibration
		
		- Clear current retention time calibration

- **Download as Spreadsheet**

	- The user can specify a path to download excel version of the dataset object.
	
	- The original data, computed statistics, normalized data, identified lipids will be in the file. 


- **Download as Binary**

	- The user can specify a path to download binary version of the python dataset object.
	
	- This file can be loaded at the start of the interactive module.

	
- **Common Errors**

	- `! ERROR: No normalization has been performed yet` This message will show if data has not been normalized yet.
	
	- `! ERROR: Unable to perform statistical analysis, check group names and try again` Make sure that the groups are pre-assigned.
	
	- `! ERROR: Failed to export statistics to file: "(Path name)"` Make sure that path entered is correct.
	
	- `! ERROR: Failed to filter data, please check your groups and try again` Make sure that the groups are pre-assigned.
	
	- `! ERROR: Could not find any matching point` The program failed to find any matching points in the dataset.
	
	- `! ERROR: unrecognized option: "(Option)"` Make sure to choose from given options

	- `! ERROR: Make sure the path specified is correct and the file exists.` message is shown upon failure. Make sure to check the file exits and is in right format specified below


