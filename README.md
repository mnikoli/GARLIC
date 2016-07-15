# GARLIC Software
# About
GARLIC (GWAS-based Prediction Toolkit for connecting Diseases and Cell Types) is a user-friendly computational toolkit that enabls users to etiologically connect human diseases/traits with relevant cell types/tissues through our novel method based on the overlaps between disease-associated genetic variants and cis-regulatory elements, without any a priori assumptions.
# Installation
## Requirements

The GARLIC software was designed in Python (2.7.x) and MySQL and tested under Ubuntu/Mint environment, but should be possible to install and run it on other OS.

In order to run GARLIC, you need to have installed additional binaries and libraries on your machine:
* R (>3.029)
* Bedtools (with exported path to $PATH variable)
* MySQL server
* Python MySQLdb

In case you are missing some of the binaries/libraries listed above, you can install them using several commands from the terminal, as shown below:
```
1. sudo apt-get update
2. sudo apt-get upgrade
3. sudo apt-get install <package_name>
```
, where <package_name> should be changed to one of the package names (corresponding to Ubuntu/Mint pacakges) below:
* r-base  (for R)
* bedtools (for Bedtools)
* mysql-server
* mysql-client
* python-mysqldb	

GARLIC also uses non-default Python packages:
* scipy
* numpy
* rpy2
* xlwt

Missing Python packages can be installed using command ```sudo pip install <package_name>```. if you don’t have pip, you can install it with ```sudo apt-get install python-pip```.

##Setup   
GARLIC scripts can be directly used after download, and although no installation is needed, minor adjustments need to be made. 
Once all GARLIC dependencies are met (binaries, packages and libraries), please follow the steps below:

1. download GARLIC.vXXX.zip from (bifacility.uni-koeln.de/viewer.php/GARLIC.v1.0.zip)
2. unzip the file
3. open /path_to_garlic/GARLIC/scripts/MySQL_connect.py (line 6) and change user and pass to match the ones you setup during MySQL server installation. In case you provided only a password, your user name is ```root```.
4. Connect to MySQL and create a new database using the following command: ```mysql -u root -p ``` (for a user ```root```). You will be asked to type in your mysql password. 
5. Once you have logged in to MySQL, the next step is to create new database called ```test_significance1``` using command: ```create database test_significance1;``` and then type ```exit;```.
6. Last step is to download mysql dump file from http://bifacility.uni-koeln.de/GARLIC/garlic_db_dump_06.16.sql.zip, unpack the file and import it to the newly created database using the following command: mysql -u root -p test_significance1 < garlic_db_dump_06.16.sql. 

That’s all! If you want to run the software use:  ```python /path_to_garlic/GARLIC/scripts/garlic.py```. Complete list of available commands and options with examples are provided further in text.

##Usage

```
garlic.py [-h] {viewTableData, testDiseases, testCellTypes, integrateGWAS, generateRegionSNPOverlaps, generateMapSNPOverlaps, analyseCombinations}
```
```
viewTableData 	galic viewTableData [-n TABLE_NAME]
Options:
  		-n TABLE_NAME  		Database table name
```
Description:
Lists the contents of a given table from the database. If table_name is not provided (by default), a list of available table names is printed as output. Otherwise, all (tab-delimited) columns and table content are printed to standard output. It is recommended to add ```> output_file_name.csv``` at the end of command to save the output in a separate file, as tables contain many rows.  
```
testDiseases 	garlic testDiseases [-n REG_MAP_NAME] [-l MIN_GRR_NUM][-t NUM_THREADS] [-s NUM_ITERS] <path>
path			path to cell type or tissue regulatory map
Options:
-n REG_MAP_NAME  	Given name for the new regulatory map
  		-l MIN_GRR_NUM   	Minimum number of GRRs in disease
  		-t NUM_THREADS   	Number of threads
  		-s NUM_ITERS     	Number of iterations for simulated score
```
Description:
		Predicts the etiological connection of all diseases/traits included in the database with the respect to a given cell type/tissue regulatory map.
```
testCellTypes 	garlic testCellTypes <did>
  		did    	     		Disease/trait id from database
```
Description:
		Predicts the etiological connection of all cell types/tissues included in the database with respect to a given disease/trait.
```
integrateGWAS	garlic integrateGWAS gwas_input ld_input
 		gwas_input 		Path to GWAS file
  		ld_input    		Path to linkage-disequilibrium (LD) file
```
Description:
Add unpublished GWAS datasets to database.
```
generateRegionSNPOverlaps	garlic generateRegionSNPOverlaps chr start end
		chr			Chromosome name
  		start     		Beginning of the region
  		end         		End of the region
```
Description:
Provides SNPs associated with a given disease/trait occurring within a genomic region of interest. If no disease/trait ids are specified (by default), SNPs associated with all diseases/traits included in the database are considered.
```
generateMapSNPOverlaps	garlic generateMapSNPOverlaps [-d did] [-m REG_MAP_ID] [-i REG_MAP_PATH] [-n REG_MAP_NAME]
Options:
  		-d did           	One or more disease ids, separated by comma
		-m REG_MAP_ID		Regulatory map id (if already in database)
  		-i REG_MAP_PATH 	Path to a regulatory map file in BED format (if not in database)
  		-n REG_MAP_NAME  	Given name for the new regulatory map
```
Description:
		Provides SNPs associated with a given disease/trait occurring within a regulatory map of interest.
```
analyseCombinations	garlic analyseCombinations [-h] [-c NUM_COMB_EL] [-t NUM_THREADS] [-s SEED] [-n NUM_COMB_TEST] [-p P_VALUE] [-m P_VALUE_IMPROVEMENT] did
  		did			Disease/trait id from database
Options:
  		-c NUM_COMB_EL		Maximum number of regulatory maps in combination [2].
  		-t NUM_THREADS		Number of threads.
  		-s SEED			Parameter to limit the number of 'seed' regulatory maps which are use d as a basis to make combinations
  		-n NUM_COMB_TEST	Number of candidate combinations to test [3].
  		-p P_VALUE		Use only those regulatory maps for which input disease has p value greater than a 		      given parameter [0.001].
  		-m P_VALUE_IMPROVE	Magnitude of p value improvement [5].
```
Description:
		Tries to identifiy two or more cell types/tissues with an increased combined etiological contribution with a given disease/trait.

##Examples

Here are some examples on how GARLIC can be used:

**Example 1. Retrieve all diseases/traits associated with TBX5 locus (extended 1kb in both directions).**

```python garlic.py generateRegionSNPOverlaps chr12 114790735 114849728```


**Exmaple 2. Retrieve SNPs associated with the Congenital Heart Malformation using DHS regulatory map from fetal heart.**

```python garlic.py generateMapSNPOverlaps -m 30 -d 336```

```python garlic.py generateMapSNPOverlaps -i /path_to_GARLIC_folder/DHS_maps/Fetal_Heart_input.csv -d 336```.


**Example 3. Test the effect of paired cell types/tissues on Type 1 diabetes autoantibodies.**

```python garlic.py analyseCombinations -d 747``` 	


**Example 4. Identify etiologically relevant cell types/tissues for a disease Atrial fibrilation.**

```python garlic.py testCellTypes 641``` 		
