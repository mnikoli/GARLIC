# GARLIC Software
## About
## Installation
###### Requirements

The GARLIC software was designed in Python (2.7.x) and MySQL and tested under Ubuntu/Mint environment, but should be possible to install and run it on other OS.

In order to run GARLIC, you need to have installed additional binaries and libraries on your machine:
* R (>3.029)
* Bedtools (with exported path to $PATH variable)
* MySQL server
* Python MySQLdb

In case you are missing some of the above listed you can install them using several commands from the terminal, as listed below:
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

######Setup   
GARLIC scripts can be used directly after download, with a few minor adjustments need applied. 
Once all GARLIC dependencies are met (binaries, packages and libraries), please follow the steps below:

1. download GARLIC.vXXX.zip from (bifacility.uni-koeln.de/viewer.php/GARLIC.v1.0.zip)
2. unzip the file
3. open /path_to_garlic/GARLIC/scripts/MySQL_connect.py (line 6) and change user and pass to match the ones you setup during MySQL server installation. In case you provided only a password, your user name is ```root```.
4. Connect to MySQL and create new database with a command: ```mysql -u root -p ```. You will be asked to provide your mysql password. 
5. Once you have logged in to MySQL, the next step is to create new database called ```test_significance1``` using command: ```create database test_significance1;``` and then type ```exit;```.
6. Last step is to import data (located in folder of the downloaded package, needs to be unzipped first) to the newly created database using a command: mysql -u root -p test_significance1 < garlic_db_dump_xx.xxx.sql

That’s all! If you want to run the software use:  ```python /path_to_garlic/GARLIC/scripts/garlic.py```. Complete list of available commands and options with examples are provided further in text.

##Usage

garlic.py [-h] {viewTableData, testDiseases, testCellTypes, integrateGWAS, generateRegionSNPOverlaps, generateMapSNPOverlaps, analyseCombinations}



viewTableData 	galic viewTableData [-n TABLE_NAME]

Options:
  		-n TABLE_NAME  	Database table name

Description:
Lists the contents of the given table from the database. If table_name is not provided, a list of available table names is printed as output.

testDiseases 	garlic testDiseases [-n REG_MAP_NAME] [-l MIN_GRR_NUM][-t NUM_THREADS] [-s NUM_ITERS] <path>
		
path			path to regulatory map

Options:
-n REG_MAP_NAME  	Given name for the new regulatory map
  		-l MIN_GRR_NUM   	Minimum number of GRRs in disease
  		-t NUM_THREADS   	Number of threads
  		-s NUM_ITERS     	Number of iterations for simulated score

		Description:
		Used to etiologically connect human complex diseases and cell type-specific CRE maps.


testCellTypes 	garlic testCellTypes <did>

  		did    	     		Disease/trait id from database

Description:
		Determines how “relevant” regulatory maps from different cell types are for a given disease/trait did. 
integrateGWAS	garlic integrateGWAS gwas_input ld_input

 		gwas_input 		Path to GWAS file
  		ld_input    		Path to LD file

Description:
Add unpublished GWAS datasets to database.

generateRegionSNPOverlaps	garlic generateRegionSNPOverlaps chr start end

chr         			Chromosome name
  		start       			Beginning of the region
  		end         		End of the region

Description:
Generates overlaps with a given list of diseases or traits or with whole disease/trait DB dataset (if no disease/trait ids are given –default). Overlaps an input regulatory map with a given set of diseases or traits. If no disease or trait ids are given, overlaps are performed with whole DB.


generateMapSNPOverlaps	garlic generateMapSNPOverlaps [-d did] [-m REG_MAP_ID] [-i REG_MAP_PATH] [-n REG_MAP_NAME]

Options:
  		-d did           		One or more disease ids, separated by comma
-m REG_MAP_ID    	Regulatory map id (if already in database)
  		-i REG_MAP_PATH 	Path to a regulatory map file (if not in database)
  		-n REG_MAP_NAME  	Given name for the new regulatory map

Description:
		Reports overlaps between given regulatory map and given disease/trait-associated SNPs. 

analyseCombinations	garlic analyseCombinations [-h] [-c NUM_COMB_EL] [-t NUM_THREADS] [-s SEED] [-n NUM_COMB_TEST] [-p P_VALUE] [-m P_VALUE_IMPROVEMENT] did

  		did			Disease/trait id from database
	
Options:
  		-c NUM_COMB_EL	Maximum number of regulatory maps in combination [2].
  		-t NUM_THREADS	Number of threads.
  		-s SEED			Parameter to limit the number of 'seed' regulatory maps which are use d as a basis to make combinations
  		-n NUM_COMB_TEST	Number of candidate combinations to test [3].
  		-p P_VALUE		Use only those regulatory maps for which input disease has p value greater than a given parameter [0.001].
  		-m P_VALUE_IMPROVE	Magnitude of p value improvement [5].

Description:
		Algorithm tries to identifiy groups of cell types from different cell types with an increased etiological contribution to a given disease/trait.  


