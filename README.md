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
```, where **<package_name>**

should be changed to one of the packaged names (corresponding to Ubuntu/Mint pacakges) below:
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

Missing Python packages can be installed using command ```sudo pip install <package_name>``` 
(if you don’t have pip, use ```sudo apt-get install python-pip```)

######Setup   
GARLIC scripts can be directly used directly after download, no installation is needed, but minor adjustments need applied. 
Once all GARLIC dependencies (binaries, packages and libraries) are installed:
1. download GARLIC.vXXX.zip from (bifacility.uni-koeln.de/viewer.php/GARLIC.v1.0.zip)
2. unzip the file
3. open /path_to_garlic/GARLIC/scripts/MySQL_connect.py
in line 6: change user and pass to match the ones you setup during MySQL server installation. In case you provided only a password, your user name is ```root```
.
4. Connect to MySQL and create new database using a command: mysql -u root -p and providing you mysql password. After execution, you will have to type your password by hand. 
5. Once you have logged in to MySQL create new database called test_significance1 using command: ‘create database test_significance1;’ and then type ‘exit;’.
6. Last step is to import data (located in folder of the downloaded package, but needs to be uncompressed as well) to the newly created database using a command: mysql -u root -p test_significance1 < garlic_db_dump_xx.xxx.sql
7. And that’s all! If you want to run the software use:  python /path_to_garlic/GARLIC/s
