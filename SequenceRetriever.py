# The following things needs to be edited when used:
# 1. SRA Toolkit directory path 
# 2. Output directory path/name
# 3. MySQL database credentials
# 4. SRA accession number retrieval query string

import re
import os
import sys
import datetime
import subprocess
import pymysql

class SequenceRetriever:
    '''
    Class responsible for retrieving the sequence data and confirming the validity of the .sra files
    downloaded. Both the file conversion (.sra to .fastq) and the validation process use SRA Toolkit's 
    VDB module. 
    '''

    def __init__(self, SRARetriever):
        self.dir_path = os.path.dirname(os.path.realpath(__file__)) + '/SRA-Numbers' # log file path
        self.sra_toolkit_path = './sratoolkit.2.11.0-centos_linux64/bin' # SRA Toolkit path

        self.output_dir = "Output" # sequence data output directory name
        self.errorLogFilename = "Error-Log"

        self.SRARetr = SRARetriever


    def download_data(self, logFileName: str, verifyData = False):
        '''
        Download the sequence data from the correspounding SRA accession numbers using the SRA Toolkit.

        :param logFileName: name of the file containing the SRA numbers
        :param verifyData: boolean to whether if the SRA accession number format will be verified 
        '''
        sra_data = []
        incorrect_sra = []

        with open(os.path.join(self.dir_path, logFileName), "r", encoding="utf-8") as file_ptr:
            sra_data = file_ptr.readlines()

        # verifies data if verifyData is true
        if verifyData == True:
            sra_data = self.verify_sra_input(sra_data)

        for sra_num in sra_data:
            sra_num = sra_num.strip().replace('\n', '')

            # download the sequence data using the SRA Toolkit
            subprocess.run(f"{self.sra_toolkit_path}/prefetch {sra_num}", shell=True)
            print(f"Converting {sra_num}.sra to {sra_num}.fastq...")
            subprocess.run(f"{self.sra_toolkit_path}/fastq-dump --split-3 --fasta 60 --outdir {sra_num} {sra_num}", shell=True)
        
        print("\nStarting Data Validation...")
        self.verify_sra_data(sra_data)

        print("\nDownload Completed")


    def verify_sra_input(self, sra_input: list) -> list:
        '''
        Verifies the input SRA accession numbers to ensure that their format is valid. SRA numbers with invalid format will be
        removed from the list of correct SRA numbers. Errors will be logged.

        :param sra_data: a list of SRA numbers to be verified
        :return: a list containing the correctly formatted SRA numbers
        :rtype: list
        '''
        correct_data = []
        incorrect_data = []

        for sra_num in sra_input:
            sra_num = sra_num.strip().replace('\n', '')
            
            if re.match("^[a-zA-Z]{3}[0-9]+$", sra_num): # i.e. SRR1568808
                correct_data.append(sra_num)
            else:
                incorrect_data.append(sra_num)
                errorMess = "Incorrect Formatting"
                print(f"SRA number {sra_num} is incorrect.")

                self.log_error(sra_num, errorMess)

        return correct_data


    def verify_sra_data(self, sra_input_list: list) -> bool:
        '''
        Validates the downloaded .sra files using the SRA Toolkit. Errors will be logged.

        :param sra_input_list: a list of SRA numbers to be validated
        :return: false if any data failed to be validated and true otherwise
        :rtype: boolean
        '''
        validation = True
        
        for sra_num in sra_input_list:
            result = subprocess.getoutput(f'{self.sra_toolkit_path}/vdb-validate {sra_num}')
            result = result.split('\n')

            for line in result:
                if re.search("[a-zA-Z0-9] ok$", line, re.IGNORECASE) != None \
                    or re.search(" reads$", line, re.IGNORECASE) != None \
                    or re.search("is consistent$", line, re.IGNORECASE) != None:
                    continue
                else:
                    errorMess = ""
                    if len(line) > 0:
                        errorMess = line
                        print(f'Data {sra_num} validation failed.\nError: {line}')
                    else:
                        errorMess = "Download Failed"
                        print(f'Data {sra_num} validation failed.\nError: Download Failed.')

                    validation = False
                    self.log_error(sra_num, errorMess)
            
        return validation


    def log_error(self, sra_num: str, error_message: str):
        '''
        Log errors to the correspounding file. 
        Error log format: error_time|sra_num|project_id|attribute_id|sample_id|file_id

        :param sra_num: the SRA accession number where the error originated from
        '''
        with open(os.path.join(self.dir_path, self.errorLogFilename), "a", encoding="utf-8") as file_ptr: 
            errorTime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

            file_ptr.seek(0)
            file_ptr.write(f"{errorTime}|{sra_num}|")
            for data in self.SRARetr.SRAInfo[sra_num]:
                file_ptr.write(f"{data}|")
            file_ptr.write(f'{error_message}\n')



class SRARetriever:
    '''
    Class responsible for the retrieval and storage of SRA accession numbers from the MySQL database.
    '''

    def __init__(self):
        self.query = 'select * from sample_metadata where ATT_ID = 234;'
        self.SRAInfo = {}


    def retrieve_sra_num(self):
        '''
        Retrieves SRA accession numbers from the MySQL database

        :return: a list containing the SRA accession numbers retrieved
        :rtype: list
        '''

        database = pymysql.connect(
            host='localhost',
            user='root',
            passwd='password',
            database='sra_test'
        )

        mycursor = database.cursor()
        try:
            mycursor.execute(self.query)

        except pymysql.Error as e:
            print("Error reading data from MySQL table: ", e)

        for sra in mycursor:
            try:
                self.SRAInfo[sra[4]] = [sra[0], sra[1], sra[2], sra[3]]
            except Exception as error:
                print(f"SRA info hash table write error for {sra[4]}: {error}")


    def store_sra_num(self, log_file_name: str):
        '''
        Stores the retrieved SRA accession numbers in the designated log file

        :param log_file_name: name for the txt file storing the SRA accession numbers
        :param sra_list: a list of SRA accession numbers
        '''
        seqRetr = SequenceRetriever(self)
        file_path = seqRetr.dir_path
        
        with open(os.path.join(file_path, log_file_name), "w", encoding="utf-8") as file_ptr:
            for sra, sra_info in self.SRAInfo.items():
                file_ptr.write(f"{sra}\n")
            print(f"{len(self.SRAInfo.items())} SRA accession number(s) has been successfully written into the log file.")


# Driver Code
if __name__ == "__main__":
    SRARetr = SRARetriever()
    SRARetr.retrieve_sra_num()
    SRARetr.store_sra_num("SRA-Log")

    seqRetr = SequenceRetriever(SRARetr)
    seqRetr.download_data("SRA-Log", verifyData = True)
    

