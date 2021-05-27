import re
import os
import sys
import subprocess
import mysql.connector

class SequenceRetriever:
    def __init__(self):
        self.dir_path = os.path.dirname(os.path.realpath(__file__)) + '/SRA-Numbers' #directory
        self.sra_toolkit_path = './sratoolkit.2.11.0-centos_linux64/bin'
        self.output_dir = "Output"



    def download_data(self, logFileName: str, verifyData = False):
        '''
        Downloads the sequence data from the correspounding SRA accession numbers using the SRA Toolkit

        :param logFileName: name of the file containing the SRA numbers
        :param verifyData: boolean to whether if the SRA accession number format needs to be verified 
        '''
        sra_data = []
        incorrect_sra = []

        with open(os.path.join(self.dir_path, logFileName), "r", encoding="utf-8") as file_ptr:
            sra_data = file_ptr.readlines()

        # verifies data if verifyData is true
        if verifyData == True:
            verified_input = self.verify_sra_input(sra_data)
            sra_data = verified_input[0]
            incorrect_sra = verified_input[1]

        for sra_num in sra_data:
            sra_num = sra_num.strip().replace('\n', '')

            #print(f"Downloading SRA Sequence {sra_num}...")
            subprocess.run(f"{self.sra_toolkit_path}/prefetch {sra_num}", shell=True)
            print(f"Converting {sra_num}.sra to {sra_num}.fastq...")
            subprocess.run(f"{self.sra_toolkit_path}/fastq-dump --split-files --fasta 60 --outdir {sra_num} {sra_num}", shell=True)
        
        print("\nStarting Data Validation...")
        self.verify_sra_data(sra_data)

        print("\nDownload Completed")

        # if verifyData == True:
        #     print("The data to the following SRA numbers has not been downloaded:")
        #     for incorrect_num in incorrect_sra:
        #         print(incorrect_num)



    def verify_sra_input(self, sra_input: list) -> list:
        '''
        Verifies the input SRA accession numbers to ensure their format validity.

        :param sra_data: a list of SRA numbers to be verified
        :return: a list containing two lists. List one contains the correct SRA numbers. List two contains the incorrect SRA numbers.
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
                print(f"SRA number {sra_num} is incorrect.")

        return [correct_data, incorrect_data]



    def verify_sra_data(self, sra_input_list: list) -> bool:
        '''
        Function used to validate the downloaded data using the SRA Toolkit.

        :param sra_input_list: a list of SRA numbers to be validated
        :return: a boolean that represent whether if any data failed to be validated
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
                    if len(line) > 0:
                        print(f'Data {sra_num} validation failed.\nError: {line}')
                    else:
                        print(f'Data {sra_num} validation failed.\nError: Download Failed.')
                    validation = False
            
        return validation



class MySQLRetriever:
    def __init__(self):
        self.query = 'SELECT sra_num FROM sra_data'

    def retrieve_sra_num(self) -> list:
        '''
        Retrieves SRA accession numbers from the MySQL database

        :return: a list containing the SRA accession numbers retrieved
        :rtype: list
        '''
        database = mysql.connector.connect(
            host='localhost',
            user='root',
            passwd='000000',
            database='sra_test'
        )

        mycursor = database.cursor()
        try:
            mycursor.execute(self.query)

        except mysql.connector.Error as e:
            print("Error reading data from MySQL table: ", e)

        sra_list = []
        for sra in mycursor:
            sra_list.append(sra[0])

        return sra_list
    


    def store_sra_num(self, log_file_name: str, sra_list: list):
        '''
        Stores the retrieved SRA accession numbers in the designated txt file

        :param log_file_name: name for the txt file storing the SRA accession numbers
        :param sra_list: a list of SRA accession numbers
        '''
        seqRetr = SequenceRetriever()
        file_path = seqRetr.dir_path
        
        with open(os.path.join(file_path, log_file_name), "w", encoding="utf-8") as file_ptr:
            for sra in sra_list:
                file_ptr.write(f"{sra}\n")


# mysqlRetr = MySQLRetriever()
# mysqlRetr.store_sra_num("test.txt", mysqlRetr.retrieve_sra_num())

# instance = SequenceRetriever()
# instance.download_data("SRA-Log", verifyData = True)


#  subprocess.run(f"{self.sra_toolkit_path}/fastq-dump --split-files --fasta 60 --outdir {self.output_dir} {sra_num}", shell=True)