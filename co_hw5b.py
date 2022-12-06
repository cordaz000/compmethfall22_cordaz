import numpy as np
import argparse
from functions_hw5 import dataFrameMaker
from functions_hw5 import fourier_n_inverse_transform

dow = np.loadtxt("dow.txt")
dow2 = np.loadtxt("dow2.txt") #has 1024 elements

first_percentage = 0.10
second_percentage = 0.02

if __name__ == "__main__":
        #Use a dictionary to map how you'll use the function fourier_n_inverse_transform
        y_label = 'Dow Values'
        function_mapper = {"DFT_10_DOW": fourier_n_inverse_transform(dow,first_percentage,\
        'Discreet FT w/ 1st 10% values kept',y_label,dataFrameMaker(dow),DCT = False),

                        "DFT_02_DOW": fourier_n_inverse_transform(dow,second_percentage,\
        'Discreet FT w/ 1st 2% values kept',y_label,dataFrameMaker(dow),DCT = False, plotting=False),

                        "DFT_10_DOW2": fourier_n_inverse_transform(dow2,first_percentage,\
        'Discreet FT w/ 1st 10% values kept',y_label,dataFrameMaker(dow2),DCT = False),
        
                        "DFT_02_DOW2": fourier_n_inverse_transform(dow2,second_percentage,\
        'Discreet FT w/ 1st 2% values kept',y_label,dataFrameMaker(dow2),DCT = False),
        
                        "DCT_10_DOW": fourier_n_inverse_transform(dow,first_percentage,\
        'Discreet cosine transform w/ 1st 10% values kept',y_label,dataFrameMaker(dow),DCT = True),
        
                        "DCT_02_DOW": fourier_n_inverse_transform(dow,second_percentage,\
        'Discreet cosine transfor w/ 1st 2% values kept',y_label,dataFrameMaker(dow),DCT = True),

                        "DCT_10_DOW2": fourier_n_inverse_transform(dow2,first_percentage,\
        'Discreet cosine transfor w/ 1st 10% values kept',y_label,dataFrameMaker(dow2),DCT = True),

                        "DCT_02_DOW2": fourier_n_inverse_transform(dow2,second_percentage,\
        'Discreet cosine transfor w/ 1st 2% values kept',y_label,dataFrameMaker(dow2),DCT = True)
                         }

        script_does = 'This portion of the script asks a user to select a case to run.'
        parser = argparse.ArgumentParser(description=script_does)
        parser.add_argument("scriptName", type= str)
        parser.add_argument("command", choices=function_mapper.keys())
        args = parser.parse_args()
        func = function[args.command]
