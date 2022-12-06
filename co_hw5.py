import numpy as np
from functions_hw5 import dataFrameMaker
from functions_hw5 import fourier_n_inverse_transform


#exercise 7.4 part a
dow = np.loadtxt("dow.txt")
dow2 = np.loadtxt("dow2.txt") #has 1024 elements


df1 = dataFrameMaker(dow)


first_percentage = 0.10
#parts b thorough c
b2d = fourier_n_inverse_transform(dow,first_percentage,'Initial_Ten_Percent_Kept',\
        'Dow Values',df1, DCT = False)

print("When we set the first 10% of the elements for the Fourier Transformed Dow array to zero\n"
        "and take the inverse Fourier Transform, the original plot\n"
        "gets smoothed out.")

print('\n')

#part d through e
df2 = dataFrameMaker(dow)
second_percentage = 0.02
pt_de = fourier_n_inverse_transform(dow,second_percentage,'Initial_Two_Percent_Kept',\
        'Dow Values',df2, DCT = False)
print("When we set the first 2% of the elements for the Fourier Transformed Dow array to zero\n"
        "and take the inverse Fourier Transform, the original plot\n"
        "gets smoothed out even more than the previous 10% retention of terms.")

#now exercise 7.6 begins
df3 = dataFrameMaker(dow2)
buldge_FT = fourier_n_inverse_transform(dow2,second_percentage,'Initial 2% Kept',\
        'Dow Values', df3, DCT = False)

df4 = dataFrameMaker(dow2)
rid_buldge_FT = fourier_n_inverse_transform(dow2,second_percentage,'Initial 2% Kept',\
        'Dow Values', df3, DCT = True)
