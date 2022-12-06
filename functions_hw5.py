#companion to hw5

import numpy as np
import datetime
import pandas as pd
import matplotlib.pyplot as plt


def dataFrameMaker(inputarray):
    """
    Return DataFrame for plotting FT and IVTs
    """
    #set up index for DataFrame
    start  = datetime.datetime(2006, 12,25)
    t_indexes = pd.date_range(start, periods = len(inputarray), freq = 'D')

    d = {"Dow": pd.Series(inputarray, index = t_indexes)}
    return pd.DataFrame(d)


def fourier_n_inverse_transform(inputarray,percent,labelLegend,ytitle,\
    DataFrame,DCT):
    """
    Provided an input file with dow.txt info and a desired
    percentage of the first elements to keep, some labels,
    and whether you're doing regular vs cosine FT, return a plot of
    the original data along with the inverse transformed data.
    
    """
    if DCT is not True:
        #compute the 1-d discrete Fourier Transform
        dft_dow = np.fft.rfft(inputarray)

        percenter = int(percent*len(dft_dow))

        #mask elements, keep initial percentage specified in percent
        filterr = np.array([0 if i <= percenter else 1 for i in range(len(dft_dow)) ])
        filtered_initial_percentage = np.copy(dft_dow)

        #where filterr is True, set zero
        filtered_initial_percentage[np.where(filterr)] = 0 

        #find the inverse Fourier Transform of the filtered set
        inv_FT = np.fft.irfft(filtered_initial_percentage)
        
        #plot it
        
        fig,ax = plt.subplots()
        DataFrame[labelLegend] = inv_FT
        DataFrame.plot(ax=ax, ylabel = ytitle)
        plt.show()
    else:
        #Type-II discrete cosine transform (DCT) of real data y
        N = len(inputarray)
        y2 = np.empty(2*N,float)
        y2[:N] = inputarray[:]
        y2[N:] = inputarray[::-1]

        c = np.fft.rfft(y2)
        phi = np.exp(-1j*np.pi*np.arange(N)/(2*N))
        DCT = np.real(phi*c[:N])

        percenter = int(percent*len(c))

        #mask elements given specified initial percentage to be kept
        filterr = np.array([0 if i <= percenter else 1 for i in range(len(c)) ])
        filtered_initial_percentage = np.copy(c)
        filtered_initial_percentage[np.where(filterr)] = 0 

        # 1D inverse DCT Type-II
        #Type-II inverse DCT of a
        N = len(filtered_initial_percentage)
        c = np.empty(N+1,complex)

        phi = np.exp(1j*np.pi*np.arange(N)/(2*N))
        c[:N] = phi*filtered_initial_percentage
        c[N] = 0.0
        IDCT = np.delete(np.fft.irfft(c)[:N],50)
        
        
        fig,ax = plt.subplots()
        DataFrame[labelLegend] = IDCT
        DataFrame.plot(ax=ax, ylabel = ytitle)
        plt.show()



