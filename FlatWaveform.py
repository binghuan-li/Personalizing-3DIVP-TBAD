'''
ECG tunned flow waveform
'''


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

# path to generic waveform
GenDir = r'./Inlet Flow Waveform/generic-waveform.xlsx' 

# path for saving resampled .vtp files
identifier = 'ZS_post1_IVP'
working_dir = r'D:/ZhongShan_Data/post_setup'  
OutDir = os.path.join(working_dir, identifier) # path to save Fourier coefficients for flat IVP

Name = identifier # name of the patient

filename = Name+'_coefficient.xlsx'
filename_csv = Name+'_coefficient_CFX.csv'
filename2 = Name+'_flowrate.csv'
filename3 = Name+ '_mean_flowrate.csv'

filepath = os.path.join(OutDir, filename)
filepath2 = os.path.join(OutDir, filename2)
filepath3= os.path.join(OutDir, filename3)
filepath4 = os.path.join(OutDir, filename_csv)
os.makedirs(OutDir, exist_ok=True)

# patient-specific values
SV_ECO = 72 # round to integer
RCSD_ECG = 1. #round to 1 decimal places
# Tcycle = 0.70 # cardiac cycle period
Tcycle=0.571;

Frequency_ratio = 0.62 ############ adjust this value to match RCSD, up to 4 decimal places!!!!
Scaling_factor = 5.5 ############### adjust this value to mathc SV (ml)


siz = 301
generic = pd.read_excel(GenDir,engine='openpyxl', header = None, names=['Coefficient', '=', 'Value', 'unit1', 'unit2'])
all_coeff = generic['Value']

# Extract coefficients
Fcoeff = []
Fcoeff[:31]= all_coeff[:31]
Fcoeff[32:62] = all_coeff[32:62]
Fcoeff = np.array(Fcoeff)

# generic waveform
v0 = 2 * np.pi / Tcycle
t = np.linspace(0,Tcycle,siz)
Flow = 0
for n in range(31):
  cos_term = Fcoeff[n] * np.cos(n * v0 * t)
  sin_term = Fcoeff[n+30] * np.sin(n * v0 * t)
  Flow += cos_term + sin_term

# Tuning RCSD
v0_new = Frequency_ratio * 2 * np.pi / Tcycle
Flow_scaled = 0
for n in range(31):
  cos_term_scaled = Fcoeff[n] * np.cos(n * v0_new * t)
  sin_term_scaled = Fcoeff[n+30] * np.sin(n * v0_new * t)
  Flow_scaled += cos_term_scaled + sin_term_scaled

# Measure RCSD
Flow_copy = Flow_scaled.copy()
zero_crossings = np.diff(np.sign(Flow_copy)) != 0
zero_location = np.where(zero_crossings == True)[0]  # np.where returns a tuple, [0] gets the array
index_max = np.argmax(Flow_copy)
index_min = np.argmin(Flow_copy)

for i in range(len(zero_location)):
  if zero_location[i] > index_max and zero_location[i] < index_min:
    index_end_systole = zero_location[i]
if not index_end_systole:
  print('Please check the flow plot!')

RCSD = t[index_end_systole] / (t[-1] - t[index_end_systole])  # calculate the systole/diastole ratio
RCSD = np.round(RCSD, decimals=1)

# Tuning RCSD and SV
if np.isclose(RCSD - RCSD_ECG, 0):
    ## Measuring stroke volume
    ###  FFT
    x1 = Flow_scaled[0::3]
    N = len(x1) // 2;
    X = np.fft.fft(x1) / N;
    t_fft = np.linspace(0, Tcycle, len(x1))
    #### Calculate coefficients
    a0 = np.real(X[0])
    number = N // 2 + 1
    an = 2 * np.real(X[1:number])
    bn = -2 * np.imag(X[1:number])
    n = np.arange(len(an))

    #### Truncate the series
    nmax = 10  # use the first 10 terms
    a0 = a0 / nmax
    an = an / nmax
    bn = bn / nmax
    #### Construct Fourier series
    Flow_fft = a0
    for n in range(1, len(an)):
        cos_term_fft = an[n - 1] * np.cos(n * v0 * t_fft)
        sin_term_fft = bn[n - 1] * np.sin(n * v0 * t_fft)
        Flow_fft += cos_term_fft + sin_term_fft

    Flow_fft = Scaling_factor * Flow_fft
    SV_fft = np.zeros(len(t_fft) - 1)
    SV_fft = [(Flow_fft[i + 1] + Flow_fft[i]) / 2 * Tcycle / len(x1) for i in range(len(t_fft) - 1)]
    SV_fft_total = np.round(np.sum(SV_fft) * 10 ** 6)
    #plt.plot(t_fft, Flow_fft * 10 ** 6 * 0.06, label='fft')
    plt.plot(t_fft, Flow_fft, label='inlet flow')
    plt.xlabel('Time / s')
    plt.ylabel('Flowrate (L/min)')

    plt.legend()
    plt.show()
    if np.isclose(SV_fft_total - SV_ECO, 0):
        data0 = {f'a0': a0* Scaling_factor}
        data1 = {f'a{i + 1}': [an[i]* Scaling_factor] for i in range(len(an))}
        data2 = {f'b{i + 1}': [bn[i]* Scaling_factor] for i in range(len(bn))}
        combined_data = {**data0, **data1, **data2}
        df = pd.DataFrame(combined_data)
        df = df.transpose()
        df.columns = ['Value']
        df.to_excel(filepath, index=True, engine='openpyxl')
        df2 = pd.DataFrame(Flow_fft)
        df2.to_csv(filepath2, index=False, header=False)
        mean_flowrate = Flow_fft.mean()
        df3= pd.DataFrame([mean_flowrate])
        df3.to_csv(filepath3, index=False, header=False, mode='a')
        ############# Write to CFX format for copying and pasting ##############
        data0_cfx = {f'a0= {a0 * Scaling_factor} [m^3 s^-1]': ['']}
        data1_cfx = {f'a{i + 1}= {an[i]*Scaling_factor} [m^3 s^-1]': [''] for i in range(len(an))}
        data2_cfx = {f'b{i + 1}= {bn[i]*Scaling_factor} [m^3 s^-1]': [''] for i in range(len(bn))}
        combined_data_cfx = {**data0_cfx, **data1_cfx, **data2_cfx}
        df4 = pd.DataFrame(combined_data_cfx)
        df4 = df4.transpose()
        df4.to_csv(filepath4, index=True)
        print('Finish writing.')
    else:
        message = f"Stroke volume is {SV_fft_total} and ECO measurement is {SV_ECO}. Please adjust scaling factor"
        print(message)


else:
    message = f"RCSD is {np.round(RCSD, decimals=2)} and ECG value is {RCSD_ECG:.2f}. Please adjust Frequency_ratio."
    print(message)

######## Fourier series used in CFX ##################
#     Flow_fft = a0 + an[0]*np.cos(v0*t_fft) + bn[0]*np.sin(v0*t_fft) + an[1]*np.cos(2*v0*t_fft) + bn[1]*np.sin(2*v0*t_fft)  +\
# an[2]*np.cos(3*v0*t_fft) + bn[2]*np.sin(3*v0*t_fft) + an[3]*np.cos(4*v0*t_fft) + bn[3]*np.sin(4*v0*t_fft)  +\
# an[4]*np.cos(5*v0*t_fft) + bn[4]*np.sin(5*v0*t_fft) + an[5]*np.cos(6*v0*t_fft) + bn[5]*np.sin(5*v0*t_fft)  +\
# an[6]*np.cos(7*v0*t_fft) + bn[6]*np.sin(7*v0*t_fft) + an[7]*np.cos(8*v0*t_fft) + bn[7]*np.sin(8*v0*t_fft)  +\
# an[8]*np.cos(9*v0*t_fft) + bn[8]*np.sin(9*v0*t_fft) +\
# an[9]*np.cos(10*v0*t_fft) + bn[9]*np.sin(10*v0*t_fft) + an[10]*np.cos(11*v0*t_fft) + bn[10]*np.sin(11*v0*t_fft)  +\
# an[11]*np.cos(12*v0*t_fft) + bn[11]*np.sin(12*v0*t_fft) + an[12]*np.cos(13*v0*t_fft) + bn[12]*np.sin(13*v0*t_fft)  +\
# an[13]*np.cos(14*v0*t_fft) + bn[13]*np.sin(14*v0*t_fft) + an[14]*np.cos(15*v0*t_fft) + bn[14]*np.sin(15*v0*t_fft)  +\
# an[15]*np.cos(16*v0*t_fft) + bn[15]*np.sin(16*v0*t_fft) + an[16]*np.cos(17*v0*t_fft) + bn[16]*np.sin(17*v0*t_fft)  +\
# an[17]*np.cos(18*v0*t_fft) + bn[17]*np.sin(18*v0*t_fft) + an[18]*np.cos(19*v0*t_fft) + bn[18]*np.sin(19*v0*t_fft)  +\
# an[19]*np.cos(20*v0*t_fft) + bn[19]*np.sin(20*v0*t_fft) + an[20]*np.cos(21*v0*t_fft) + bn[20]*np.sin(21*v0*t_fft)  +\
# an[21]*np.cos(22*v0*t_fft) + bn[21]*np.sin(22*v0*t_fft) + an[22]*np.cos(23*v0*t_fft) + bn[22]*np.sin(23*v0*t_fft)  +\
# an[23]*np.cos(24*v0*t_fft) + bn[23]*np.sin(24*v0*t_fft) + an[24]*np.cos(25*v0*t_fft) + bn[24]*np.sin(25*v0*t_fft)

