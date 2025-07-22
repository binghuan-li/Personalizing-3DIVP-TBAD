'''
ECG tunned flow waveform, 

Use SI units.
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

    plt.figure();
    plt.plot(t_fft, Flow_fft*60000, label='inlet flow');
    plt.xlabel('Time [s]');
    plt.ylabel('Volumetric Flow Rate [L/min]');
    plt.legend();
    plt.show();

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
        data0_cfx = {f'a0in= {a0 * Scaling_factor} ': ['']}
        data1_cfx = {f'a0{i + 1}in= {an[i]*Scaling_factor}': [''] for i in range(len(an))}
        data2_cfx = {f'b0{i + 1}in= {bn[i]*Scaling_factor}': [''] for i in range(len(bn))}
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
