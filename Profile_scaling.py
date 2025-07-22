import os
import os.path as osp
import numpy as np
from glob import glob
import pandas as pd
import pyvista as pv
pv.set_plot_theme("Document")
import matplotlib.pyplot as plt
import descriptors_utils as dut
from scipy import interpolate
from collections import deque
import utils as ut


def scaling(step:int, 
            working_dir:str, 
            profile_path:str, 
            output_dir:str=None, 
            patient_specific:bool=True, 
            identifier:str=None):
    
    '''
    Args:
        - step:int, adjust this number [0-29] to ensure the same SDR between profile and synthetic flowwaveform
        - mapping_dir:str, path to mapped .vtp files from inlet_mapping.py
        - profile_path:str, path to the selected inlet flow waveform .csv file
        - output_dir:str, 
        - patient_specific:bool=True, 
        - identifier:str=None

    Return: None
    '''
    print('>> STEP 2/3: Scaling... ', end='')

    time_intp_options = {
        'T4df': 1,   # adjust to the real period.  for patients with 4D-flow, it should be the value from 4D-flow.
                                                    # For patients without 4D-flow, it should be the same value as expected period
        'Tfxd': 1,   #
        'num_frames_fxd': 20}

    Number = step #### Adjust this number [0-29] to ensure the same SDR between profile and synthetic flowwaveform ####

    # mapping directory

    preprocDir = os.path.join(working_dir, 'inlet_mapping') #r'C:\Users\lbing\Desktop\test_IVP' # Path to mapped .vtp files from inlet_mapping.py
    
    # profile path
    csv_path = profile_path  # path to the selected inlet flowwaveform .csv file

    if output_dir is None:
        output_dir = working_dir; 
    outputDir =  os.path.join(output_dir, 'prof_scaling');
    
    # saving_path
    if identifier is not None:
        filename = identifier + '_scaled_mean_flowrate.csv';
    else:
        filename = 'scaled_mean_flow_rate.csv';
    filepath= os.path.join(outputDir, filename)
    os.makedirs(outputDir, exist_ok=True)

    patient_specific_4D = patient_specific   # True: point by point scaling; False: general scaling
    ##---------------------------------------------------------Do Not Change--------------------------------------------------------

    
    synthetic_planes = [pv.read(fn) for fn in sorted(glob(osp.join(preprocDir, '*.vtp')))]  ## select the mapped file folder


    ## read flow waveform from Matlab output
    tuned_flowrate_csv = pd.read_csv(csv_path, header=None)
    tuned_flowrate_csv.rename(columns={0:'Velocity'},inplace=True)
    tuned_velocity = tuned_flowrate_csv['Velocity']
    t_tuned = np.linspace(0,1,len(tuned_velocity))
    t_new = np.linspace(0,1,20)
    tinterp_tuned_velocity = interpolate.interp1d(t_tuned,tuned_velocity,kind = 'cubic', axis=0)(t_new)
    flowrate_tuned = tinterp_tuned_velocity
    SV_tuned =[flowrate_tuned[0]/20]
    for k in range(len(t_new)-1):
        SV_tuned.append((flowrate_tuned[k]+flowrate_tuned[k+1])/2 * 1/20)

    total_SV_tuned = np.sum(SV_tuned) * 1000 * 1000 # unit: ml/s
    tuned_peak = np.argmax(abs(flowrate_tuned)) # used for following shifting
    tuned_min = np.argmin(flowrate_tuned) # used for following shifting


    ## interpolation on systole and diastole
    extended_planes = ut.time_interpolation_extend(synthetic_planes)  # extend to 30 timeframes in case of shorter RCSD (<0.7)
    synthetic_plannes_aligned = ut.time_interpolation(extended_planes[:Number],time_intp_options) # interpolate to the specific number of timeframes to scale the SDR

    # Peak systole alignment
    SV_synthetic_plannes_aligned = dut.compute_flowrate(synthetic_plannes_aligned)['Q(t)']
    synthetic_peak = np.argmax(SV_synthetic_plannes_aligned)
    q= deque(np.arange(len(synthetic_plannes_aligned)))
    circshift = tuned_peak - synthetic_peak
    q.rotate(circshift)
    synthetic_plannes_aligned = [synthetic_plannes_aligned[int(k)] for k in q]


    # Compute SV from the synthetic IVPs
    SV_syn = dut.compute_flowrate(synthetic_plannes_aligned)['Q(t)']
    SV_area = [SV_syn[0]/20]
    for k in range(len(t_new)-1):
        SV_area.append((SV_syn[k] + SV_syn[k+1])/2 * 1/20)

    ratio = np.max(SV_tuned) / np.max(SV_area)

    ## scale
    n=[]

    ######################   general scaling  ############################
    if patient_specific_4D == False:
        for k in range(len(synthetic_plannes_aligned)):
            if k !=0 and k!= len(synthetic_plannes_aligned)-1:
                if  (SV_syn[k]>0 and flowrate_tuned[k] > 0) or (SV_syn[k]<0 and flowrate_tuned[k] < 0):
                    synthetic_plannes_aligned[k]['Velocity'] *= ratio
                    #synthetic_plannes_aligned[k]['Velocity'] *= 1
                    n.append(k)

            else:
                synthetic_plannes_aligned[k]['Velocity'] *= 1
                n.append(k)
            synthetic_plannes_aligned[k].save(osp.join(outputDir, 'scaledProf_{:02d}.vtp'.format(k)))

    ######################   4D-flow scaling  ############################  
    if patient_specific_4D == True:
        for k in range(len(SV_area)):
            if SV_area[k] != 0 and k !=19:
                ratio_test = SV_tuned[k] / SV_area[k]
                synthetic_plannes_aligned[k]['Velocity'] *= ratio_test
                n.append(k)
            synthetic_plannes_aligned[k].save(osp.join(outputDir, identifier+'scaledProf_{:02d}.vtp'.format(k)))

    flow = dut.compute_flowrate(synthetic_plannes_aligned )['Q(t)']
    Q_meanflow = dut.compute_flowrate(synthetic_plannes_aligned )['Q_mean']

    fig = plt.figure();
    ax = fig.add_subplot(111);
    ax.plot(flowrate_tuned, label='tuned flat')
    ax.plot(flow, label='synthetic scaled')
    ax.axvline(0, linestyle ="--")
    ax.axhline(0, linestyle ="--")
    ax.set_ylabel('flowrate [L/min]')
    plt.legend()
    plt.show()

    mean_flowrate = flow.mean()
    print('\tThe mean flowrate of the synthetic flowwaveform is:', mean_flowrate)
    df3= pd.DataFrame([Q_meanflow])
    df3.to_csv(filepath, index=False, header=False, mode='a')
    print('Done!')
    


