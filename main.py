"""
main.py

Driver to generate synthetic 4D flow inlet velocity profile for CFD simulation.
Date: 22/07/2025
"""
import os
from inlet_mapping import mapping
from Profile_scaling import scaling
from write_for_solver import write_to_solver

def main():
    # unique identifier of the case
    identifier = 'ZS_post1_IVP'

    # working directory
    working_dir = r'D:/ZhongShan_Data/post_setup'

    # path to the selected synthetic IVP folder 
    source_profile_dir = r'./Profiles/Case01'

    # path to the target inlet plane geometry, in format .stl
    target_profile_fn = r'D:/ZhongShan_Data/post_mesh/inlet_plane/in.stl'

    # path to the aorta geometry for direction alignment, in format .stl
    aorta_geo = r'D:\ZhongShan_Data\scaled_for_flow.stl'

    # flip the velocity direction. Usually, this is set to True, 
    # but might have to change depending on target plane orientation
    flip_normals = False

    # path to the selected 0D inlet flow waveform, in format .csv
    profile_path = r'./Inlet Flow Waveform/SV_60/Waveform_60_01.csv'

    # [0-29] to ensure the same SDR between profile and synthetic flowwaveform
    step=20

    # time 
    tcycle=0.571

    #-------------------------------------------------------------------------
    working_dir = os.path.join(working_dir, identifier)
    mapping(
        output_dir=working_dir, 
        source_profile_dir=source_profile_dir, inlet_plane=target_profile_fn, 
        flip_normals=flip_normals, 
        identifier=identifier)
    
    while True:
        usr_check = input("Proceed? [Y/N] ").strip().upper()
        if usr_check == 'Y':
            break
        elif usr_check == 'N':
            print("Exited. Please change parameters.")
            return
        else:
            print("Invalid input.")
    
    scaling(
        step=step, 
        working_dir=working_dir, 
        profile_path=profile_path, 
        output_dir=working_dir, 
        patient_specific=True, 
        identifier=identifier)

    while True:
        usr_check = input("Proceed? [Y/N] ").strip().upper()
        if usr_check == 'Y':
            break
        elif usr_check == 'N':
            print("Exited. Please change parameters.")
            return
        else:
            print("Invalid input.")

    write_to_solver(
        solver='cfx',
        identifier=identifier, 
        working_dir=working_dir, 
        aorta_geo=aorta_geo, 
        tcycle=tcycle);

if __name__ == '__main__':
    main();