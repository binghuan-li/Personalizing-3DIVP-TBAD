"""
main.py
"""


from inlet_mapping import mapping
from Profile_scaling import scaling
from write_for_solver import write_to_solver

import os
import argparse

if __name__ == '__main__':

    # parser = argparse.ArgumentParser(description='Process some integers.')
    # parser.add_argument('num1', type=int, help='The first number to add.')
    # parser.add_argument('num2', type=int, help='The second number to add.')
    # parser.add_argument("-d", "--debug", action="store_true", help="set debug mode");
    # args = parser.parse_args();


    identifier = 'patientxx'
    working_dir = r'C:\Users\lbing\Desktop\test_IVP'  # path for saving resampled .vtp files
    source_profile_dir = r'./Profiles/Case01'  # path to the selected synthetic IVP folder
    target_profile_fn = r'C:\Users\lbing\Desktop\test_IVP\TF18mo.stl'  #  path to the stl file of the target plane
    flip_normals = True # usually set to True, but might have to change depending on target plane orientation
    profile_path = r'.\Inlet Flow Waveform\SV_60\Waveform_60_01.csv'
    tcycle=0.83;

    working_dir = os.path.join(working_dir, identifier)
    mapping(output_dir=working_dir, source_profile_dir=source_profile_dir, inlet_plane=target_profile_fn, flip_normals=True, identifier=identifier)
    scaling(step=20, mapping_dir=working_dir, profile_path=profile_path, output_dir=working_dir, patient_specific=True, identifier=identifier)

    write_to_solver(solver='cfx', identifier=identifier, profile_dir=working_dir, inlet_plane=target_profile_fn, tcycle=0.983);
