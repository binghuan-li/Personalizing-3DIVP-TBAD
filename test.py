"""
main.py
"""

from inlet_mapping import mapping
from Profile_scaling import scaling
from write_for_solver import write_to_solver

import os
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='');

    parser.add_argument('--identifier', type=str, default='TEST_CASE', help='unique identifier of the patient.');

    parser.add_argument('--save_dir', type=str, default='', help='path to the working/saving directory', required=True);

    parser.add_argument('--source_dir', type=str, help='path to the selected synthetic IVP folder', required=True);

    parser.add_argument('--in_plane', type=str, default='', help='path to the stl file of the target plane, in .stl format', required=True);

    parser.add_argument('--tcycle', type=float, help='cardiac cycle.', required=True)
    parser.add_argument('--flip_normals', action='store_true', help='usually set to True, but might have to change depending on target plane orientation.');

    args = parser.parse_args()



    identifier = 'ZS_post1_IVP'
    working_dir = r'D:/ZhongShan_Data/post_setup'
    working_dir = os.path.join(working_dir, identifier)

    source_profile_dir = r'./Profiles/Case01'
    target_profile_fn = r'D:/ZhongShan_Data/post_mesh/inlet_plane/in.stl'  
    flip_normals = True 
    
    working_dir = os.path.join(working_dir, identifier)
    mapping(output_dir=working_dir, source_profile_dir=source_profile_dir, inlet_plane=target_profile_fn, flip_normals=True, identifier=identifier)
    
    # profile_path = r'D:\ZhongShan_Data\post_setup\ZS_post1_IVP\ZS_post1_IVP_flowrate.csv'
        
    profile_path = r'./Inlet Flow Waveform/SV_70/Waveform_70_02.csv'
    tcycle=0.571;

    # working_dir = r'D:\ZhongShan_Data\post_setup\ZS_post1_IVP'
    # working_dir = r'D:\ZhongShan_Data\post_setup'
    scaling(step=20, mapping_dir=working_dir, profile_path=profile_path, output_dir=working_dir, patient_specific=True, identifier=identifier)

    write_to_solver(solver='cfx', identifier=identifier, profile_dir=working_dir, inlet_plane=target_profile_fn, tcycle=tcycle);
