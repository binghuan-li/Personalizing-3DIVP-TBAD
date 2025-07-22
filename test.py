"""
main.py
"""
import os
from inlet_mapping import mapping
from Profile_scaling import scaling
from write_for_solver import write_to_solver

if __name__ == '__main__':

    identifier = 'ZS_post1_IVP'
    working_dir = r'D:/ZhongShan_Data/post_setup'
    working_dir = os.path.join(working_dir, identifier)

    source_profile_dir = r'./Profiles/Case01'
    target_profile_fn = r'D:/ZhongShan_Data/post_mesh/inlet_plane/in.stl'  
    flip_normals = False
    profile_path = r'./Inlet Flow Waveform/SV_60/Waveform_60_01.csv'
    tcycle=0.571;
    
    mapping(
        output_dir=working_dir, 
        source_profile_dir=source_profile_dir, inlet_plane=target_profile_fn, 
        flip_normals=flip_normals, 
        identifier=identifier)
    
    # profile_path = r'D:\ZhongShan_Data\post_setup\ZS_post1_IVP\ZS_post1_IVP_flowrate.csv'

    scaling(
        step=20, 
        working_dir=working_dir, 
        profile_path=profile_path, 
        output_dir=working_dir, 
        patient_specific=True, 
        identifier=identifier)

    write_to_solver(
        solver='cfx',
        identifier=identifier, 
        working_dir=working_dir, 
        aorta_geo=r'D:\ZhongShan_Data\scaled_for_flow.stl', 
        tcycle=tcycle);
