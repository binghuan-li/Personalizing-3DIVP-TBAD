"""
main.py

Driver to generate synthetic 4D flow inlet velocity profile for CFD simulation.

Ref. https://doi.org/10.1016/j.compbiomed.2025.110158

22/07/2025
"""

import os
from inlet_mapping import mapping
from Profile_scaling import scaling
from write_for_solver import write_to_solver

if __name__ == '__main__':
    # unique identifier of the case
    identifier = ''

    # working directory
    working_dir = r''

    # path to the selected synthetic IVP, e.g., './Profiles/Case01'
    source_profile_dir = r''

    # path to the target inlet plane geometry, in format .stl
    target_profile_fn = r''

    # path to the aorta geometry for direction alignment, in format .stl
    aorta_geo = r'';

    # flip the velocity direction. Usually, this is set to True, 
    # but might have to change depending on target plane orientation
    flip_normals = False

    # path to the selected 0D inlet flow waveform, in format .csv
    profile_path = r''

    # time 
    tcycle=0.;


    working_dir = os.path.join(working_dir, identifier)
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
        aorta_geo=aorta_geo, 
        tcycle=tcycle);
