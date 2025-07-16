import os
import argparse

from inlet_mapping import mapping
from Profile_scaling import scaling
from write_for_solver import write_to_solver

## ======================================================================

parser = argparse.ArgumentParser(description='');
parser.add_argument('--identifier', type=str, default='TEST_CASE', help='unique identifier of the patient.');
parser.add_argument('--save_dir', type=str, help='path to the working/saving directory', required=True);
parser.add_argument('--source_dir', type=str, help='path to the selected synthetic IVP folder', required=True);
parser.add_argument('--in_plane', type=str, help='path to the stl file of the target plane, in .stl format', required=True);
parser.add_argument('--profile_path', type=str, help='path to the standard flowrate waveform, in .csv format', required=True);
parser.add_argument('--tcycle', type=float, help='cardiac cycle.', required=True);
parser.add_argument('--flip_normal', action='store_true', help='usually set to True, but might have to change depending on target plane orientation.');
parser.add_argument('--solver', type=str, default='cfx', help='usually set to True, but might have to change depending on target plane orientation.');

args = parser.parse_args();

## ======================================================================

working_dir = os.path.join(args.save_dir, args.identifier)

mapping(output_dir=working_dir, 
        source_profile_dir=args.source_dir, 
        inlet_plane=args.in_plane, 
        flip_normals=args.flip_normal, 
        identifier=args.identifier)

scaling(step=20, 
        mapping_dir=working_dir, 
        profile_path=args.profile_path, 
        output_dir=working_dir, 
        patient_specific=True, 
        identifier=args.identifier);

write_to_solver(solver=args.solver, 
                identifier=args.identifier, 
                profile_dir=working_dir, 
                inlet_plane=args.in_plane, 
                tcycle=args.tcycle);