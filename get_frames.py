from squid import files, lammps
from copy import deepcopy
import os.path as osp
import numpy as np

# Script 3:
# Read relaxed structures of initial and final frames
# create interpolated frames
# recorner on origin
base_path = '/scratch16/pclancy3/divya/titanium'
data_path = '/scratch16/pclancy3/divya/titanium/data'
xyz_path = '/scratch16/pclancy3/divya/titanium/xyz'
ff_file = '/scratch16/pclancy3/divya/titanium/ti2004.eam'
frame_path = '/scratch16/pclancy3/divya/titanium/frames'

def create_frames(initial_frame, final_frame, num):
    zero_coord = -1 * initial_frame[0].flatten()
    inter_frames = [deepcopy(initial_frame) for _ in range(num)]\
    
    starting_coords = [
        atom.flatten()
        for atom in initial_frame
    ]

    final_coords = [
        atom.flatten()
        for atom in final_frame
    ]

    for i, (init, fin) in enumerate(zip(starting_coords, final_coords)):
        inter_coords = np.linspace(init, fin, num)
        for j, coords in enumerate(inter_coords):
            inter_frames[j][i].set_position(list(coords))
            inter_frames[j][i].translate(list(zero_coord))
    
    return inter_frames

def write_frames_data(frames, name):
    frame_name = name + '.image.$NUM$'
    for f, frame in enumerate(frames):
        temp = []
        if f == 0:
            iname = name + '.init.data'
            image = 'LAMMPS Description\n\n$ATOMS$ atoms\n\nAtoms\n\n'
        else:
            iname = frame_name.replace('$NUM$', str(f))
            image = '$ATOMS$\n'
        for a, atom in enumerate(frame):
            xyz = ['{:.7f}'.format(xyz) for xyz in atom.flatten()]
            if f == 0:
                image += str(a + 1) + '\t1\t' + '\t'.join(xyz) + '\n'
            else:
                image += str(a + 1) + '\t' + ' '.join(xyz) + '\n'
        image = image.replace('$ATOMS$', str(len(frame)))
        with open(osp.join(frame_path, iname), 'w+') as fobj:
            fobj.write(image[:-1])
    files.write_xyz(frames, osp.join(frame_path, name + '.xyz'))
    

if __name__ == '__main__':
    pairs = [
        ['ti_555_1.init.dump', 'ti_555_1.final.dump'],
        ['ti_555_2.init.dump', 'ti_555_2.final.dump']
    ]
    nums = [5, 5]
    names = ['ti_555_2nn',
             'ti_555_ex']
    for r, rxn in enumerate(pairs):
        start = lammps.read_dump(osp.join(xyz_path, rxn[0]),
                                 coordinates=['xu', 'yu', 'zu'],
                                 extras=['id', 'type'])[-1]
        end = lammps.read_dump(osp.join(xyz_path, rxn[1]), coordinates=['xu', 'yu', 'zu'])[-1]
        frames = create_frames(start, end, nums[r])
        write_frames_data(frames, names[r])


