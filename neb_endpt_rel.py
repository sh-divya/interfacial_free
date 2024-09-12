from squid import files, lammps, structures
import os.path as osp
import copy

# Script 2:
# given an initial structure w/o vacancy,
# and the indices of the vacancy atoms + moving atoms
# create the ideal inital frame.
# ALL ATOM INDICES PROVIDED SHOULD
# BE ZERO-INDEXED

# from ideal inital frame, create ideal final frame
# create data files for both frames
# submit geometry optimizations, of init and final frames
# unwrapped coordinates

base_path = '/scratch16/pclancy3/divya/titanium'
data_path = '/scratch16/pclancy3/divya/titanium/data'
xyz_path = '/scratch16/pclancy3/divya/titanium/xyz'
ff_file = '/scratch16/pclancy3/divya/titanium/ti2004.eam'

job_str= '''units metal
pair_style eam/fs

boundary p p p
lattice bcc 3.29194 origin 0 0 0 orient x 1 0 0 orient y 0 1 0 orient z 0 0 1
region cell block 0 $A$ 0 $B$ 0 $C$ units lattice
create_box 1 cell
read_data $DATA$ add append


pair_coeff * * $FILE$ Ti
neighbor 2.0 bin
neigh_modify delay 10 check yes

dump 1 all custom 1 $NAME$.dump id type xu yu zu
dump_modify 1 sort id
dump 2 all xyz 1 $NAME$.xyz
dump_modify 2 element Ti

timestep 0.001
thermo 10
thermo_style custom pe temp

min_style cg 
minimize 1e-25 1e-25 5000 10000
'''

def write_input_data(struc, abc, name):
    data_str = 'LAMMPS Description\n\n$ATOMS$ atoms\n\nAtoms\n\n'
    for a, atom in enumerate(struc):
        xyz = ['{:.7f}'.format(xyz) for xyz in atom.flatten()]
        data_str += str(a + 1) + '\t1\t' + '\t'.join(xyz) + '\n'
    data_str = data_str.replace('$ATOMS$', str(len(struc)))
    inp_str = job_str.replace('$DATA$', osp.join(data_path, name + '.data'))
    inp_str = inp_str.replace('$NAME$', osp.join(xyz_path, name))
    inp_str = inp_str.replace('$FILE$', ff_file)
    for x, y in zip(abc, ['$A$', '$B$', '$C$']):
        inp_str = inp_str.replace(y, str(x))
    with open(osp.join(data_path, name + '.data'), 'w+') as fobj:
        fobj.write(data_str)
    return inp_str

def create_vacancy(initial_file, vacancy_idx, moving_idx):
    init_str = files.read_xyz(initial_file)
    for i, idx in enumerate(moving_idx):
        if idx < vacancy_idx:
            pass
        if idx > vacancy_idx:
            moving_idx[i] -= 1
    removed_atom = init_str.pop(vacancy_idx)

    return init_str, moving_idx

def create_final(str_with_vacancy, new_moving_idx, moving_final_coords):
    for i, idx in enumerate(new_moving_idx):
        str_with_vacancy[idx].set_position(moving_final_coords[i])

    return str_with_vacancy

if __name__ == '__main__':

    init_file = 'ti_555.xyz'
    abc = [5, 15, 5]
    vacancy = 147
    movers = [[148], [196, 148]]
    final_coords = [[(8.22985, 14.81370, 8.22985)],
                    [(8.22985, 14.81370, 8.22985), (9.87582, 13.16780, 9.87582)]]
    for i in range(2):
        init_str, moving_atoms = create_vacancy(init_file, vacancy, movers[i])
        print(moving_atoms)
        # final_str = create_final(copy.deepcopy(init_str), moving_atoms, final_coords[i])
        # name = init_file.split('.')[0]
        # init = name + '_{}.init'.format(i + 1)
        # final = name + '_{}.final'.format(i + 1)
        # # print(init, final)
        # inp1 = write_input_data(init_str, abc, init)
        # inp2 = write_input_data(final_str, abc, final)
        # job1 = lammps.job(init, inp1, queue='defq', allocation='pclancy3',
        #                  walltime='1:0:0', nprocs=1, pair_coeffs_in_data_file=False)
        # job2 = lammps.job(final, inp2, queue='defq', allocation='pclancy3',
        #             walltime='1:0:0', nprocs=1, pair_coeffs_in_data_file=False)
