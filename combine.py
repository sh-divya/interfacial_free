from squid import files, structures
import os.path as osp

struct_path = "/scratch16/pclancy3/divya/interfacial_free/structs"

niti_surface = files.read_xyz(osp.join(struct_path, "niti110.xyz"))
max_y = max([a.y for a in niti_surface])
max_z = max([a.z for a in niti_surface])
niti_idx = [i for i in range(len(niti_surface))]
# for i, a in enumerate(niti_surface):
#     if a.y < max_y:
#         if a.z < max_z:
#             niti_idx.append(i)
print(len(niti_surface))
cr_surface = files.read_xyz(osp.join(struct_path, "cr110.xyz"))
max_y = max([a.y for a in cr_surface])
max_z = max([a.z for a in cr_surface])
cr_idx = [i for i in range(len(cr_surface))]
# for i, a in enumerate(cr_surface):
#     if a.y < max_y:
#         if a.z < max_z:
#             cr_idx.append(i)
print(len(cr_surface))
combined_surface = []

min_niti_x = min([atom.x for atom in niti_surface])
max_cr_x = max([atom.x for atom in cr_surface])

for atom in cr_surface:
    atom.translate([2 * min_niti_x - 26, -2.0137, 0.0])

combined_surface += [atom for a, atom in enumerate(niti_surface) if a in niti_idx]
combined_surface += [atom for a, atom in enumerate(cr_surface) if a in cr_idx]
files.write_xyz(combined_surface, osp.join(struct_path, "crniti_apart.xyz"))
