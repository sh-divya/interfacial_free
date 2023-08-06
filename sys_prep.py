from squid import files, structures
import os.path as osp
from pathlib import Path
import numpy as np


struct_path = Path(__file__).parent / "structs"


def combine(name, surf1, surf2, distance, boundary=False):
    surf1 = files.read_xyz(str(struct_path / surf1))

    surf2 = files.read_xyz(str(struct_path / surf2))

    combined_surface = []
    max_y = max([a.y for a in surf1])
    max_z = max([a.z for a in surf2])
    min_x1 = min([atom.x for atom in surf1])
    max_x2 = max([atom.x for atom in surf1])

    if boundary:
        surf1_index = []
        surf2_index = []
        for i, a in enumerate(surf1):
            if a.y < max_y:
                if a.z < max_z:
                    surf1_index.append(i)
        for i, a in enumerate(surf2):
            if a.y < max_y:
                if a.z < max_z:
                    surf2_index.append(i)
    else:
        surf1_index = [i for i in range(len(surf1))]
        surf2_index = [i for i in range(len(surf2))]

    if isinstance(distance, int):
        distance = [distance, 0.0, 0.0]

    for atom in surf2:
        atom.translate([2 * min_x1 - distance[0], distance[1], distance[2]])

    combined_surface += [atom for a, atom in enumerate(surf1) if a in surf1_index]
    combined_surface += [atom for a, atom in enumerate(surf2) if a in surf2_index]
    files.write_xyz(combined_surface, str(struct_path / name))


def box_info(system):
    system = files.read_xyz(str(struct_path / system))

    max_x = max([a.x for a in system])
    min_x = min([a.x for a in system])

    max_y = max([a.y for a in system])
    min_y = min([a.y for a in system])

    max_z = max([a.z for a in system])
    min_z = min([a.z for a in system])

    print("Lower bounds for XYZ")
    print(min_x, min_y, min_z)

    print("Upper bounds for XYZ")
    print(max_x, max_y, max_z)


def overlap(system):
    system = files.read_xyz(str(struct_path / system))

    coords = [atom.flatten() for atom in system]
    print(np.array(coords).shape)
    u, c = np.unique(np.array(coords), return_counts=True, axis=0)
    dup = u[c > 1]

    print(dup)


if __name__ == "__main__":
    # box_info("ni100.xyz")
    # box_info("ni110.xyz")
    # box_info("Cr_110_128atoms.xyz")
    box_info("ni111.xyz")
