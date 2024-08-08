from squid import files, structures
import os.path as osp
from pathlib import Path
import numpy as np
from pymatgen.core import Structure
from pymatgen.transformations.advanced_transformations import SupercellTransformation
from pymatgen.io.xyz import XYZ
from ase.build.surface import fcc111
from ase.io import read, write

struct_path = Path(__file__).parent / "structs"


def combine(name, surf1, surf2, distance, boundary=False):
    surf1 = files.read_xyz(str(struct_path / surf1))

    surf2 = files.read_xyz(str(struct_path / surf2))

    combined_surface = []
    max_x = max([a.x for a in surf1])
    max_y = max([a.y for a in surf2])
    min_z1 = min([atom.z for atom in surf1])
    max_z2 = max([atom.z for atom in surf1])

    if boundary:
        surf1_index = []
        surf2_index = []
        for i, a in enumerate(surf1):
            if a.y < max_y:
                if a.x < max_x:
                    surf1_index.append(i)
        for i, a in enumerate(surf2):
            if a.y < max_y:
                if a.x < max_x:
                    surf2_index.append(i)
    else:
        surf1_index = [i for i in range(len(surf1))]
        surf2_index = [i for i in range(len(surf2))]

    if isinstance(distance, int):
        distance = [distance, 0.0, 0.0]

    for atom in surf2:
        atom.translate([distance[0], distance[1], distance[2]])

    # print("Surf1 bounds")
    # box_info(surf1)
    # print("Surf2 bounds")
    # box_info(surf2)
    combined_surface += [atom for a, atom in enumerate(surf1) if a in surf1_index]
    combined_surface += [atom for a, atom in enumerate(surf2) if a in surf2_index]
    files.write_xyz(combined_surface, str(struct_path / name))


def box_info(system):
    if isinstance(system, str):
        system = files.read_xyz(str(struct_path / system))

    max_x = max([a.x for a in system])
    min_x = min([a.x for a in system])

    max_y = max([a.y for a in system])
    min_y = min([a.y for a in system])

    max_z = max([a.z for a in system])
    min_z = min([a.z for a in system])

    # print("Lower bounds for XYZ")
    # print(min_x, min_y, min_z)

    # print("Box Dims")
    # print(max_x - min_x, max_y - min_y, max_z - min_z)
    return (max_x - min_x, max_y - min_y, max_z - min_z)


def overlap(system):
    system = files.read_xyz(str(struct_path / system))

    coords = [atom.flatten() for atom in system]
    print(np.array(coords).shape)
    u, c = np.unique(np.array(coords), return_counts=True, axis=0)
    dup = u[c > 1]

    print(dup)


def supercell(cif, replicas=[3, 3, 3]):
    cif_path = struct_path / cif
    unit = Structure.from_file(str(cif_path))
    matrix = []
    for i, r in enumerate(replicas):
        tmp = [0 for _ in range(3)]
        tmp[i] = r
        matrix.append(tmp)
    transform = SupercellTransformation(matrix)
    cell = transform.apply_transformation(unit)
    cif = cif.split(".")[0]
    replicas = "".join(map(str, replicas))
    xyz = XYZ(cell)
    xyz.write_file(filename=str(struct_path / f"{cif}_{replicas}.xyz"))


def fcc_surface(prefix, element, latLen, size=[3, 3, 3]):
    slab = fcc111(element, a=latLen, size=size, vacuum=None)
    size = "".join(map(str, size))
    write(struct_path / f"{prefix}_{size}.xyz", slab, "xyz")
    return slab


def unwrapz(atoms, box_dims):
    for i, a1 in enumerate(atoms):
        d1 = a1.flatten()
        for j, a2 in enumerate(atoms[i + 1 :]):
            d2 = a2.flatten()
            delta = (d1 - d2).tolist()
            if abs(delta[2]) > box_dims[2] / 2:
                d = delta[2]
                lz = box_dims[2]
                if d > 0:
                    d2[2] = d2[2] + lz
                if d < 0:
                    d2[2] = d2[2] - lz
                atoms[j + i + 1].set_position(d2)

    return atoms


def add_strain(atoms, box_dims, strain, direction):
    pass


if __name__ == "__main__":
    # fcc_surface("ni_eam", "Ni", 3.5181211, [7, 7, 9])
    # supercell("ag_eam.cif", [12, 12, 11])
    # supercell("ni_eam.cif", [14, 14, 13])
    # box_info("ag_eam_121211.xyz")
    # box_info("ni_eam_141413.xyz")
    # box_info("ni_unit.xyz")
    # box_info("ag_unit.xyz")
    # box_info("Cr_110_128atoms.xyz")
    # box_info("ni111.xyz")
    # combine("agni111_v2.xyz", "ag111_555.xyz", "ni111_666.xyz", [0, 0, 39])
    print(box_info("agni111_v2.xyz"))
    # combine("apart_new.xyz", "cr110_9rel.xyz", "niti110_9rel.xyz", [20, 0, 0])
    # combine("interface_new.xyz", "Cr_110_9layers.xyz", "niti_110_9layers.xyz", [2, 0, 0])
    # combine("interface0_scf.xyz", "cr110_qe_rel.xyz", "niti110_qe_rel.xyz", [8.26, 2, 0])
    # box_info(files.read_xyz(str(struct_path / "niti_110_9layers.xyz")))
    # box_info(files.read_xyz(str(struct_path / "apart_smd_init.xyz")))
    # box_info(files.read_xyz(str(struct_path / "interface0_scf.xyz")))
