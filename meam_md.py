from squid import files, lammps, structures
from pathlib import Path
import os.path as osp
import click
import yaml

BASE_PATH = Path(__file__).parent
xyz_path = BASE_PATH / "structs"
data_path = BASE_PATH / "data"
ff_path = BASE_PATH / "force_fields"
config_path = BASE_PATH / "config"

# PREBASH = "ml cl-lammps\nsource /data/apps/go.sh\n"

DELTA = 0.5


def dump_atoms_data(system, boundary=True):
    """
    Generate a formatted string representation of atom data and the number of atoms within a specified boundary.

    Args:
        system (structures.System): A molecular system containing atoms.
        boundary (bool, optional): Flag indicating whether to consider the boundary when extracting atom indices. Defaults to True.

    Returns:
        tuple: A tuple containing a formatted string representation of atom data and the number of atoms within the boundary.
    """
    max_x = max([a.x for a in system.atoms])
    min_x = min([a.x for a in system.atoms])
    max_y = max([a.y for a in system.atoms])
    max_z = max([a.z for a in system.atoms])
    atoms_idx = []
    atom_subset = []
    if boundary:
        count = 0
        for i, a in enumerate(system.atoms):
            if a.y < max_y - DELTA:
                if a.z < max_z - DELTA:
                    if a.x < max_x - DELTA:
                        atoms_idx.append(i)
                        a.index = count
                        atom_subset.append(a)
                        count = count + 1
    else:
        atoms_idx = list(range(len(system.atoms)))
        atom_subset = system.atoms[:]
    system.atoms = atom_subset
    return (
        "\n".join(
            [
                "%d %d %f %f %f" % (a.index + 1, system.i2t(a.label), a.x, a.y, a.z)
                for i, a in enumerate(system.atoms)
            ]
        ),
        len(atoms_idx),
    )


def write_input_data(system, xyz_low, xyz_high, nt, overlap):
    """
    Generate input data for a LAMMPS simulation by writing atom, atom type, and boundary information to a data file.

    Args:
        system (structures.System): A molecular system containing atoms.
        ff_file (list): A list of two force field file names.
        xyz_low (tuple): A tuple of three floats representing the lower bounds of the simulation box.
        xyz_high (tuple): A tuple of three floats representing the upper bounds of the simulation box.
        job_str (str): A string template for the LAMMPS input file.
        overlap (bool): Flag indicating whether to consider atom overlap when extracting atom indices.

    Returns:
        str: The modified input string for the LAMMPS simulation.
    """
    atoms_data, num_atoms = dump_atoms_data(system, overlap)

    data_file = data_path / (system.name + ".data")
    f = open(data_file, "w")

    f.write("LAMMPS Description\n\n%d atoms\n" % num_atoms)
    f.write("%d atom types\n" % nt)

    f.write("%3.5f %3.5f xlo xhi\n" % (xyz_low[0], xyz_high[0]))
    f.write("%3.5f %3.5f ylo yhi\n" % (xyz_low[1], xyz_high[1]))
    f.write("%3.5f %3.5f zlo zhi\n" % (xyz_low[2], xyz_high[2]))

    f.write("\n\nAtoms\n\n")
    f.write(atoms_data)
    f.write("\n\n")
    f.close()

    return data_file


def job(job_config, struc, box_size):
    """
    Prepare and execute a LAMMPS simulation job.

    Args:
        job_config (dict): A dictionary containing job-specific configuration parameters.
        struc (structures.Molecule): A molecular structure representing the system.
        box_size (tuple): A tuple of two tuples representing the lower and upper bounds of the simulation box.
        atom_dix (dict): A dictionary mapping element names to atom labels.

    Returns:
        lammps.job: A LAMMPS job object.
    """
    job_name = job_config["job_name"]
    job_str = job_config["job_str"]
    job_str = job_str.replace("$FF_ELEM$", job_config["elems"])
    ff_file = job_config["ff"]
    atom_dix = {}
    count = 1
    for atom in struc:
        el = atom.element
        if el not in atom_dix.keys():
            atom_dix[el] = count
            count += 1
        atom.label = atom_dix[el]
    job_str = job_str.replace("$ELEM$", " ".join(list(atom_dix.keys())))
    tmp = structures.Molecule(struc)
    box = structures.System(job_name)
    box.add(tmp)
    box.atom_labels = list(atom_dix.values())

    if job_config.get("data"):
        data_file = data_path / job_config["data"]
    else:
        data_file = write_input_data(
            box,
            box_size[0],
            box_size[1],
            len(atom_dix),
            job_config["overlap"],
        )

    inp_str = job_str.replace("$DATA$", str(data_file.resolve()))
    inp_str = inp_str.replace("$NAME$", str((xyz_path / job_name).resolve()))
    inp_str = inp_str.replace("$THERMO$", " ".join(job_config["output"]))
    if len(ff_file) > 1:
        inp_str = inp_str.replace("$FILE1$", str((ff_path / ff_file[0]).resolve()))
        inp_str = inp_str.replace("$FILE2$", str((ff_path / ff_file[1]).resolve()))
    else:
        inp_str = inp_str.replace("$FILE$", str((ff_path / ff_file[0]).resolve()))

    print("Submitting...")
    resources = job_config["resources"]
    # print(inp_str)
    # raise Exception
    jobj = lammps.job(
        job_name,
        inp_str,
        queue=resources["queue"],
        allocation="pclancy3",
        walltime=resources["walltime"],
        nprocs=resources["cores"],
        pair_coeffs_in_data_file=False,
        prebash=job_config["prebash"],
    )

    return jobj


@click.command()
@click.option("--task", default="md_apart_rel")
@click.option("--name", default=None)
def run(task, name):
    run_base(task, name)


def run_base(task="md_apart_rel", name=None):
    """
    Execute a LAMMPS simulation job based on the provided task and name.

    Args:
        task (str, optional): The task or YAML configuration file name. Defaults to "md_apart_rel".
        name (str, optional): The job name override. Defaults to None.
    """
    # atom_dix = {"Ni": 1, "Ti": 2, "Cr": 3, "Ag": 4}

    systems_path = config_path / "systems/md"
    task = systems_path / (task + ".yaml")

    inputs_path = config_path / "inputs"
    all_sys = list(systems_path.iterdir())
    if task in all_sys:
        yaml_file = task
        with open(yaml_file, "r") as y:
            config = yaml.safe_load(y)
            if name is not None:
                config["job_name"] = name
            struc = files.read_xyz(str(xyz_path / config["xyz"]))
            inp = config["sim"]
            config["job_str"] = yaml.safe_load(
                open(str(inputs_path / (inp + ".yaml")))
            )["job_str"]
            job(
                config,
                struc,
                config["system_bounds"],
            )
    else:
        return


if __name__ == "__main__":
    run()
    # tasks = [
    #     'meam_cr_bulk',
    #     'meam_cr_free1',
    #     'meam_cr_free2',
    #     'meam_cr_unit',
    #     'meam_cr110_free2',
    #     'meam_cr110_free1',
    #     'meam_cr111_free2',
    #     'meam_cr111_free1',
    # ]
    # tasks = [
    #     'eam_niti_bulk',
    #     'eam_niti_free1ni',
    #     'eam_niti_free1ti',
    #     'eam_niti_free2ti',
    #     'eam_niti_free2ni',
    #     'eam_niti_unit',
    #     'eam_niti110_free2',
    #     'eam_niti110_free1',
    #     'eam_niti111ni_free2',
    #     'eam_niti111ti_free2',
    #     'eam_niti111ni_free1',
    #     'eam_niti111ti_free1'
    # ]
    # for t in tasks:
    #     run_base(t)