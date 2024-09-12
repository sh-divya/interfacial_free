from squid import files, lammps, structures
# from meam_md import dump_atoms_data, write_input_data
from meam_md import job
from sys_prep import box_info

from pathlib import Path
from copy import deepcopy
import click
import yaml

BASE_PATH = Path(__file__).parent
xyz_path = BASE_PATH / "structs"
data_path = BASE_PATH / "data"
ff_path = BASE_PATH / "force_fields"
config_path = BASE_PATH / "config"


def ti_jobs(job_config, struc):
    """
    Execute Thermodynamics Integration (TI) jobs using LAMMPS.

    Args:
        job_config (dict): Configuration parameters for the job.
        struc (list): List of frames or structures.
        box_size (list): Box size boundaries.
        atom_dix (dict): Dictionary mapping atom elements to labels.

    Returns:
        lammps.job: A LAMMPS job object for last TI frame.
    """
    job_name_base = job_config["job_name"]
    job_str = job_config["job_str"]
    overlap = job_config["overlap"]
    for f, frame in enumerate(struc):
        job_name = job_name_base + "_" + str(f)
        iconfig = deepcopy(job_config)
        iconfig["job_name"] = job_name
        # atom.label = atom_dix[atom.element]
        _, _, lz = box_info(frame)
        lz = lz + 10
        iconfig["system_bounds"][1][-1] = lz
        job(iconfig, frame, iconfig["system_bounds"])


def lambda_jobs(job_config, struc):
    """
    Execute Lambda Integration (LMDA) jobs using LAMMPS.

    Args:
        job_config (dict): Configuration parameters for the job.
        struc (list): List of structures.
        box_size (list): Box size boundaries.
        atom_dix (dict): Dictionary mapping atom elements to labels.

    Returns:
        lammps.job: A LAMMPS job object.
    """
    raise Exception
    job_name = job_config["job_name"]
    job_str = job_config["job_str"]
    overlap = job_config["overlap"]
    for atom in struc:
        atom.label = atom_dix[atom.element]

    tmp = structures.Molecule(struc)
    box = structures.System(job_name)
    box.add(tmp)
    box.atom_labels = list(atom_dix.values())
    lmp_inp_str = write_input_data(
        box,
        job_config["ff"],
        box_size[0],
        box_size[1],
        job_config["elems"],
        job_str,
        overlap,
    )

    jobj = lammps.job(
        job_name,
        lmp_inp_str,
        queue="defq",
        allocation="pclancy3",
        walltime="24:0:0",
        nprocs=8,
        pair_coeffs_in_data_file=False,
    )


@click.command()
@click.option("--task", default="ti_apart")
@click.option("--name", default=None)
def run(task, name):
    run_base(task, name)

def run_base(task, name=None):
    """
    Execute a LAMMPS simulation job based on the provided task and name.

    Args:
        task (str, optional): The task or YAML configuration file name. Defaults to "md_apart_rel".
        name (str, optional): The job name override. Defaults to None.
    """
    atom_dix = {"Ni": 1, "Ti": 2, "Cr": 3}

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
            if config["pmf_type"] == "ti":
                config["job_str"] = yaml.safe_load(open(inputs_path / f"{config['sim']}.yaml", "r"))[
                    "job_str"
                ]
                ti_jobs(config, struc)
            elif config["pmf_type"] == "lmda":
                config["job_str"] = yaml.safe_load(
                    open(inputs_path / "lmda.yaml", "r")
                )["job_str"]
                lambda_jobs(config, struc)
    else:
        return
    

if __name__ == "__main__":
    run_base("meam_crniti_110_ti")
