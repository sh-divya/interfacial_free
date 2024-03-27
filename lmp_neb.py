from squid import files, structures
from squid.jobs import submit_job
from meam_md import dump_atoms_data, write_input_data
from pathlib import Path
import os.path as osp
import shutil
import click
import yaml
import os
import itertools
import time

BASE_PATH = Path(__file__).parent.resolve()
xyz_path = BASE_PATH / "structs"
data_path = BASE_PATH / "data"
ff_path = BASE_PATH / "force_fields"
config_path = BASE_PATH / "config"
job_path = BASE_PATH / "neb"

PREBASH = "ml cl-lammps\n"


def generate_parameter_combinations(config, avoid_keys=[]):
    # Extract keys and their corresponding lists of values
    keys = [k for k in config.keys() if k not in avoid_keys]
    value_lists = [config[key] for key in keys]

    # Generate all possible combinations of values
    all_combinations = list(itertools.product(*value_lists))

    # Create dictionaries for each combination
    param_dicts = []
    for combo in all_combinations:
        param_dict = {keys[i]: combo[i] for i in range(len(keys))}
        for key in avoid_keys:
            param_dict[key] = config[key]
        param_dicts.append(param_dict)

    return param_dicts


def neb_slurm(job_name, neb_config, inp_str, prebash):
    if job_path.exists():
        pass
    else:
        job_path.mkdir()
    base_cmd = prebash
    add_modules = []
    cores = neb_config["nframes"]
    walltime = neb_config["walltime"]

    for mod in add_modules:
        base_cmd += "ml " + mod + "\n"

    base_cmd += "mpiexec -np $NFRAMES$ lmp_mpi -partition $NFRAMES$x1 -in $INPUT$"
    inp_path = job_path / job_name

    if inp_path.exists():
        shutil.rmtree(inp_path)
    inp_path.mkdir()
    new_inp = inp_path / f"{job_name}.in"
    with open(str(new_inp), "w+") as fobj:
        fobj.write(inp_str)
    sub_cmd = base_cmd.replace("$NFRAMES$", str(neb_config["nframes"]))
    sub_cmd = sub_cmd.replace("$INPUT$", str(new_inp))
    os.chdir(str(inp_path))
    job = submit_job(
        job_name,
        sub_cmd,
        queue="defq",
        allocation="pclancy3",
        walltime=walltime,
        cpus_per_task=cores,
    )
    # job.wait()
    time.sleep(1)
    os.chdir(BASE_PATH)

    return job


def job(job_config, struc, box_size, atom_dix):
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
    job_str = job_str.replace("$FF_ELEM$", " ".join(job_config["elems"]))
    abc = box_size.split()
    neb_config = job_config["neb"]
    img_path = data_path / job_config["image_base"]
    for atom in struc:
        atom.label = atom_dix[atom.element]

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
            job_config["elems"],
            job_config["overlap"],
        )

    inp_str = job_str.replace("$DATA$", str(data_file))
    inp_str = inp_str.replace("$LATTICE$", str(job_config["lattice"]))
    inp_str = inp_str.replace("$FILE$", str((ff_path / job_config["ff"])))
    inp_str = inp_str.replace("$A$", abc[0])
    inp_str = inp_str.replace("$B$", abc[1])
    inp_str = inp_str.replace("$C$", abc[2])
    inp_str = inp_str.replace("$ID$", " ".join(job_config["id"]))
    inp_str = inp_str.replace("$STEP$", str(neb_config["dx"]))
    inp_str = inp_str.replace("$K1$", str(neb_config["spring"]))
    inp_str = inp_str.replace("$DIST$", neb_config["neigh"])
    inp_str = inp_str.replace("$NSTEP$", str(neb_config["nsteps"]))
    inp_str = inp_str.replace("$NCI$", str(neb_config["ci_steps"]))
    inp_str = inp_str.replace("$IMAGE$", f"{img_path}")

    print("Submitting...")
    jobj = neb_slurm(
        job_name,
        neb_config,
        inp_str,
        prebash=PREBASH,
    )

    return jobj


@click.command()
@click.option("--task", default="ti_100_ex")
@click.option("--name", default=None)
@click.option("--para", default="base")
def run(task, name, para):
    """
    Execute a LAMMPS simulation job based on the provided task and name.

    Args:
        task (str, optional): The task or YAML configuration file name. Defaults to "md_apart_rel".
        name (str, optional): The job name override. Defaults to None.
    """
    atom_dix = {"Ni": 1, "Ti": 2, "Cr": 3}

    systems_path = config_path / "systems" / "neb"
    inputs_path = config_path / "inputs"
    para_path = config_path / "neb"

    task = systems_path / (task + ".yaml")
    setting = para_path / (para + ".yaml")
    inp = inputs_path / "eam_neb.yaml"

    all_sys = list(systems_path.iterdir())
    all_sets = list(para_path.iterdir())
    if task in all_sys:
        yaml_file = task
        with open(yaml_file, "r") as y:
            config = yaml.safe_load(y)
            settings = yaml.safe_load(str(para_path))
            if name is not None:
                config["job_name"] = name
            config["job_str"] = yaml.safe_load(open(inp, "r"))["job_str"]
            struc = files.read_xyz(str(xyz_path / config["xyz"]))[-1]
            if setting in all_sets:
                if para == "base":
                    config["neb"] = yaml.safe_load(open(setting, "r"))
                elif para == "sweep":
                    neb_config = yaml.safe_load(open(setting, "r"))
                    name = config["job_name"]
                    for i, c in enumerate(
                        generate_parameter_combinations(
                            neb_config, ["nframes", "walltime"]
                        )
                    ):
                        config["neb"] = c
                        config["job_name"] = name + f"_sweep{i}"
                        # print(config["job_name"])
                        job(config, struc, config["system_bounds"], atom_dix)
                        # raise Exception
                    return
            else:
                return
            job(
                config,
                struc,
                config["system_bounds"],
                atom_dix,
            )
    else:
        return


if __name__ == "__main__":
    run()
