from squid import files, structures
from squid.jobs import slurm

import os
from pathlib import Path
import click
import yaml

BASE_PATH = Path(__file__).parent.resolve()
xyz_path = BASE_PATH / "structs"
data_path = BASE_PATH / "data"
ff_path = BASE_PATH / "force_fields/pseudo"
config_path = BASE_PATH / "config"

DELTA = 0.25
DELTA = 0

MASSES = {1: 58.6934, 2: 47.867, 3: 51.9961}
TYPES = ["Ni", "Ti", "Cr"]


def dump_atoms_data(system, boundary="xyz"):
    max_x = max([a.x for a in system.atoms])
    min_x = min([a.x for a in system.atoms])
    max_y = max([a.y for a in system.atoms])
    min_y = min([a.y for a in system.atoms])
    max_z = max([a.z for a in system.atoms])
    min_z = min([a.z for a in system.atoms])
    atoms_idx = []
    elems = set([a.label for a in system.atoms])
    if boundary is not None:
        for i, a in enumerate(system.atoms):
            if a.y < max_y - DELTA:
                if a.z < max_z - DELTA:
                    atoms_idx.append(i)
    else:
        atoms_idx = list(range(len(system.atoms)))
    atoms_data = []
    for i, a in enumerate(system.atoms):
        if i in atoms_idx:
            atoms_data.append(f"{a.label}\t{a.x: 12.10f}\t{a.y: 12.10f}\t{a.z: 12.10f}")
    return (
        "\n".join(atoms_data),
        len(atoms_idx),
        elems,
        [min_x, max_x, min_y, max_y, min_z, max_z],
    )


def write_input_script(system, job_str, abc, seconds=258400):
    atoms_data, num_atoms, num_elems, box_dims = dump_atoms_data(system, None)
    inp_str = job_str.replace("$NAME$", system.name)
    inp_str = inp_str.replace("$TIME$", str(seconds))
    inp_str = inp_str.replace("$NUM$", str(num_atoms))
    inp_str = inp_str.replace("$TYPE$", str(len(num_elems)))
    inp_str = inp_str.replace("$POSITIONS$", atoms_data)
    a = box_dims[1] - box_dims[0]
    b = box_dims[3] - box_dims[2]
    c = box_dims[5] - box_dims[4]
    if abc:
        a, b, c = list(map(float, abc.split(" ")))
    abc = f"{a: <12f}\t0.0000000000\t0.0000000000\n0.0000000000\t{b: <12f}\t0.0000000000\n0.0000000000\t0.0000000000\t{c: <12f}"
    inp_str = inp_str.replace("$BOX$", abc)
    pseudo = ""
    for e in num_elems:
        pp_path = list(ff_path.glob(f"{TYPES[e - 1].lower()}_pbe_v1*.uspp.F.UPF"))[0]
        pseudo += f"{e}\t{MASSES[e]}\t{pp_path.name}\n"
    inp_str = inp_str.replace("$PSEUDO$", pseudo)
    return inp_str


def create_job_dir(job_name, inp_str, clean=False):
    inp_dir = BASE_PATH / "qe" / job_name
    try:
        inp_dir.mkdir(parents=True, exist_ok=not clean)
    except FileExistsError:
        for f in inp_dir.iterdir():
            f.unlink()
    inp_path = inp_dir / (job_name + ".in")
    with open(inp_path, "w") as f:
        f.write(inp_str)

    return inp_dir, inp_path


def job(sim, job_name, struc, atom_dix, abc, kp, start):
    for atom in struc:
        atom.label = atom_dix[atom.element]

    tmp = structures.Molecule(struc)
    box = structures.System(job_name)
    box.add(tmp)
    box.atom_labels = list(atom_dix.values())
    inputs_path = config_path / "inputs"
    slurm_str = yaml.safe_load(open(inputs_path / "qe_sub.yaml", "r"))["slurm_str"]
    qe_inp_temp = inputs_path / (sim + ".yaml")
    job_str = yaml.safe_load(open(qe_inp_temp, "r"))["job_str"]

    job_str = job_str.replace("$KPOINTS$", kp)
    job_str = job_str.replace("$START$", f"'{start}'")
    qe_inp_str = write_input_script(box, job_str, abc)
    job_dir, job_inp = create_job_dir(job_name, qe_inp_str, clean=True)

    slurm_inp = slurm_str.replace("$INPUT$", str(job_inp))
    job_out = job_dir / (job_name + ".out")
    slurm_inp = slurm_inp.replace("$OUTPUT$", str(job_out))

    return job_dir, slurm_inp


@click.command()
@click.option("--task", default="cr_unit_vc")
@click.option("--name", default=None)
@click.option("--start", default="from_scratch")
def sub(task, name, start):
    systems_path = config_path / "systems/qe"
    task = systems_path / (task + ".yaml")

    all_sys = list(systems_path.iterdir())
    if task in all_sys:
        yaml_file = task
        with open(yaml_file, "r") as y:
            config = yaml.safe_load(y)
            if name is not None:
                config["job_name"] = name
            struc = files.read_xyz(str(xyz_path / config["xyz"]))
            atom_dix = config["atoms"]
            abc = config["abc"]
            kp = config["kpoints"]
            job_dir, slurm_inp = job(
                config["sim"], config["job_name"], struc, atom_dix, abc, kp, start
            )
            kp_mult = 1
            for k in kp.split(" "):
                kp_mult = kp_mult * int(k)
            os.chdir(job_dir)
            ncores = config["cores"]
            if ncores >= kp_mult:
                kp_mult = ncores // kp_mult
            else:
                kp_mult = kp_mult // ncores
            nodes = (ncores // 48) + 1
            slurm_inp = slurm_inp.replace("$CORES$", str(ncores))
            slurm_inp = slurm_inp.replace("$NK$", str(kp_mult))
            print(f"Submitting...{config['job_name']}")
            # print(BASE_PATH)
            # print(slurm_inp)
            slurm.submit_job(
                config["job_name"],
                slurm_inp,
                queue="defq",
                walltime=config["walltime"],
                ntasks=ncores,
                allocation="pclancy3",
            )
            os.chdir(BASE_PATH)


if __name__ == "__main__":
    sub()
