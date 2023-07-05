from squid import files, structures
from squid.jobs import slurm

from config.qe_sub import slurm_str

# from config.bulk_rel import job_str

import os
from pathlib import Path
import click


base_path = Path("/scratch16/pclancy3/divya/interfacial_free/qe")
xyz_path = Path("/scratch16/pclancy3/divya/interfacial_free/structs")
data_path = Path("/scratch16/pclancy3/divya/interfacial_free/data")
ff_path = Path("/scratch16/pclancy3/divya/interfacial_free/force_fields/pseudo")


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
            if a.x < max_x:
                if a.x > min_x:
                    if a.y < max_y:
                        if a.y > min_y:
                            if a.z < max_z:
                                if a.z > min_z:
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


def write_input_script(system, job_str, abc, seconds=259000):
    atoms_data, num_atoms, num_elems, box_dims = dump_atoms_data(system, None)
    inp_str = job_str.replace("$NAME$", system.name)
    inp_str = job_str.replace("$TIME$", str(seconds))
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
    inp_dir = base_path / job_name
    try:
        inp_dir.mkdir(parents=True, exist_ok=not clean)
    except FileExistsError:
        for f in inp_dir.iterdir():
            f.unlink()
    inp_path = inp_dir / (job_name + ".in")
    with open(inp_path, "w") as f:
        f.write(inp_str)

    return inp_dir, inp_path


def job(job_name, struc, atom_dix, abc):
    for atom in struc:
        atom.label = atom_dix[atom.element]

    tmp = structures.Molecule(struc)
    box = structures.System(job_name)
    box.add(tmp)
    box.atom_labels = list(atom_dix.values())
    if job_name in ["niti110", "cr110", "niti110_stopped", "cr110_stopped"]:
        from config.inputs.surf_rel import job_str
    elif job_name in ["niti", "cr"]:
        from config.inputs.bulk_rel import job_str
    elif job_name in ["crniti_interface", "interface_stopped"]:
        from config.inputs.interface import job_str
    qe_inp_str = write_input_script(box, job_str, abc)
    job_dir, job_inp = create_job_dir(job_name, qe_inp_str, clean=True)

    slurm_inp = slurm_str.replace("$INPUT$", str(job_inp))
    job_out = job_dir / (job_name + ".out")
    slurm_inp = slurm_inp.replace("$OUTPUT$", str(job_out))

    return job_dir, slurm_inp


@click.command()
@click.option("--system", default="0")
@click.option("--ncores", default="64")
@click.option("--kp", default="8")
@click.option("--walltime", default="72:00:00")
@click.option("--abc", default=None)
def sub(system, ncores, kp, walltime, abc):
    xyz_list = [
        "niti110.xyz",
        "cr110.xyz",
        "crniti_interface.xyz",
        "niti.xyz",
        "cr.xyz",
        "niti110_stopped.xyz",
        "cr110_stopped.xyz",
        "interface_stopped.xyz",
    ]
    for s in system:
        xyz = xyz_list[int(s)]
        struc = files.read_xyz(str(xyz_path / xyz))
        job_name = xyz.split(".")[0]
        atom_dix = {
            "0": {"Ni": 1, "Ti": 2},
            "1": {"Cr": 3},
            "2": {"Cr": 3, "Ni": 1, "Ti": 2},
            "3": {"Ni": 1, "Ti": 2},
            "4": {"Cr": 3},
            "5": {"Ni": 1, "Ti": 2},
            "6": {"Cr": 3},
            "7": {"Cr": 3, "Ni": 1, "Ti": 2},
        }
        job_dir, slurm_inp = job(job_name, struc, atom_dix[s], abc)
        os.chdir(job_dir)
        ncores = int(ncores)
        nodes = (ncores // 24) + 1
        slurm_inp = slurm_inp.replace("$CORES$", str(ncores))
        slurm_inp = slurm_inp.replace("$NK$", kp)
        slurm.submit_job(
            job_name,
            slurm_inp,
            queue="defq",
            walltime=walltime,
            ntasks=ncores,
            allocation="pclancy3",
        )
        os.chdir(base_path.parent)


if __name__ == "__main__":
    sub()
    # niti 110 vcrel: "6.35681 + 15 14.8326 10.4882"
    # cr 110 vc-rel: "6.0411 + 15, 14.0959, 9.9673"
    # niti unit vc-rel: "3.01051 3.01051 3.01051"
    # cr unit vc-rel:  "2.96890 2.96890 2.96890"
    # interface vc-rel: "15.3431 + 15, 14.8326, 10.4882"

    # cr -   20.706259434  -0.006745908  -0.017097752
    # -0.004443549  14.115987471  -0.010719601
    # -0.015684556  -0.006835575  11.068669984

    # niti -   20.063879787   0.000029337   0.006241442
    # 0.000019306  14.751828299   0.000003000
    # 0.001221396   0.000002150  10.929892144

    # interface -   30.632191943   0.015812917  -0.012826124
    #               0.007728923  15.292166742  -0.005533737
    #              -0.004415996  -0.003905910  10.820216095
