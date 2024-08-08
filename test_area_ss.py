from pathlib import Path
from copy import deepcopy
import random
import yaml

from squid import files
from sys_prep import box_info, unwrapz
from meam_md import run_base
from post import read_one_dir

BASE_PATH = Path(__file__).parent
struct_path = BASE_PATH / "structs"
systems_path = BASE_PATH / "config/systems/md"
lammps_path = BASE_PATH / "lammps"

fptr = "agni111_interface_v2_unwrapz.xyz"
log = "agni111_interface_v2"
thermo = ["step", "etotal", "press", "temp", "ke", "pe", "lx", "ly", "lz"]
steps, data = read_one_dir(job_name=log, thermo_format=thermo, indices=[0, 6, 7, 8])
vdw_rad = {"Ni": 1.63, "Ag": 1.72}
frames = files.read_xyz(str(struct_path / fptr))
idx = list(range(1, len(frames)))[-100:]
idx = random.sample(idx, 5)
# idx = [439, 486, 482, 466, 447]
idx = [482, 466, 447]
frames = [frames[i] for i in idx]
frame_start_step = 0
time_start_step = 1250
idx = [(i - 1) * 20000 + time_start_step for i in idx]
atom_pops = [3260, 5184]

for f, frm in enumerate(frames):
    i = steps.index(idx[f])
    a, b, c = [d[i] for d in data]
    # print(a, b, c)

    ag_atoms = frm[: atom_pops[0]]
    ni_atoms = frm[atom_pops[0] :]
    boxes = {}
    boxes["Ag"] = list(box_info(ag_atoms))
    boxes["Ni"] = list(box_info(ni_atoms))
    ag_str = "-".join(f"{boxes['Ag'][-1]:.2f}".split("."))
    ni_str = "-".join(f"{boxes['Ni'][-1]:.2f}".split("."))

    name_agni = f"agni_{ag_str}ag_{ni_str}ni.xyz"
    name_ag = f"stress_{ag_str}ag.xyz"
    name_ni = f"stress_{ni_str}ni.xyz"

    files.write_xyz(frm, str(struct_path / name_agni))
    files.write_xyz(ag_atoms, str(struct_path / name_ag))
    files.write_xyz(ni_atoms, str(struct_path / name_ni))

    for inp in ["minus", "plus"]:
        agni_cfg = yaml.safe_load(
            open(systems_path / "eam_agni_strain_template.yaml", "r")
        )
        inp_yml = f"agni_deform_{inp}"
        name = name_agni.split(".")[0] + f"_{inp}"
        agni_cfg["xyz"] = name_agni
        agni_cfg["job_name"] = name
        agni_cfg["system_bounds"] = [[0, 0, 0], [a, b, c]]
        agni_cfg["sim"] = inp_yml
        fobj = open(systems_path / f"{name}.yaml", "w")
        yaml.dump(agni_cfg, fobj)
        fobj.close()
        run_base(task=name)

    for i, j in zip(["Ag", "Ni"], [name_ag, name_ni]):
        inp_yml = "agni_single_stress"
        name = j.split(".")[0] + "_surf"
        cfg = yaml.safe_load(open(systems_path / "eam_agni_stress_template.yaml", "r"))
        cfg["job_name"] = name
        cfg["xyz"] = j
        cfg["system_bounds"] = [[0, 0, 0], [a, b, boxes[i][-1]]]
        cfg["sim"] = inp_yml
        cfg["elems"] = i
        fobj = open(systems_path / f"{name}.yaml", "w")
        yaml.dump(cfg, fobj)
        fobj.close()
        run_base(task=name)
