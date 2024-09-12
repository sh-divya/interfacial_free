from post import read_one_dir
import numpy as np
from pathlib import Path
from squid import files
import yaml
from sys_prep import box_info

BASE_PATH = Path(__file__).parent
sys_path = BASE_PATH / "config" / "systems" / "md"
structs_path = BASE_PATH / "structs"
lammps_path = BASE_PATH / "lammps"



if __name__ == "__main__":
    tasks = [
        # "meam_cr_unit",
        # "dmd_niti_unit",
        # "dmd_niti_bulk",
        # "meam_cr_bulk",
        "dmd_niti_free1ni",
        "dmd_niti111ni_free2",
        "dmd_niti_free2ni",
        "dmd_niti111ti_free2",
        "dmd_niti_free2ti",
        "dmd_niti_free1ti",
        "dmd_niti110_free2",
        "dmd_niti110_free1",
        "dmd_niti111ni_free1",
        "dmd_niti111ti_free1",
        # "meam_cr111_free2",
        # "meam_cr110_free1",
        # "meam_cr111_free1",
        # "meam_cr_free2",
        # "meam_cr110_free2",
        # "meam_cr_free1",
        # "meam_crniti_001ni",
        # "meam_crniti_110",
        # "meam_crniti_111ni",
        # "meam_crniti_111ti",
        # "meam_crniti_001ti"
    ]
    for t in tasks:
        yml_fptr = sys_path / f"{t}.yaml"
        sys = yaml.safe_load(open(str(yml_fptr)))
        name = sys["job_name"]
        thermo = sys["output"]
        xyz = structs_path / f"{name}.xyz"
        frames = files.read_xyz(str(xyz))[-100:]
        areas = []
        for f in frames:
            lx, ly, _ = box_info(f)
            areas.append(lx * ly)
        area = np.mean(areas)
        _, e = read_one_dir(name, thermo)
        energy = np.mean(e[0][-1000:])
        print(t, area, energy)