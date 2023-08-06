from supplements.time_energy_gen import read_log_or_out
from squid.utils.units import convert_energy
from squid import files

import matplotlib.pyplot as plt
from pathlib import Path
import os.path as osp
import numpy as np
import click
import yaml
import os
import re

BASE_PATH = Path(__file__).parent
LAMMPS_DIR = BASE_PATH / "lammps"
STRUCT_DIR = BASE_PATH / "structs"
CONFIG_PATH = BASE_PATH / "config"


def read_one_dir(
    job_name, thermo_format, indices=[0, 1], unit_conversion=None, plot_log=True
):
    job_dir = osp.join(LAMMPS_DIR, job_name)
    if plot_log:
        all_out_fptr = [f for f in os.listdir(job_dir) if f.endswith(".log")]
    else:
        all_out_fptr = [f for f in os.listdir(job_dir) if re.match(r".+\.o\d+", f)]

    data_list = {idx: [] for idx in indices}

    for f, fptr in enumerate(all_out_fptr):
        with open(osp.join(job_dir, fptr)) as fobj:
            column_indices = [thermo_format.index(idx) for idx in indices]
            step_and_vals = [
                [float(row[i]) for i in column_indices]
                for row in read_log_or_out(fobj, thermo_format)
            ]

        step, yvals = zip(*step_and_vals)

        if unit_conversion is not None:
            yvals = [[unit_conversion(val) for val in yval] for yval in yvals]

        for idx, yval in zip(indices, zip(*yvals)):
            data_list[idx].append(yval)

    return data_list


def read_ti(job_base, thermo_format, nframes):
    pmf = []
    errors = []

    for i in range(nframes):
        dir_name = job_base + f"_{i}"
        data_list = read_one_dir(dir_name, thermo_format, indices=[2, 5])
        f1x, f2x = data_list[2], data_list[5]
        m_f1x, m_f2x = np.mean(f1x[-10000:]), np.mean(f2x[-10000:])
        s_f1x, s_f2x = np.std(f1x[-10000:]), np.std(f2x[-10000:])
        pmf_val = -1 * (m_f1x + m_f2x)
        error_val = np.sqrt(s_f1x**2 + s_f2x**2)

        # Conversion to kcal
        pmf_val = convert_energy("eV", "kcal", pmf_val)
        error_val = convert_energy("eV", "kcal", error_val)

        pmf.append(pmf_val)
        errors.append(error_val)

    return pmf, errors


def process(task):
    systems_path = CONFIG_PATH / "systems"
    simulation = systems_path / (task + ".yaml")
    all_sys = list(systems_path.iterdir())
    if simulation in all_sys:
        yaml_file = simulation
        with open(yaml_file, "r") as y:
            config = yaml.safe_load(y)
            sim = config["sim"]
            thermo_format = config["output"]

            if sim == "md":
                data_list = read_one_dir(task, thermo_format, indices=[0, 1])
                results = data_list
                errors = []
            elif sim == "ti":
                results, errors = read_ti(task, thermo_format, nframes=23)
            elif sim == "smd":
                data_list = read_one_dir(task, thermo_format, indices=[2, 3])
                results = data_list
                errors = []
            else:
                raise ValueError(f"Unsupported simulation type: {sim}")

        return results, errors


if __name__ == "__main__":
    task = "apart_rel_metal"
    results, errors = process(task)
    print("Mean", np.mean(results[-10000:]))
    print("Error", np.std(results[-10000:]))

    # 4.336, 0.04336, 0.004336, 4.336e-9
    # 1e-5 x 3, 1e-7

    # legends = [
    #     "k=4.336;v=1e-5",
    #     "k=0.04336;v=1e-5",
    #     "k=0.004336;v=1e-5",
    #     "k=4.336e-9;v=1e-9"
    # ]
