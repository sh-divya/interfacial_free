from squid.utils.units import convert_energy
from squid import files
from supplements.time_energy_gen import read_log_or_out
import numpy as np

import matplotlib.pyplot as plt
import os.path as osp
import numpy as np
import os

LAMMPS_DIR = "/scratch16/pclancy3/divya/interfacial_free/lammps"
STRUCT_DIR = "/scratch16/pclancy3/divya/interfacial_free/structs"

# 4.336, 0.04336, 0.004336, 4.336e-9
# 1e-5 x 3, 1e-7


# legends = [
#     "k=4.336;v=1e-5",
#     "k=0.04336;v=1e-5",
#     "k=0.004336;v=1e-5",
#     "k=4.336e-9;v=1e-9"
# ]

plt.figure()


def read_one_dir(job_name, thermo_format):
    job_dir = osp.join(LAMMPS_DIR, job_name)
    all_out_fptr = [f for f in os.listdir(job_dir) if f.split(".")[-1][0] == "l"]
    print(all_out_fptr)

    for f, fptr in enumerate(all_out_fptr):
        with open(osp.join(job_dir, fptr)) as fobj:
            # dcom = [float(row[-2]) for row in read_log_or_out(fobj, thermo_format)]
            step, energy = zip(
                *[
                    (float(row[0]), convert_energy("eV", "kcal/mol", float(row[1])))
                    for row in read_log_or_out(fobj, thermo_format)
                ]
            )
            # energy = [convert_energy("eV", "kcal/mol", x) for x in energy]
        plt.plot(list(energy)[-500:], label=job_name)

    return list(energy)


def read_ti(job_base, thermo_format, nframes):
    pmf = []
    errors = []
    for i in range(nframes):
        dir_name = job_base + f"_{i}"
        log_file = osp.join(LAMMPS_DIR, osp.join(dir_name, dir_name + ".log"))
        last = 0
        with open(log_file) as fobj:
            f1x, f2x = zip(
                *[
                    (float(row[2]), float(row[5]))
                    for row in read_log_or_out(fobj, thermo_format)
                ]
            )
            m_f1x, m_f2x = np.mean(f1x[-10000:]), np.mean(f2x[-10000:])
            s_f1x, s_f2x = np.std(f1x[-10000:]), np.std(f2x[-10000:])
            pmf.append(convert_energy("eV", "kcal", -1 * (m_f1x + m_f2x)))
            errors.append(convert_energy("eV", "kcal", -1 * (m_f1x + m_f2x)))

    return pmf, errors


if __name__ == "__main__":

    thermo_format = ["Step", "TotEng", "PotEng", "Temp"]
    job_name = "apart_rel_metal"
    init_energy = read_one_dir(job_name, thermo_format)[-500:]
    init_mean = np.mean(init_energy)
    init_std = np.std(init_energy)
    print(init_std)

    job_name = "interface_rel_metal"
    final_energy = read_one_dir(job_name, thermo_format)[-500:]
    final_mean = np.mean(final_energy)
    final_std = np.std(final_energy)
    print(final_std)

    print("Enthalpy of interface formation (kcal/mol)", final_mean - init_mean)
    print("Error", init_std + final_std)

    # plt.figure()
    # job_base = "apart_ti"
    # thermo_format = [
    #     "PotEng",
    #     "Temp",
    #     "f_1[1]",
    #     "f_1[2]",
    #     "f_1[3]",
    #     "f_2[1]",
    #     "f_2[2]",
    #     "f_2[3]",
    # ]

    # ti_frames = []

    # ti_pmf, errs = read_ti(job_base, thermo_format, 23)

    # job_name = "apart_smd"
    # thermo_format = ["PotEng", "Temp", "f_push[6]", "f_push[7]"]
    # dcom = read_one_dir(job_name, thermo_format)
    # plt.plot(dcom[::10000], lammps_pmf[::10000], label="LAMMPS estimate")
    # dcom = dcom[54000:67200:600] + [18.4]
    # ti_pmf = [np.trapz(ti_pmf[:i], dcom[:i]) for i in range(len(ti_pmf))]
    # err_bars = [np.trapz(errs[:i], dcom[:i]) for i in range(len(errs))]
    # # print(len(ti_pmf))
    # plt.scatter(dcom, ti_pmf, color="r", label="TI calculation")
    # plt.errorbar(dcom, ti_pmf, yerr=err_bars, color="r", capsize=5.0)
    # plt.title("SMD with varying pulling parameters")
    # plt.ylabel("Dist between COMs")
    # plt.xlabel("Timesteps (0.001ps)")
    # plt.ylabel("PMF (kcal/mol)")
    plt.legend()
    plt.savefig("nvt.png")
    # plt.savefig("initial_pmf_estimate.png")

    # 4.336, 0.04336, 0.004336, 4.336e-9
    # 1e-5 x 3, 1e-7

