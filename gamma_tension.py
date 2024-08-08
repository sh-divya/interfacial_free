from post import read_one_dir
import numpy as np
from pathlib import Path
from squid.lammps.io import read_dump

BASE_PATH = Path(__file__).parent
structs_path = BASE_PATH / "structs"
lammps_path = BASE_PATH / "lammps"


def deltaU(job_name, Uinit, kb, T):
    eV_to_Joule = 1.60218e-19
    try:
        thermo = ["step", "etotal", "press", "temp", "ke", "pe"]
        _, data = read_one_dir(job_name=job_name, thermo_format=thermo, indices=[0, 5])
    except IndexError:
        thermo = ["step", "etotal", "press", "temp", "ke", "pe", "lx", "ly", "lz"]
        _, data = read_one_dir(job_name=job_name, thermo_format=thermo, indices=[0, 5])
    # delU = [
    #     np.exp(-1 * ((d - Uinit) * eV_to_Joule) / (kb * T)) for d in data[0][-1000:]
    # ]
    delU = [(d - Uinit) * eV_to_Joule for d in data[0][-1000:]]

    return np.mean(delU)


def test_area_method(plus, minus, kb, T, delA):
    # kbT_delA = kb * T / (4 * delA)
    # return -1 * kbT_delA * (np.log(plus) - np.log(minus))
    return (plus - minus) / (delA * 4)


def bulk_stress_correction(job_name, vol):
    frames = read_dump(
        str(structs_path / f"{job_name}_wbox.dump"),
        extras=["c_allStress[1]", "c_allStress[2]", "c_allStress[3]"],
    )[-50:]
    stress = []

    for frm in frames:
        stress_per_frame = [0.0, 0.0, 0.0]
        for atom in frm:
            stress_per_atom = atom.extras
            stress_per_frame[0] += float(stress_per_atom["c_allStress[1]"])
            stress_per_frame[1] += float(stress_per_atom["c_allStress[2]"])
            stress_per_frame[2] += float(stress_per_atom["c_allStress[3]"])
        stress.append(stress_per_frame)

    stress = np.mean(stress, axis=0)
    stress = 0.5 * (stress[0] + stress[1]) - stress[2]

    return 1e5 * stress / vol


def gamma(g_ta, l1, l2, corr1, corr2):
    return g_ta - 0.5 * l1 * corr1 - 0.5 * l2 * corr2


if __name__ == "__main__":
    # test = bulk_stress_correction("stress_37-85ag_surf", 5954970.78)
    # print(test)
    # raise Exception

    Lz_ag = [37.85, 38.31, 38.25, 38.09, 38.01]
    Lz_ni = [35.70, 35.80, 35.83, 35.99, 35.85]
    ab = [[30.14, 52.2], [30.11, 52.14], [30.14, 32.2], [30.13, 52.19], [30.11, 52.17]]

    plus_minus = "agni_$lag$ag_$lni$ni"
    stress_ag = "stress_$lag$ag_surf"
    stress_ni = "stress_$lni$ni_surf"

    refU = -31213.14  # eV
    boltzmann = 1.38064852e-23  # In Joules
    Temperature = 823
    # eV_to_Joule = 1.60218e-19

    for i, (l1, l2) in enumerate(zip(Lz_ag, Lz_ni)):
        l1_str = "-".join(f"{l1:.2f}".split("."))
        l2_str = "-".join(f"{l2:.2f}".split("."))
        plus_sim = (
            plus_minus.replace("$lag$", l1_str).replace("$lni$", l2_str) + "_plus"
        )
        minus_sim = (
            plus_minus.replace("$lag$", l1_str).replace("$lni$", l2_str) + "_minus"
        )
        ag_sim = stress_ag.replace("$lag$", l1_str)
        ni_sim = stress_ni.replace("$lni$", l2_str)
        plus_deltaU = deltaU(plus_sim, refU, boltzmann, Temperature)  # * eV_to_Joule
        minus_deltaU = deltaU(minus_sim, refU, boltzmann, Temperature)  # * eV_to_Joule
        ogA = ab[i][0] * ab[i][1]
        delA = 0.005 * ogA * 1e-20  # m^2
        gamma_ta = test_area_method(
            plus_deltaU, minus_deltaU, boltzmann, Temperature, delA
        )
        print(plus_sim, minus_sim, ag_sim, ni_sim)
        print("G_ta", gamma_ta)
        vol_ag = l1 * ogA
        bulk_ag = bulk_stress_correction(ag_sim, vol_ag)
        vol_ni = l2 * ogA
        bulk_ni = bulk_stress_correction(ni_sim, vol_ni)
        gamma_ss = gamma(gamma_ta, l1 * 1e-10, l2 * 1e-10, bulk_ag, bulk_ni)
        print("G_ss", gamma_ss)
