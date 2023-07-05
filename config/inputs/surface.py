job_str = """
&CONTROL
  calculation = 'vc-relax'
  etot_conv_thr =   2.7000000000d-04
  forc_conv_thr =   1.0000000000d-04
  outdir = './out/'
  prefix = '$NAME$'
  pseudo_dir = '/scratch16/pclancy3/divya/interfacial_free/force_fields/pseudo'
  tprnfor = .true.
  tstress = .true.
  verbosity = 'high'
/
&SYSTEM
  degauss =   2.2049585400d-02
  ecutrho =   3.2000000000d+02
  ecutwfc =   4.0000000000d+01
  ibrav = 0
  nat = $NUM$
  nosym = .true.
  nspin = 1
  ntyp = $TYPE$
  occupations = 'smearing'
  smearing = 'cold'
/
&ELECTRONS
  conv_thr =   4.0000000000d-10
  electron_maxstep = 200
  mixing_beta =   4.0000000000d-01
  mixing_mode = 'local-TF'
/
&IONS
  ion_positions = 'default'
  ion_dynamics = 'bfgs'
/
&CELL
  cell_dofree='ibrav'
/
ATOMIC_SPECIES
$PSEUDO$
ATOMIC_POSITIONS angstrom
$POSITIONS$
K_POINTS automatic
4 2 2 0 0 0
CELL_PARAMETERS angstrom
$BOX$
"""

