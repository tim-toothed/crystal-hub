nja_cfs/
├── __init__.py
├── additional/
│   ├── full_basis.py       # According to pdf
├── magnetic/
│   ├── __init__.py
│   ├── properties.py       # Magnetics class
│   └── numba_functions.py  # JIT-compiled functions
├── crystal_field/
│   ├── __init__.py
│   ├── parameters.py       # Bkq, AOM conversions
│   ├── rotations.py        # LF rotation functions
│   └── pcm.py              # Point charge model
├── io/
│   ├── __init__.py
│   ├── readers.py          # read_xyz, read_AILFT_orca6, etc.
│   └── writers.py          # w_inp_dat, etc.
├── utils/
│   ├── __init__.py
│   ├── coordinates.py      # Cartesian/spherical conversions
│   ├── projections.py      # projection functions
│   └── helpers.py          # utility functions
├── tables/
│   ├── __init__.py
│   ├── data.py             # Reference data and tables
│   └── cfp_*.txt           # Data files
└── visualization/
    ├── __init__.py
    ├── plotting.py         # plot functions
    └── figures.py          # fig_* functions