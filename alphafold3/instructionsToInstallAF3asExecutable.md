Instructions on how to install Alphafold3 as an executable package locally so we can modify it (adapted from https://github.com/pyDock/AlphaFold3-Conda-Install):

### 0. Install Miniconda

    # Download the Miniconda installer
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh

    # Run the installer
    bash miniconda.sh

    # Source the .bashrc to update your environment
    source ~/.bashrc

### 1. Create a Conda Environment with Python 3.11

    # Create the environment
    conda create -n Alphafold3 python=3.11

    # Activate the environment
    conda activate Alphafold3

    # Prevent Python from using ~/.local for user-installed packages
    # This ensures the environment is fully isolated ("closed")
    conda env config vars set PYTHONUSERBASE=intentionally-disabled

    # Reactivate the environment for the change to take effect
    conda deactivate
    conda activate Alphafold3

### 2. Install Development Tools and Dependencies

    # Install Boost for Python 3.11 and numpy
    conda install -c conda-forge boost boost-cpp numpy -y # Boost for Python 3.11

    # Install compression libraries
    conda install -c conda-forge bzip2 zstd -y

    # Install git and zlib
    conda install -c conda-forge zlib -y

    # Install HMMER
    conda config --add channels bioconda 
    conda install -c conda-forge hmmer gsl=2.6 -y

    # Install and Upgrade pip within the Alphafold3 environment
    conda install pip -y
    pip install --upgrade pip  # Update pip (specific to the AF3 environment)

### 3. Install Required Python Packages with pip

    pip install pandas==2.2.3 matplotlib==3.10.0 absl-py==2.1.0 chex==0.1.87 \
        dm-haiku==0.0.13 dm-tree==0.1.8 filelock==3.16.1 \
        "jax[cuda12]==0.4.34" jax-cuda12-pjrt==0.4.34 jax-triton==0.2.0 \
        jaxlib==0.4.34 jaxtyping==0.2.34 jmp==0.0.4 ml-dtypes==0.5.0 \
        numpy==2.1.3 nvidia-cublas-cu12==12.6.3.3 \
        nvidia-cuda-cupti-cu12==12.6.80 nvidia-cuda-nvcc-cu12==12.6.77 \
        nvidia-cuda-runtime-cu12==12.6.77 nvidia-cudnn-cu12==9.5.1.17 \
        nvidia-cufft-cu12==11.3.0.4 nvidia-cusolver-cu12==11.7.1.2 \
        nvidia-cusparse-cu12==12.5.4.2 nvidia-nccl-cu12==2.23.4 \
        nvidia-nvjitlink-cu12==12.6.77 opt-einsum==3.4.0 pillow==11.0.0 \
        rdkit==2024.3.5 scipy==1.14.1 tabulate==0.9.0 toolz==1.0.0 \
        tqdm==4.67.0 triton==3.1.0 typeguard==2.13.3 \
        typing-extensions==4.12.2 zstandard==0.23.0

### 4. Install AlphaFold 3

#### 4.1 Clone AF inference repository

    # Set the desired application directory
    export APPDIR="/home/user/Programs"  # Replace "/home/user/Programs" with your desired path

    # Create the directory and navigate to it
    mkdir -p $APPDIR
    cd $APPDIR

    # Clone the AlphaFold 3 repository
    git clone https://github.com/google-deepmind/alphafold3.git

    # Define the AlphaFold 3 directory variable
    export ALPHAFOLD3DIR="$APPDIR/alphafold3"
    cd ${ALPHAFOLD3DIR}

#### 4.2 Download db (optional, we already downloaded them so just copy)

    # Modify the download path in the script
    sed -i 's|$HOME|$ALPHAFOLD3DIR|g' fetch_databases.sh

    # Make the script executable
    chmod +x fetch_databases.sh

    # Run the script to download the databases
    ./fetch_databases.sh

#### 4.3 Obtain Model Parameters and Place Them in models (require them using AF formulair)

#### 4.4 Install AlphaFold 3 from the Repository

    cd ${ALPHAFOLD3DIR}

    # Export paths for zlib
    export CXXFLAGS="-I$(dirname $(find ${CONDA_PREFIX} -name zlib.h | head -n 1))"
    export LDFLAGS="-L$(dirname $(find ${CONDA_PREFIX} -name libz.so | head -n 1)) -lz"

    # Install AlphaFold 3 without additional dependencies
    python -m pip install --no-deps -e .

#### 4.5 Build additional complements

    cd ${CONDA_PREFIX}/bin

    # Execute the build script
    ./build_data  # Execute

#### 4.6 Test the installation

    cd ${ALPHAFOLD3DIR}

    # Display the help message
    python run_alphafold.py --help

### 5.Fix JAX error (if it appears)

If a segmentation fault occurs when initializing pjrt_plugin, it may be due to import order conflicts between JAX and other C++/CUDA extensions (e.g. alphafold3.cpp, SciPy, RDKit). A simple workaround is to ensure jax and jnp are imported before any other modules.

This command moves the imports to the top (before from collections.abc ...) and removes duplicate occurrences later in the file:

    sed -i".bk" -e 's|from collections\.abc import Callable, Sequence|import jax\nfrom jax import numpy as jnp\n\nfrom collections.abc import Callable, Sequence|' \
    -e '/^import jax$/d' \
    -e '/^from jax import numpy as jnp$/d' ${PWD}/run_alphafold.py