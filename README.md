# CADET-Equations

[![DOI](https://zenodo.org/badge/936276911.svg)](https://doi.org/10.5281/zenodo.18339286)
[![CI](https://github.com/cadet/cadet-equations/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/cadet/cadet-equations/actions/workflows/ci.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/cadet/CADET-Equations/graph/badge.svg?token=CU85FHDKOO)](https://codecov.io/gh/cadet/CADET-Equations)
[![Hosted on Streamlit](https://img.shields.io/badge/Hosted%20on-Streamlit-FF4B4B?logo=streamlit&logoColor=white)](https://cadet-equations-74ko8eryoxmsbqspggxrj2.streamlit.app)

CADET-Equations is a Python tool designed to generate modeling equations for packed-bed chromatography.  
It provides a simple user interface that allows users to configure the model, and outputs the corresponding mathematical equations in LaTeX format.
CADET-Equations is [hosted on Streamlit](https://cadet-equations-74ko8eryoxmsbqspggxrj2.streamlit.app).

## Installation

To get started with CADET-Equations, you need to clone the repository and install the
required Python dependencies and a LaTeX compiler.

### Step 1: Clone the Repository

```
git clone https://github.com/cadet/CADET-Equations.git
```

```
cd CADET-Equations
```

### Step 2: Install Python Dependencies

Create a new Conda environment with Python 3.10:

```
conda create -n cadet-equations python=3.10
```

Activate the environment:

```
conda activate cadet-equations
```

Install `pip`:

```
conda install pip
```

Install the required packages:

```
pip install -r requirements.txt
```

### Step 3: Install LaTeX Compiler

CADET-Equations requires a LaTeX compiler to generate mathematical equations in LaTeX format.

#### Linux

Install the necessary LaTeX packages:

```
sudo apt-get install texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra
```

#### Windows

1. Download the MiKTeX installer from https://miktex.org/download and choose the
   **"Net Installer"** for Windows (64-bit).
2. Run the downloaded `.exe` installer.
3. Choose **"Install for just me"** (unless you need a system-wide installation).
4. Enable **"Install missing packages on-the-fly"**. This lets MiKTeX fetch any required
   LaTeX packages automatically when the app needs them.

Then open a **new Command Prompt** (or Miniforge Prompt) and verify the installation:

```
pdflatex --version
```

## Usage

Once you have installed the necessary dependencies, run the tool from the repository root
to generate packed-bed chromatography modeling equations:

```
streamlit run Equation-Generator.py
```

### Notes
- Ensure that the LaTeX compiler is correctly installed and accessible from your terminal.
- If you encounter any issues with LaTeX rendering, check the LaTeX installation or adjust configurations as needed.

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0) - see the [LICENSE](LICENSE) file for details.
