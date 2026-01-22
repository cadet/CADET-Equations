# CADET-Equations

[![DOI](https://zenodo.org/badge/936276911.svg)](https://doi.org/10.5281/zenodo.18339286)
[![CI](https://github.com/cadet/cadet-equations/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/cadet/cadet-equations/actions/workflows/ci.yml?query=branch%3Amaster)

CADET-Equations is a Python tool designed to generate modeling equations for packed-bed chromatography.  
It provides a simple user interface that allows users to configure the model, and outputs the corresponding mathematical equations in LaTeX format.

## Installation

To get started with CADET-Equations, you need to install the required Python dependencies and a LaTeX compiler.

### Step 1: Install Python Dependencies

Create a new Conda environment with Python 3.10:

```bash
conda create -n cadet-equations python=3.10

Next, install the necessary Python packages using `pip` by running the following commands:

```bash
conda activate cadet-equations

conda install pip

pip install -r requirements.txt
```

### Step 2: Install LaTeX Compiler

CADET-Equations requires a LaTeX compiler to generate mathematical equations in LaTeX format. 
1. If you are using a Linux-based system, you can install the necessary LaTeX packages using the following command:

```bash
# Install LaTeX compiler
sudo apt-get install texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra
```
2. If you use windows, do as follows:
### 1. Download MiKTeX

Download the MiKTeX installer for Windows from the official site:

👉 https://miktex.org/download

Choose the **"Net Installer"** for Windows (64-bit).

---

### 2. Install MiKTeX

- Run the downloaded `.exe` installer.
- Choose **"Install for just me"** (unless you need system-wide installation).
- During installation, make sure to enable:

  ✅ **"Install missing packages on-the-fly"**

> This is important: it allows MiKTeX to automatically fetch any required LaTeX packages when needed by your app or script.

---

### 3. Verify the Installation

After installation is complete:

- Open a **new Command Prompt** (or Miniforge Prompt)
- Run the following command:

```bash
pdflatex --version

This will install the required LaTeX components needed to render equations properly.

## Usage

Once you have installed the necessary dependencies, you can use the tool to generate packed-bed chromatography modeling equations.

```bash
# Run the tool (replace with your script or command)
streamlit run Equation-Generator.py
```

### Notes
- Ensure that the LaTeX compiler is correctly installed and accessible from your terminal.
- If you encounter any issues with LaTeX rendering, check the LaTeX installation or adjust configurations as needed.

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0) - see the [LICENSE](LICENSE) file for details.
