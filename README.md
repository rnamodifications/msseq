# MS-Seq. (Anchor-based and Bi-directional RNA Sequencing Algorithms)

## Description
A set of LC-MS-based RNA sequencing algorithms.

## System Requirements

### Software requirements
The project has been tested on macOS(Sonoma 14.0).

#### Python dependencies
1. Python version 3.7+ is required. Go to the [Python Website]( https://www.python.org/downloads/) to download python and follow the instruction for installation.
2. The rest dependencies are listed inside requirements.txt, please refer to the [Prepare the Environment](#prepare-the-environment) section on how to install them. They should install within 1 minute, depending on your network speed.

## Download
Clone this repository. Click the green button "Code" and click the "Download ZIP" button. Now a zip file named "msseq-main.zip" is downloaded. Unzip this file. In macOS/Linux system, usually the path of the project is ~/Downloads/msseq-main. 

## Prepare the Environment 
1. Open the Terminal.
2. Create an virtual environment by entering the following command
```Bash
python3 -m venv <your_virtual_workspace_path>
```
To make it simple, we can setup the virtual environment "vir_env" in "Downloads" folder by entering the command
```Bash
python3 -m venv ~/Downloads/vir_env
```
3. Activate the virtual environment
```Bash
source <your_virtual_workspace_path>/bin/activate
```
In our case, enter the command
```Bash
source ~/Downloads/vir_env/bin/activate
```
4. Go to the root directory of msseq project and install all the required libraries. Enter the command
```Bash
cd <root_directory_path>
pip install -r requirements.txt
```
In our case, enter
```Bash
cd ~/Downloads/msseq-main
pip install -r requirements.txt
```

## Sequencing
Run the following command to do Anchor-based sequencing.
```Bash
python seq/process.py -d <LC-MS_dataset> -o <3|5> -l <anchor_mass>
```
Run the following command to do Bi-directional sequencing.
```Bash
python seq/process.py -d <LC-MS_dataset> -o <3|5>
```

## Help
Run the following command to get instructions.
```Bash
python seq/process.py --help
```

