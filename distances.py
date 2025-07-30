import sys
import subprocess
import os
import pandas as pd
import re
import matplotlib.pyplot as plt
import shutil
import glob

def get_time_step(xtc):
    """
    Summary:
        Get the time step per farme in a trajectory and calculates the time setp of the first 100 frames.

    Parameters:
        xtc: GROMACS trajectory file.
    
    Returns:
        time_step: Time step value that represents the first 100 frames of the trajectory.
    """

    # Run a GROMCAS command-line to check the number of frames of the input trajectory
    result = subprocess.run(f"gmx check -f {xtc}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Get the stdout and stderr from the previous command-line
    output = result.stdout.decode("utf-8") + result.stderr.decode("utf-8")

    # Find the timestep (ps) of 1 single frame from the stdout and stderr
    time_match = re.search(r"^Step\s+\d+\s+(\d+)", output, re.MULTILINE)

    # If it finds, it will calculate the time step needed to extract the desired first number of frames, which is 100, so we multiply the time step per 100 and round it to the units
    if time_match:
        ps = int(time_match.group(1))
        time_step = (ps * 100) + 1
        return time_step
    else:
        # If it doesnt find anything, it will print a debugging message and leave the code
        print("Values not found!!!")
        sys.exit(1)


def extract_first_100_frames(xtc, tpr, time_step):
    """
    Summary: 
        Extract a structure file for each of the first 100 frames of a trajectory.

    Parameters: 
        xtc: GROMACS trajectory file.
        tpr: GROMACS trajectory topology file.
        skip_value: Skip value thats gonna be used in order to determine the number of frames to be extracted.
    """

    # Get the basename of the input topology file, without the type identifier 
    basename = os.path.basename(tpr).split(".")[0]

    # Run a GROMACS + shell script command-line in order to the first 100 frames of the trajectory, converting into .gro files each and making sure the selected option is 'Protein' (option 1)
    subprocess.run(f"echo '1' | gmx trjconv -s {tpr} -f {xtc} -o {basename}_100frames_.gro -sep -b 0 -e {time_step}", shell=True)
    
def calculate_distances(tpr, gro, resid_one, resid_two):
    """
    Summary:
        Creates a dataframe containing every distance value between two residues of every frame extracted.

    Paramaters:
        tpr: GROMACS trajectory topology file.
        gro: GROMACS structure file.
        resid_one: Number of the first residue.
        resid_two: Number of the second residue.

    Returns:
        df: Pandas dataframe containing every distance between the two residues of every frames ectracted before.
        ylims : Limits of the distance values.
    """

    # Run a GROMACS + shell script command-line that creates a customize index that include the group with the two residues e want to analyze
    subprocess.run(f"echo 'r {resid_one}\nr {resid_two}\nq' | gmx make_ndx -f {gro} -o index_{resid_one}_{resid_two}.ndx", shell=True)

    # Create a variable to store the frames number and their respective distance
    data =[]

    # Get the basename of the input topology file, without the type identifier 
    basename = os.path.basename(tpr).split(".")[0]
    
    # Run a GROMACS command-line on every .gro frame file created previously and get the average distance between the two residues we want to analyze
    for i in range(101):
        results = subprocess.run(f"gmx distance -f {basename}_100frames_{i}.gro -s {tpr} -n index_{resid_one}_{resid_two}.ndx -select 'com of group 'r_{resid_one}' plus com of group 'r_{resid_two}''", shell=True, capture_output=True, text=True)
        
        match = re.search(r"Average distance:\s+([\d\.]+)", results.stdout)
        
        # Add the frame number and respective distance to the variable and make sure if the distance was not found it will add a Nonetype 
        if match:
            distance = float(match.group(1))
        else:
            distance = None
        data.append([i, distance])

    # Convert the variable into a Pandas dataframe
    df = pd.DataFrame(data, columns=["Frame", f"Distance between residues {resid_one} and {resid_two} (Å)"])
    
    return df

def get_plots(frame_plot_path, time_plot_path, df, color):
    """
    Summary:
        Generate line plots of the distances by frames and by time (ps).
    Parameters:
        frame_plot_path: Absolute path of the plot containing the distances by frames.
        time_plot_path: Absolute path of the plot containing the distances by time (ps).
        df: Pandas dataframe cotaining the data containing the distances of the first 100 frames.
        color: Plot line color.
    """

    # Create a plot with every registered distance from the two residues during the first 100 frames
    # Define the size of the plot
    plt.figure(figsize=(10, 6))
    # Define which data is gonna be represented by the X value and the Y value, such as the line color, marker size and linewidth
    plt.plot(df["Frame"], df[f"Distance between residues {first_res} and {second_res} (Å)"], marker="o", color=f"{color}", markersize=5, linewidth=1)
    # Define a title for the plot
    plt.title(f"Distances of {first_res} and {second_res} on first 100 frames")
    # Define a label for the X values
    plt.xlabel("Frames")
    # Define a label for the Y values
    plt.ylabel(f"Distance between residues {first_res} and {second_res} (Å)")
    # Rotate the X values tick labels by 45 degrees for better readability
    plt.xticks(rotation=45)
    # Adjust subplot parameters to give specified padding and prevent label overlap
    plt.tight_layout()
    # Define custom plot settings for font size, size of title and size of labels
    plt.rcParams.update({
        "font.size": 12,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
    })
    # Save the plot as a PNG file
    plt.savefig(frame_plot_path, dpi=300)

    # Create a plot with every registered distance from the two residues during the first 10000 ps
    # Define the size of the plot
    plt.figure(figsize=(10, 6))
    # Define which data is gonna be represented by the X value and the Y value, such as the line color, marker size and linewidth
    plt.plot(df["Frame"]*100, df[f"Distance between residues {first_res} and {second_res} (Å)"], marker="o", color=f"{color}", markersize=5, linewidth=1)
    # Define a title for the plot
    plt.title(f"Distances of {first_res} and {second_res} on first 100 frames")
    # Define a label for the X values
    plt.xlabel("Time (ps)")
    # Define a label for the Y values
    plt.ylabel(f"Distance between residues {first_res} and {second_res} (Å)")
    # Rotate the X values tick labels by 45 degrees for better readability
    plt.xticks(rotation=45)
    # Adjust subplot parameters to give specified padding and prevent label overlap
    plt.tight_layout()
    # Define custom plot settings for font size, size of title and size of labels
    plt.rcParams.update({
        "font.size": 12,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
    })
    # Save the plot as a PNG file
    plt.savefig(time_plot_path, dpi=300)



if __name__=="__main__":

    # Make sure every parameter is written
    if len(sys.argv)!=7:
        print("Usage: python3 distances.py {dynamic.xtc} {dynamic.tpr} {dynamic.gro} {first residue} {second residue} {plot line color}")
    else:

        # Define the input parameters to be used
        xtc = sys.argv[1]
        tpr = sys.argv[2]
        gro = sys.argv[3]
        first_res = sys.argv[4]
        second_res = sys.argv[5] 
        color = sys.argv[6]

        # Get the current working directory and basename of the input topology file, without the type identifier
        cwd = os.getcwd()
        basename = os.path.basename(tpr).split(".")[0]

        # Create a counter variable
        i = 1

        # If the storing directory is not created yet, it will create the directory
        if not os.path.exists(f"{basename}_{first_res}_{second_res}_distance"):

            # Create the storing directory
            os.mkdir(f"{basename}_{first_res}_{second_res}_distance")
            # Get the time step value for the frames extraction
            tsetp_value = get_time_step(xtc)
            # Extract the first 100 frames of the trajectory
            extract_first_100_frames(xtc, tpr, tsetp_value)
            # Get a Pandas dataframe with the distances between the two residues over the first 100 frames
            df = calculate_distances(tpr, gro, first_res, second_res)
            # Print the dataframe for debugging
            print(df)

            # Move every frame file to the new directory + the custom index file
            for file in glob.glob(f"{basename}_100frames_*.gro"):
               shutil.move(file, f"{basename}_{first_res}_{second_res}_distance/{file}")
            os.rename(f"{cwd}/index_{first_res}_{second_res}.ndx", f"{cwd}/{basename}_{first_res}_{second_res}_distance/index_{first_res}_{second_res}.ndx")

            # Create the paths of the plot files
            frame_plot_path = f"{cwd}/{basename}_{first_res}_{second_res}_distance/{basename}_{first_res}_{second_res}_distance_by_frames.png"
            time_plot_path = f"{cwd}/{basename}_{first_res}_{second_res}_distance/{basename}_{first_res}_{second_res}_distance_by_time.png"

            # Build the plots with the values f the dataframe
            get_plots(frame_plot_path, time_plot_path, df, color)

        else:
            
            # If the directory was already created, it will create another one with a numeration aside
            while os.path.exists(f"{basename}_{first_res}_{second_res}_distance({i})"):
                i += 1

            # Create the storing directory
            os.mkdir(f"{basename}_{first_res}_{second_res}_distance({i})")
            # Get the time step value for the frames extraction
            tstep_value = get_time_step(xtc)
            # Extract the first 100 frames of the trajectory
            extract_first_100_frames(xtc, tpr, tstep_value)
            # Get a Pandas dataframe with the distances between the two residues over the first 100 frames
            df = calculate_distances(tpr, gro, first_res, second_res)
            # Print the dataframe for debugging
            print(df)

            # Create the path of the plot file
            frame_plot_path = f"{cwd}/{basename}_{first_res}_{second_res}_distance/{basename}_{first_res}_{second_res}_distance_by_frames.png({i})"
            time_plot_path = f"{cwd}/{basename}_{first_res}_{second_res}_distance/{basename}_{first_res}_{second_res}_distance_by_time.png({i})"

            # Move every frame file to the new directory + the custom index file
            for file in glob.glob(f"{basename}_100frames_*.gro"):
               shutil.move(file, f"{basename}_{first_res}_{second_res}_distance({i})/{file}")
            os.rename(f"{cwd}/index_{first_res}_{second_res}.ndx", f"{cwd}/{basename}_{first_res}_{second_res}_distance({i})/index_{first_res}_{second_res}.ndx")

            # Build the plots with the values f the dataframe
            get_plots(frame_plot_path, time_plot_path, df, color)
