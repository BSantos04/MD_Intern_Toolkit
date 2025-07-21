import sys
import os
import matplotlib.pyplot as plt
import pandas as pd

def csv_to_df(csv):
    """
    Summary:
        This function converts a .csv file containing the distances between the frames 403 and 446 of the all-atom structure during the first 10 ns (10000 ps)
    Parameters:
        csv: Input .csv file
    Returns:
        df: Pandas dataframe with every distance between the two residues during time (ps)
    """

    # Create a variable to store the initial data
    data = []

    # Open and parse the file, and then store the data into the variable, creating a Python dataframe
    with open(csv) as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            if line.endswith(","):
                line = line[:-1]
            data.append(line.split(","))

    # Convert the dataframe into a pandas dataframe
    df = pd.DataFrame(data, columns=["Time (ps)", "Distance between residues 446 and 403 (Å)"])
    # Multiply the time per 10000, so it converts from ns to ps, and convert the values from strings to floats (not in that respective order)
    df["Time (ps)"] = df["Time (ps)"].astype(float) * 1000
    df["Distance between residues 446 and 403 (Å)"] = df["Distance between residues 446 and 403 (Å)"].astype(float)
    # Print the dataframe for debugging
    print(df)
    return df
            
if __name__ == "__main__":

    # Make sure that every input parameter is selected
    if len(sys.argv)!=2:
        print("Usage: python3 csv_to_plot.py {.csv file}")
    else:

        # Define input parameter
        csv = sys.argv[1]

        # Get the basename of the .csv input file wihtout the file type identifier
        name = os.path.basename(csv).split(".")[0]

        # Convert the .csv file into a pandas dataframe
        df = csv_to_df(csv)

        # Convert the dataframe into a .PNG plot
        plt.figure(figsize=(10, 6))
        plt.plot(df["Time (ps)"], df["Distance between residues 446 and 403 (Å)"], marker="o", color="red")
        plt.title("Distances of 446 and 403 on first 100 frames")
        plt.xlabel("Time (ps)")
        plt.ylabel("Distance between residues 403 and 446 (Å)")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f"{name}.png")
