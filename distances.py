import sys
import subprocess
import os
import pandas as pd
import re
import matplotlib.pyplot as plt
import shutil
import glob

def extract_first_100_frames(xtc, tpr):
    """

    """
    basename = os.path.basename(tpr).split(".")[0]

    subprocess.run(f"echo '1' | gmx trjconv -s {tpr} -f {xtc} -o {basename}_100frames_.gro -skip 80 -sep", shell=True)
    
def calculate_distances(tpr, gro, resid_one, resid_two):
    
    subprocess.run(f"echo 'r {resid_one}\nr {resid_two}\nq' | gmx make_ndx -f {gro} -o index_{resid_one}_{resid_two}.ndx", shell=True)

    data =[]

    basename = os.path.basename(tpr).split(".")[0]
    
    for i in range(101):
        results = subprocess.run(f"gmx distance -f {basename}_100frames_{i}.gro -s {tpr} -n index_{resid_one}_{resid_two}.ndx -select 'com of group 'r_{resid_one}' plus com of group 'r_{resid_two}''", shell=True, capture_output=True, text=True)
        
        match = re.search(r"Average distance:\s+([\d\.]+)", results.stdout)
        
        if match:
            distance = float(match.group(1))
        else:
            distance = None
        data.append([i, distance])
    
    return pd.DataFrame(data, columns=["Frame", f"Distance between residues {resid_one} and {resid_two} (Å)"])

if __name__=="__main__":
    if len(sys.argv)!=7:
        print("Usage: python3 distances.py {dynamic.xtc} {dynamic.tpr} {dynamic.gro} {first residue} {second residue} {line color}")
    else:
        xtc = sys.argv[1]
        tpr = sys.argv[2]
        gro = sys.argv[3]
        first_res = sys.argv[4]
        second_res = sys.argv[5] 
        color = sys.argv[6]

        cwd = os.getcwd()
        i = 1

        basename = os.path.basename(tpr).split(".")[0]

        if not os.path.exists(f"{basename}_{first_res}_{second_res}_distance"):

            os.mkdir(f"{basename}_{first_res}_{second_res}_distance")
            extract_first_100_frames(xtc, tpr)
            df = calculate_distances(tpr, gro, first_res, second_res)
            print(df)

            for file in glob.glob(f"{basename}_100frames_*.gro"):
               shutil.move(file, f"{basename}_{first_res}_{second_res}_distance/{file}")
            os.rename(f"{cwd}/index_{first_res}_{second_res}.ndx", f"{cwd}/{basename}_{first_res}_{second_res}_distance/index_{first_res}_{second_res}.ndx")

            plt.figure(figsize=(10, 6))
            plt.plot(df["Frame"], df[f"Distance between residues {first_res} and {second_res} (Å)"], marker="o", color=f"{color}")
            plt.title(f"Distances of {first_res} and {second_res} on first 100 frames")
            plt.xlabel("Frames")
            plt.ylabel(f"Distance between residues {first_res} and {second_res}")
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(f"{cwd}/{basename}_{first_res}_{second_res}_distance/{basename}_{first_res}_{second_res}_distance_by_frames.png")

            plt.figure(figsize=(10, 6))
            plt.plot(df["Frame"]*100, df[f"Distance between residues {first_res} and {second_res} (Å)"], marker="o", color=f"{color}")
            plt.title(f"Distances of {first_res} and {second_res} on first 100 frames")
            plt.xlabel("Time (ps)")
            plt.ylabel(f"Distance between residues {first_res} and {second_res}")
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(f"{cwd}/{basename}_{first_res}_{second_res}_distance/{basename}_{first_res}_{second_res}_distance_by_time.png")

        else:
            while os.path.exists(f"{basename}_{first_res}_{second_res}_distance({i})"):
                i += 1

            os.mkdir(f"{basename}_{first_res}_{second_res}_distance({i})")
            extract_first_100_frames(xtc, tpr)
            df = calculate_distances(tpr, gro, first_res, second_res)
            print(df)

            for file in glob.glob(f"{basename}_100frames_*.gro"):
               shutil.move(file, f"{basename}_{first_res}_{second_res}_distance({i})/{file}")
            os.rename(f"{cwd}/index_{first_res}_{second_res}.ndx", f"{cwd}/{basename}_{first_res}_{second_res}_distance({i})/index_{first_res}_{second_res}.ndx")

            plt.figure(figsize=(10, 6))
            plt.plot(df["Frame"], df[f"Distance between residues {first_res} and {second_res} (Å)"], marker="o", color=f"{color}")
            plt.title(f"Distances of {first_res} and {second_res} on first 100 frames")
            plt.xlabel("Frames")
            plt.ylabel(f"Distance between residues {first_res} and {second_res} (Å)")
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(f"{cwd}/{basename}_{first_res}_{second_res}_distance/{basename}_{first_res}_{second_res}_distance_by_frames.png")

            plt.figure(figsize=(10, 6))
            plt.plot(df["Frame"]*100, df[f"Distance between residues {first_res} and {second_res} (Å)"], marker="o", color=f"{color}")
            plt.title(f"Distances of {first_res} and {second_res} on first 100 frames")
            plt.xlabel("Time (ps)")
            plt.ylabel(f"Distance between residues {first_res} and {second_res} (Å)")
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(f"{cwd}/{basename}_{first_res}_{second_res}_distance/{basename}_{first_res}_{second_res}_distance_by_time.png")


