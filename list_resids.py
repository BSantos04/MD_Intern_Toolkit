import sys
import re
import pandas as pd


def list_res(file, res_1, res_2, outfile):
    """
    Summary:
        List the location of two choosen residues in a GRO file and creates a .csv file and a .xlsx file containing the atom index, residue number and respective residue name listed on a table
    Parameters:
        file: Input GRO file
        res_1: Abreviation of the first residue name
        res_2: Abreeviation of the second residue name
        outfile: Base name of the output file
    """

    # Define a counter variable in order to register the index of each stored residue
    counter = 0
    # Define a list varibale to store the residue index, residue number and residue abreviation, creating a Python matrix
    res = []

    # Open the input file
    with open(file) as infile:

        # Run on every line of the file to find the desired residues, making sure to reomve every space at the start and end of every line
        for line in infile:
            line = line.strip()
            # Use regular expressions to find for the desired residues
            match = re.match(r'\s*(\d+)([A-Z]{3})\s+\S+\s+\d+', line)
            # If the line finds a match for the pattern, it will also find for a match for the desired residue
            if match:
                counter += 1
                aa = match.group(2).upper()
                # If it finds a match for the desired residue, it will add the index, number and name of that to the Python matrix previously created
                if aa == res_1.upper() or aa == res_2.upper():
                    res.append([counter, match.group(1), aa])
    
    # Creates a pandas dataframe from the Python matrix created, define the name of each column, and print it for debugging
    df = pd.DataFrame(res, columns=["Atom Index", "Number of the residue", "Name of the aminoacid"])
    print(df)

    # Creates a .csv and a .xlsx files based on the dataframe and print a debugging message
    df.to_csv(f"{outfile}.csv", index=False)
    df.to_excel(f"{outfile}.xlsx", index=False)
    print("Excel and csv files saved sucessfully!!!")


if __name__ == "__main__":
    # Make sure every needed parameter is met, and if not, it will exit the program
    if len(sys.argv)!=5:
        print("Usage: python3 list_resids.py {path of the input file} {name of the first residue (Ex.: SER)} {name of the second residue (Ex.: TYR)} {basename of the output file}")
        sys.exit(1)
    else:
        # Define every input parameter
        infile = sys.argv[1]
        res_1 = sys.argv[2]
        res_2 = sys.argv[3]
        outfile = sys.argv[4]

        # Generate the .csv and .xlsx files with desired aminoacids listed
        list_res(infile, res_1, res_2, outfile)


