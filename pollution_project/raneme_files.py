import os
import re

# Function to rename files in a directory
def rename_files(directory):
    # List all files in the directory
    files = os.listdir(directory)

    # Define the pattern to match
    pattern = re.compile(r'(out_)(\d{1,5})(\.png)')

    # Iterate through each file
    for filename in files:
        # Match the pattern
        match = pattern.match(filename)
        if match:
            # Extract the parts of the filename
            prefix = match.group(1)
            number = match.group(2)
            extension = match.group(3)

            # Pad the number with zeros if necessary
            new_number = number.zfill(5)

            # Construct the new filename
            new_filename = prefix + new_number + extension

            # Rename the file
            os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))

# Specify the directory containing the files
directory = './results'

# Call the function to rename files in the directory
rename_files(directory)
