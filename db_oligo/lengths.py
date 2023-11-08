
# Import necessary libraries
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

# Set the filename to be read
filename = "HK5_He_out.full_wRE_noPrim.genes"

# Initialize lists to store sequence lengths and sequence names
seq_lengths = []
seq_names = []

# Use the SeqIO.parse() function to iterate over the sequences in the file
for record in SeqIO.parse(filename, "fasta"):
    # Add the sequence length to the list
    seq_lengths.append(len(record.seq))
    # Add the sequence name to the list
    seq_names.append(record.id)

# Convert the list of sequence lengths to a numpy array
seq_lengths = np.array(seq_lengths)

# Calculate the median sequence length
median_length = np.median(seq_lengths)

# Calculate the maximum sequence length
max_length = np.max(seq_lengths)

# Calculate the minimum sequence length
min_length = np.min(seq_lengths)

# Create a histogram of sequence lengths
plt.hist(seq_lengths, bins=50, range=(900, 1100))
plt.title("Sequence Lengths")
plt.xlabel("Length")
plt.ylabel("Count")
plt.show()

# Print the results
print("Median sequence length: ", median_length)
print("Maximum sequence length: ", max_length)
print("Minimum sequence length: ", min_length)

# Print the results
print("Median sequence length: ", median_length)
print("Maximum sequence length: ", max_length)
print("Minimum sequence length: ", min_length)
