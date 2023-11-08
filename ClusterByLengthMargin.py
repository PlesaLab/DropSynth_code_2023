import sys
#import matplotlib.pyplot as plt

visualize = False

if len(sys.argv) != 5 and len(sys.argv) != 6:
	print("ARGUMENTS NEEDED: input_file output_file min_length max_length [--visual]")
	exit()

if len(sys.argv) == 6 and sys.argv[5] != "--visual":
	print("ARGUMENTS NEEDED: input_file output_file min_length max_length [--visual]")
	print("if you are using --visual, it must be the last arguement")
	exit()

if len(sys.argv) == 6:
	visualize = True

input_file = str(sys.argv[1])
output_file = str(sys.argv[2])
min_length = int(sys.argv[3])
max_length = int(sys.argv[4])

in_file = open(input_file,"r")
out_file = open(output_file,"w")

text = in_file.readlines()

lengths = []
count = 0
identifier = ''
for t in text:
	stripped = t.strip()
	if t[0] == '>':
		identifier = t
	elif min_length <= len(stripped) and len(stripped) <= max_length:
		out_file.write(identifier)
		out_file.write(t)
		count += 1
	if visualize and t[0] != '>':
		# Add sequence length to plot
		lengths.append(len(stripped))


out_file.close()
print("NUMBER OF SEQUENCES OUTPUT: ", count)
print("Min Length: ", str(min(lengths)))
print("Max Length: ", str(max(lengths)))

# Visual
#if visualize:
#	# make a histogram with bins starting at min ending at max with intervals of 10
#	plt.ylabel('Frequency')
#	plt.xlabel('Sequence Length')
#	plt.title('Frequency of Sequence Lengths Between 400 and 480 and Sequence Identity Threshold of 85%')
#	plt.hist(lengths, bins = range(min(lengths), max(lengths) + 1, 2))
#	plt.show()
