# import numpy as np
# import matplotlib.pyplot as plt

# # Load the data from the new file
# absmL12 = np.loadtxt('/Users/wantake/Desktop/code/Model/SSEAll/SSE-M/L12n/absm.dat')
# m2L12 = np.loadtxt('/Users/wantake/Desktop/code/Model/SSEAll/SSE-M/L12n/m2.dat')
# # Reshape the data into a 2D array with 10 columns (each column represents a group)
# rabsmL12 = absmL12.reshape(-1, 30)
# rm2L12 = m2L12.reshape(-1, 30)
# # Save the reshaped data to a new .dat file
# # np.savetxt('reshaped_data.dat', rabsm, fmt='%f', delimiter=' ')

# # Calculate mean and standard deviation across the new data sets
# mabsmL12 = np.mean(rabsmL12, axis=1)
# eabsmL12 = np.std(rabsmL12, axis=1)

# mm2L12 = np.mean(rm2L12, axis=1)
# em2L12 = np.std(rm2L12, axis=1)



# # Load the data from the new file
# absmL14 = np.loadtxt('/Users/wantake/Desktop/code/Model/SSEAll/SSE-M/L14/absm.dat')
# m2L14 = np.loadtxt('/Users/wantake/Desktop/code/Model/SSEAll/SSE-M/L14/m2.dat')
# # Reshape the data into a 2D array with 10 columns (each column represents a group)
# rabsmL14 = absmL14.reshape(-1, 30)
# rm2L14 = m2L14.reshape(-1, 30)
# # Save the reshaped data to a new .dat file
# # np.savetxt('reshaped_data.dat', rabsm, fmt='%f', delimiter=' ')

# # Calculate mean and standard deviation across the new data sets
# mabsmL14 = np.mean(rabsmL14, axis=1)
# eabsmL14 = np.std(rabsmL14, axis=1)

# mm2L14 = np.mean(rm2L14, axis=1)
# em2L14 = np.std(rm2L14, axis=1)





# # Load the data from the new file
# absmL16 = np.loadtxt('/Users/wantake/Desktop/code/Model/SSEAll/SSE-M/absm.dat')
# m2L16 = np.loadtxt('/Users/wantake/Desktop/code/Model/SSEAll/SSE-M/m2.dat')
# # Reshape the data into a 2D array with 10 columns (each column represents a group)
# rabsmL16 = absmL16.reshape(-1, 30)
# rm2L16 = m2L16.reshape(-1, 30)
# # Save the reshaped data to a new .dat file
# # np.savetxt('reshaped_data.dat', rabsm, fmt='%f', delimiter=' ')

# # Calculate mean and standard deviation across the new data sets
# mabsmL16 = np.mean(rabsmL16, axis=1)
# eabsmL16 = np.std(rabsmL16, axis=1)

# mm2L16 = np.mean(rm2L16, axis=1)
# em2L16 = np.std(rm2L16, axis=1)




















# Load the temperature data
# T = np.loadtxt('/Users/wantake/Desktop/code/Model/SSEAll/SSE-M/T.dat')


    # Plot the data
# plt.errorbar(T, mabsmL12, yerr=eabsmL12, color='red', fmt='-.', label='L12 absm Error Bars', ecolor='orange', capsize=2)

# plt.errorbar(T, mabsmL14, yerr=eabsmL14, color='red', fmt='--', label='L14 absm Error Bars', ecolor='orange', capsize=2)

# plt.errorbar(T, mabsmL16, yerr=eabsmL16, color='red', fmt='-', label='L16 absm Error Bars', ecolor='orange', capsize=2)

# plt.errorbar(T, mm2L12/12/12, yerr=em2L12/12/12, color='darkred', fmt='-.', label='L12 m^2 Error Bars', ecolor='orange', capsize=2)

# # plt.errorbar(T, mm2L14, yerr=em2L14, color='darkred', fmt='--', label='L14 m^2 Error Bars', ecolor='orange', capsize=2)

# plt.errorbar(T, mm2L16/16/16, yerr=em2L16/16/16, color='darkred', fmt='-', label='L16 m^2 Error Bars', ecolor='orange', capsize=2)






# plt.xlabel('Temperature')
# # plt.ylabel('|m|')
# plt.title('Heisenberg bin=30 m=1e3 i=1e5')
# plt.legend()
# plt.show()




import numpy as np
import matplotlib.pyplot as plt

# Function to load, reshape, calculate mean and std, and return data
def process_data(file_path, reshape_size):
    data = np.loadtxt(file_path)
    reshaped_data = data.reshape(-1, reshape_size)
    mean = np.mean(reshaped_data, axis=1)
    std = np.std(reshaped_data, axis=1)
    return mean, std

# Load and process data for L=12
meanL12, stdL12 = process_data('/Users/wantake/Desktop/code/Model/SSEAll/SSE-M/L12/absm.dat', 30)

# # Load and process data for L=16 (similar to L=12, replace file paths and reshape sizes)
meanL16, stdL16 = process_data('/Users/wantake/Desktop/code/Model/SSEAll/SSE-M/L16/absm.dat', 30)

# # Load and process data for L=24
meanL24, stdL24 = process_data('/Users/wantake/Desktop/code/Model/SSEAll/SSE-M/L24/absm.dat', 30)

# # Load and process data for L=48
# meanL48, stdL48 = process_data('/path/to/L48/datafile', 30)

# # Load and process data for L=64
# meanL64, stdL64 = process_data('/path/to/L64/datafile', 30)

# # Load and process data for L=128
# meanL128, stdL128 = process_data('/path/to/L128/datafile', 30)


# Load the temperature data
T = np.loadtxt('/Users/wantake/Desktop/code/Model/SSEAll/SSE-M/T.dat')




# Plotting (update this part based on your requirements)

plt.errorbar(T, meanL12, yerr=stdL12, color='grey', fmt='-*', label='L12', ecolor='orange', capsize=2)

plt.errorbar(T, meanL16, yerr=stdL16, color='darkgrey', fmt='-o', label='L16', ecolor='orange', capsize=2)

plt.errorbar(T, meanL24, yerr=stdL24, color='black', fmt='-o', label='L24', ecolor='orange', capsize=2)

# plt.errorbar(T, mm2L12/12/12, yerr=em2L12/12/12, color='darkred', fmt='-.', label='L12 m^2 Error Bars', ecolor='orange', capsize=2)

# # plt.errorbar(T, mm2L14, yerr=em2L14, color='darkred', fmt='--', label='L14 m^2 Error Bars', ecolor='orange', capsize=2)

# plt.errorbar(T, mm2L16/16/16, yerr=em2L16/16/16, color='darkred', fmt='-', label='L16 m^2 Error Bars', ecolor='orange', capsize=2)



plt.xlabel('Temperature')
plt.ylabel('|m|')
plt.title('Heisenberg bin=30 m=1e4 i=1e4')
plt.legend()
plt.show()



# Example: plt.errorbar(x_values, meanL12, yerr=stdL12, label='L=12')
# Add similar lines for L=16, L=24, L=48, L=64, and L=128
# plt.legend()
# plt.show()
