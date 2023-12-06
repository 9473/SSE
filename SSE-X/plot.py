# import numpy as np
# import matplotlib.pyplot as plt

# # Load data files
# T = np.loadtxt('/Users/wantake/Desktop/code/Model/SSEAll/SSE-X/T.dat')

# y1 = np.loadtxt('/Users/wantake/Desktop/code/Model/SSEAll/SSE-X/m^2.dat')
# y2 = np.loadtxt('absm.dat')
# y3 = np.loadtxt('asu.dat')
# y4 = np.loadtxt('usu.dat')
# # y5 = np.loadtxt('qean5.dat')
# # y6 = np.loadtxt('qean6.dat')
# # y7 = np.loadtxt('qean7.dat')
# # y8 = np.loadtxt('qean8.dat')
# # y9 = np.loadtxt('qean9.dat')
# # y10 = np.loadtxt('qean10.dat')

# # # Check if all data arrays and temperature array are of the same length
# # if not all(len(arr) == len(Tsta) for arr in [qeasta1, qeasta2, qeasta3, qeasta4, qeasta5]):
# #     raise ValueError("Data arrays and temperature array must have the same length.")

# # Combine data into a single 2D array for easier processing
# # data = np.array([y1, y2, y3, y4, y5, y6, y7, y8, y9, y10])

# # # Calculate mean and standard deviation (error bar) across the 5 data sets for each temperature
# # mean_values = np.mean(data, axis=0)
# # error_bars = np.std(data, axis=0)

# # x = np.loadtxt('/Users/wantake/Desktop/code/week/nematicqea/Tmenatic.dat')
# # ydata = np.loadtxt('/Users/wantake/Desktop/code/ysy/qeaOutput.dat')

# # numbers = 199

# # y=[]
# # for i in range(numbers):
# #     y.append(np.mean(ydata[i, :]))

# # yerry = []
# # for i in range(numbers):
# #     yerry.append((np.std(ydata[i, :])) / ((numbers) ** (1 / 2)))

# # plt.plot(x, y)
# # plt.errorbar(x, y, yerr=yerry,color='blue' ,fmt='-', label='ysy',ecolor='white',capsize=2)

# # Plotting
# # plt.errorbar(Tsta, mean_values, yerr=error_bars,color='skyblue' ,fmt='-', label='Qea of nematic with Error Bars',ecolor='white',capsize=2)
# plt.plot(T,y1,'-o',color='lightgreen',label='<m^2>')
# plt.plot(T,y2,'-o',color='green',label='<|m|>')
# plt.plot(T,y3,'-o',color='darkgreen',label='<asu>')
# plt.plot(T,y4,'-*',color='b',label='<usu>')
# plt.xlabel('Temperature')
# # plt.ylabel('')
# plt.title('4x4 Heisenberg bin=1')
# plt.legend()
# plt.savefig("1.png",dpi=1000,bbox_inches = 'tight')
# plt.show()

import numpy as np
import matplotlib.pyplot as plt

# Load the data from the new file
total_data = np.loadtxt('/Users/wantake/Desktop/code/Model/SSEAll/SSE-X/usu.dat')

# Reshape the data into a 2D array with 10 columns (each column represents a group)
reshaped_data = total_data.reshape(-1, 10)

# Calculate mean and standard deviation across the new data sets
mean_values = np.mean(reshaped_data, axis=1)
error_bars = np.std(reshaped_data, axis=1)

# Load the temperature data
T = np.loadtxt('/Users/wantake/Desktop/code/Model/SSEAll/SSE-X/T.dat')

# Ensure the temperature array matches the size of the mean_values array
if len(T) != len(mean_values):
    print("Error: Temperature data size does not match the mean values data size.")
else:
    # Plot the data
    plt.errorbar(T, mean_values, yerr=error_bars, color='red', fmt='-o', label='uniform X Error Bars', ecolor='orange', capsize=2)
    plt.xlabel('Temperature')
    plt.ylabel('X')
    plt.title('4x4 Heisenberg bin=10 m=1e3 i=1e5')
    plt.legend()
    plt.show()
