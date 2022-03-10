import csv
from collections import defaultdict

planetary_data = {}
with open("planets.csv", 'r') as data_file:
    row1 = data_file.readline()
    row1 = row1.strip().split(",")[1:] # throw away first label
    # print(row1, len(row1))
    for row in data_file:
        row = row.strip().split(",")
        planetary_data.setdefault(row[0],{})[physical_parameter] = float(row[1])

print(planetary_data)

planetary_data = {}
with open("planets.csv", 'r') as data_file:
    row1 = data_file.readline()
    row1 = row1.strip().split(",")[1:] # throw away first label
    for row in data_file:
        row = row.strip().split(",")
        for cols in range(len(row1)):
            physical_parameter = row1[cols]
            #print(physical_parameter)
            #planetary_data.setdefault(row[0],{})[physical_parameter] = float(row[cols])
            #print(row[0], physical_parameter, row[cols+1])
            planetary_data.setdefault(row[0], {})[physical_parameter] = float(row[cols+1])

print(planetary_data)
