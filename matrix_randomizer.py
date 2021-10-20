import csv, random

n = 1000
min_value, max_value = -9, 9

with open('matrix.csv', 'w') as File:
    writer = csv.writer(File, delimiter=',', lineterminator = '\n')
    for _ in range(n):
        writer.writerow((random.randint(min_value, max_value) for i in range(n)))
