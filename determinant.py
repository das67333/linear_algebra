import csv, time

def det(M):
    M = [row[:] for row in M]
    N, sign, prev = len(M), 1, 1
    for i in range(N-1):
        if M[i][i] == 0:
            swapto = next( (j for j in range(i+1,N) if M[j][i] != 0), None )
            if swapto is None:
                return
            M[i], M[swapto], sign = M[swapto], M[i], -sign
        for j in range(i+1,N):
            for k in range(i+1,N):
                assert ( M[j][k] * M[i][i] - M[j][i] * M[i][k] ) % prev == 0
                M[j][k] = ( M[j][k] * M[i][i] - M[j][i] * M[i][k] ) // prev
        prev = M[i][i]
    return sign * M[-1][-1]

def timer(matrix):
    t1 = time.time()
    d = det(matrix)
    t2 = time.time()
    print(t2 - t1)
    return d

def reduce(matrix, n, m):
    assert(n < len(matrix) and m < len(matrix[0]))
    result = []
    for i in range(n):
        result.append(matrix[i][:m])
    return result

m = []
with open("matrix.csv") as file:
    reader = csv.reader(file)
    for i in reader:
        m.append(list(map(int, i)))
# for i in range(100, 500, 50):
print(timer(reduce(m, 150, 150)))
