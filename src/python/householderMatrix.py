import numpy as np

x = np.array([0, 2, 3, 4, 5, 6])
e1 = np.array([0, 1, 0, 0, 0, 0])
sigma = np.linalg.norm(x)
V = np.add(x, np.multiply(sigma, e1))  # X+sigma*e1
V = np.divide(V, np.linalg.norm(V))
H = np.identity(6)
H = np.subtract(H, np.multiply(2, np.outer(V, V)))
print(H)
# print(np.dot(H, np.transpose(x)))
print(np.dot(H, x))
