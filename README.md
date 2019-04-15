# ANLY550

File to run Prim's algorithm using an adjacency matrix. Computes average MST weight for given number of trials. Inputs from std.in are: 0 num_points, num_trials, dimension

dimension = 0 generates random edge weights from [0,1]

dimension = 2 generates random vertices from [0,1]^2 and computes 2D Euclideand distance

dimension = 3 generates random vertices from [0,1]^3 and computes 3D Euclideand distance

dimension = 4 generates random vertices from [0,1]^4 and computes 4D Euclideand distance


Caching is used, in addition to throwing away edges from complete graphs for n >= 128 using function k(n) = 0.5^(m-6) where n = 2^m.
