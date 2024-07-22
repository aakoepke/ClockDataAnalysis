######## Old C++ code with new C++ FFT implementation Comparison ##########

N=5000
t.vec <- 1:N  # Time vector

## tapers
setW = 8/N
setK = 5
taperObject <- get_tapers(t.vec, W = setW, K = setK)
taperMat <- taperObject$tapers


tN = 10
R_mat <- toeplitz(c(seq(1,0.1, length.out = tN), rep(0, times = N-tN))) #to start
c_vec <- c(R_mat[,1], 0, rev(R_mat[1,2:N]))
length(c_vec)

freq <- seq(0, .5, length.out = N)  # Frequency vector

V_star_mat <- t(taperMat*exp(1i*2*pi*freq[3]*t.vec))
V_mat <- taperMat*exp(-1i*2*pi*freq[505]*t.vec)

sourceCpp('est_entry_FFT.cpp')
sourceCpp('est_entry.cpp')

cpp_start_time = Sys.time()
print("Starting Cpp function")
test_entry = est_entry_FFT(V_star_mat, V_mat, c_vec, K, N)
print(Sys.time()-cpp_start_time)

print(test_entry)


cpp_start_time = Sys.time()
print("Starting Cpp function")
test_entry = est_entry(V_star_mat, V_mat, R_mat, K)
print(Sys.time()-cpp_start_time)

print(test_entry)

