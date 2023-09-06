using LegendrePolynomials
t_list = [0.036 0.072 0.108 0.1728]
Nz = 101
zgrid = range(0,1,Nz)
C = ones(Nz, length(t_list))
for n=1:50
    C += (2*n+1)*Pl.(zgri)