using Oscar
n = 8
G1 = g1_named(n)
u = magic_unitary(n)

G0 = g0_named(n)

showall(reduction_string(G0, bg(9,3,4,u=u)
                     +rinj(3,2,u=u)*sum(u[4:n,4])
                     +sum([rinj(3,i,u=u)*sum(u[2:n,4]) for i in (4,n)])
                     -u[3,2]*sum([u[i,j] for i in (3:n) for j in (3:n) if j != 4])*col_sum(4,u)
                     +wel(3,2,3,u=u)*sum(u[2:n,4])
                     -u[3,2]*sum([inj(i,4,j,u=u) for i in (2:n) for j in (2:n) if i != j])
                     +sum([wel(3,2,i,u=u)*u[j,4] for i in (4:n) for j in (2:n)])
                     -u[3,2]*sum([ip(i,4,u=u) for i in (3:n)])
                     -u[3,2]*sum([rwel(i,4,u=u) for i in (5:n)])
                     +sum([u[3,i]*row_sum(j,u)*u[k,4] for i in (3:n) for j in (2:n) for k in (3:n) if j!=3])
                     +u[3,3]*sum([row_sum(j,u) for j in (4:n)])*u[2,4]
                     -u[3,3]*sum(u[3:n,2])*col_sum(4,u)
                     +sum([rinj(3,k,u=u)*u[i,4] for i in (2:n) for k in (5:n-1)])
                     -u[3,3]*sum([rwel(i,4,u=u) for i in (2:n) if i!=3 && i!=4])
                     +sum([wel(3,j,2,u=u)*u[i,4] for i in (2:n) for j in (3:n) if j!= 4 && j!= 5])
                     -sum([inj(3,3,j,u=u)*u[i,4] for i in (2:n) for j in (2:n) if j !=3 && (i,j)!=(2,2)])
                     -u[3,3]*sum([inj(j,4,i,u=u) for i in (2:n) for j in (2:n) if i != j])
                     +sum([wel(3,3,j,u=u)*u[i,4] for i in (2:n) for j in (4:n)])
                     -u[3,3]*sum([ip(i,4,u=u) for i in (3:n)])
                     -u[3,3]*sum([u[i,j] for i in (3:n) for j in (5:n)])*col_sum(4,u)
                     +u[3,4]*sum([row_sum(i,u) for i in (4:n)])*u[2,4]
                     -sum([u[3,k]*rwel(i,4,u=u) for i in (2:n) for k in (4:n) if i != 4])
                     -sum([u[3,k]*sum(u[3:n,j])*col_sum(4,u) for j in (2:n) for k in (4:n) if j != 4])
                     +sum([wel(3,4,j,u=u)*sum(u[2:n,4]) for j in (2:3)])
                     -u[3,4]*sum([inj(2,4,i,u=u) for i in (3:n)])
                     -sum(inj(3,4,i,u=u)*u[j,4] for i in (4:n) for j in (2:n) if (i,j) != (4,4))
                     -u[3,4]*ip(4,4,u=u)
                     +sum([wel(3,4,j,u=u)*u[i,4] for j in (5:n) for i in (2:n)])
                     -sum([u[3,k]*inj(j,4,i,u=u) for i in (2:n) for j in (2:n) for k in (5:n) if j!=3 && i!=j])
                     -sum([u[3,j]*ip(i,4,u=u) for i in (4:n) for j in (5:n)])
                     +sum([ip(3,j,u=u)*u[i,4] for i in (2:n) for j in (5:n)])
                     +sum([wel(3,5,2,u=u)*u[i,4] for i in (2:n)])
                     +sum([wel(3,k,j,u=u)*u[i,4] for i in (2:n) for j in (3:n) for k in (5:n) if j!=4 && j!=k])
                     +sum([u[3,j]*row_sum(i,u) for i in (4:n) for j in (5:n)])*u[2,4]
                     -sum([row_sum(i,u)*u[j,4] for i in (2:n) for j in (3:n) if i != 3])
                     -sum([row_sum(i,u)*u[2,4] for i in (4:n)])
                     +sum([rwel(j,4,u=u) for j in (2:n) if j != 4])
                     -sum([rinj(3,j,u=u) for j in (2:n) if j != 3])
                     +sum([u[j,2]*col_sum(4,u) for j in (3:n)])
                     +(n-3)*sum([u[3,j]*col_sum(4,u) for j in (2:n)])
                     -sum([u[3,3]*row_sum(j,u) for j in (2:n) if j!=3])
                     -sum([wel(3,k,j,u=u) for j in (2:n) for k in (2:n) if j != k])
                     +sum([inj(3,3,j,u=u) for j in (2:n) if j != 3])
                     +sum([u[j,k]*col_sum(4,u) for j in (3:n) for k in (3:n) if k!= 4])
                     +sum([inj(k,4,j,u=u) for j in (2:n) for k in (2:n) if j != k && (k,j) != (3,2) && (k,j) != (3,4)])
                     -sum([u[3,k]*row_sum(j,u) for j in (2:n) for k in (4:n) if j!=3])
                     +sum([ip(k,4,u=u) for k in (4:n)])
                     -sum([ip(3,k,u=u) for k in (5:n)])
                     +sum([u[3,k]*col_sum(4,u) for k in (4:n)])
                     +sum([wel(3,k,4,u=u) for k in (5:n)])
                     +sum([row_sum(k,u) for k in (2:n) if k !=3])
                     -(n-2)*col_sum(4,u)))



s = 5
t = 6

reduction_string(G0, bg(9,s,t,u=u)
                     +rinj(s,2,u=u)*sum(u[4:n,t])
                     +sum([rinj(s,i,u=u)*sum(u[2:n,t]) for i in (3,n) if i!= s])
                     -u[s,2]*sum([u[i,j] for i in (3:n) for j in (3:n) if j != t])*col_sum(t,u)
                     +wel(s,2,3,u=u)*sum(u[2:n,t])
                     +sum([rinj(s,j,u=u)*u[i,t] for i in (2:n) for j in (4:n-1) if j!=s])
                     +sum([wel(s,2,i,u=u)*u[j,t] for i in (4:n) for j in (2:n)])
                     -u[s,2]*sum([inj(i,t,j,u=u) for i in (2:n) for j in (2:n) if i != j])
                     -u[s,2]*sum([ip(i,t,u=u) for i in (3:n)])
                     -u[s,2]*sum([rwel(i,t,u=u) for i in (4:n) if i!=t])
                     +sum([u[s,i]*row_sum(j,u)*u[k,t] for i in (3:n) for j in (2:n) for k in (3:n) if j!=s])
                     +u[s,3]*sum([row_sum(j,u) for j in (3:n) if j!=s])*u[2,t]
                     -u[s,3]*sum(u[3:n,2])*col_sum(t,u)
                     -u[s,3]*sum([rwel(i,t,u=u) for i in (2:n) if i!=s && i!=t])
                     +sum([wel(s,j,2,u=u)*u[i,t] for i in (2:n) for j in (3:n) if j!= s && j!= t])
                     -sum([u[s,k]*sum(u[3:n,j])*col_sum(t,u) for j in (3:n) for k in (3:n) if j != t])
                     +sum([ip(s,3,u=u)*u[i,t] for i in (2:n)])
                     +sum([wel(s,3,j,u=u)*u[i,t] for i in (2:n) for j in (4:n)])
                     -u[s,3]*sum([inj(j,t,i,u=u) for i in (2:n) for j in (2:n) if i != j])
                     -u[s,3]*sum([ip(i,t,u=u) for i in (3:n)])
                     +u[s,4]*sum([row_sum(i,u) for i in (4:n)])*u[2,t]

############
                     +sum([rinj(3,k,u=u)*u[i,4] for i in (2:n) for k in (5:n-1)])
                     -sum([inj(3,3,j,u=u)*u[i,4] for i in (2:n) for j in (2:n) if j !=3 && (i,j)!=(2,2)])
                     -u[3,3]*sum([u[i,j] for i in (3:n) for j in (5:n)])*col_sum(4,u)
                     -sum([u[3,k]*rwel(i,4,u=u) for i in (2:n) for k in (4:n) if i != 4])
                     +sum([wel(3,4,j,u=u)*sum(u[2:n,4]) for j in (2:3)])
                     -u[3,4]*sum([inj(2,4,i,u=u) for i in (3:n)])
                     -sum(inj(3,4,i,u=u)*u[j,4] for i in (4:n) for j in (2:n) if (i,j) != (4,4))
                     -u[3,4]*ip(4,4,u=u)
                     +sum([wel(3,4,j,u=u)*u[i,4] for j in (5:n) for i in (2:n)])
                     -sum([u[3,k]*inj(j,4,i,u=u) for i in (2:n) for j in (2:n) for k in (5:n) if j!=3 && i!=j])
                     -sum([u[3,j]*ip(i,4,u=u) for i in (4:n) for j in (5:n)])
                     +sum([ip(3,j,u=u)*u[i,4] for i in (2:n) for j in (5:n)])
                     +sum([wel(3,5,2,u=u)*u[i,4] for i in (2:n)])
                     +sum([wel(3,k,j,u=u)*u[i,4] for i in (2:n) for j in (3:n) for k in (5:n) if j!=4 && j!=k])
                     +sum([u[3,j]*row_sum(i,u) for i in (4:n) for j in (5:n)])*u[2,4]
                     -sum([row_sum(i,u)*u[j,4] for i in (2:n) for j in (3:n) if i != 3])
                     -sum([row_sum(i,u)*u[2,4] for i in (4:n)])
                     +sum([rwel(j,4,u=u) for j in (2:n) if j != 4])
                     -sum([rinj(3,j,u=u) for j in (2:n) if j != 3])
                     +sum([u[j,2]*col_sum(4,u) for j in (3:n)])
                     +(n-3)*sum([u[3,j]*col_sum(4,u) for j in (2:n)])
                     -sum([u[3,3]*row_sum(j,u) for j in (2:n) if j!=3])
                     -sum([wel(3,k,j,u=u) for j in (2:n) for k in (2:n) if j != k])
                     +sum([inj(3,3,j,u=u) for j in (2:n) if j != 3])
                     +sum([u[j,k]*col_sum(4,u) for j in (3:n) for k in (3:n) if k!= 4])
                     +sum([inj(k,4,j,u=u) for j in (2:n) for k in (2:n) if j != k && (k,j) != (3,2) && (k,j) != (3,4)])
                     -sum([u[3,k]*row_sum(j,u) for j in (2:n) for k in (4:n) if j!=3])
                     +sum([ip(k,4,u=u) for k in (4:n)])
                     -sum([ip(3,k,u=u) for k in (5:n)])
                     +sum([u[3,k]*col_sum(4,u) for k in (4:n)])
                     +sum([wel(3,k,4,u=u) for k in (5:n)])
                     +sum([row_sum(k,u) for k in (2:n) if k !=3])
                     -(n-2)*col_sum(4,u)))











