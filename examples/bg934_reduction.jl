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

# general version ??? seems so, unless s,t ? 2

#NOT((s,t)<=(4,4))  -> ok    (to be confirmed, only some examples)
#s = 4 , t >= 5 -> ok   (tbc)
#s=t=4 -> error
#s=4, t=2 -> "almost" zero, 6 lines, almost all length 3 and all inj or ip
#s=3 , t>= 5  -> ok
#s=3 , t= 4  -> error
#s=3, t=2 -> "almost" zero, 6 lines, almost all length 3 and all inj or ip

##########....
##########....
#a#eoooooo....
#a#eoooooo....
#a#aoooooo....
#a#aoooooo....
#a#aoooooo....
#a#aoooooo....


k = 5
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
                     -u[s,3]*sum([rwel(i,t,u=u) for i in (2:n) if i!=t])
                     +sum([wel(s,3,j,u=u)*u[i,t] for i in (2:n) for j in (2:n) if j !=3])
                     +sum([wel(s,k,j,u=u)*u[i,t] for i in (2:n) for j in (2:n) for k in (4:n) if j!=t && j!=k])
                     -sum([u[s,k]*u[i,j] for i in (3:n) for j in (2:n) for k in (3:n) if j != t])*col_sum(t,u)
                     +sum([ip(s,3,u=u)*u[i,t] for i in (2:n)])
                     -u[s,3]*sum([inj(j,t,i,u=u) for i in (2:n) for j in (2:n) if i != j])
                     -u[s,3]*sum([ip(i,t,u=u) for i in (3:n)])
                     +sum([u[s,j]*row_sum(i,u) for i in (3:n) for j in (4:n) if  i != s])*u[2,t]
                     -sum([u[s,k]*rwel(i,t,u=u) for i in (2:n) for k in (4:n) if i != t])
                     +sum([ip(s,j,u=u)*u[i,t] for i in (2:n) for j in (4:n) if j!=t])
                     -sum([u[s,k]*inj(j,t,i,u=u) for i in (2:n) for j in (2:n) for k in (4:n) if j!=s && i!=j])
                     -sum([u[s,j]*ip(i,t,u=u) for i in (3:n) for j in (4:n) if i !=s])
                     -sum([row_sum(i,u)*u[j,t] for i in (2:n) for j in (3:n) if i !=s])
                     -sum([row_sum(i,u)*u[2,t] for i in (3:n) if i != s])
                     +sum([rwel(j,t,u=u) for j in (2:n) if j != t])
                     +sum([u[j,2]*col_sum(t,u) for j in (3:n)])
                     -sum([rinj(s,j,u=u) for j in (2:n) if j != s])
                     -sum([wel(s,k,j,u=u) for j in (2:n) for k in (2:n) if j != k])
                     +(n-3)*sum([u[s,j]*col_sum(t,u) for j in (2:n)])
                     +sum([u[j,k]*col_sum(t,u) for j in (3:n) for k in (3:n) if t!= 4])
                     -sum([u[s,3]*row_sum(j,u) for j in (2:n) if j!=s])
                     -sum([ip(s,k,u=u) for k in (3:n)])
                     +sum([u[s,k]*col_sum(t,u) for k in (3:n) if k !=t])
                     -sum([u[s,k]*row_sum(j,u) for j in (2:n) for k in (4:n) if j != s])
                     +sum([wel(s,k,t,u=u) for k in (4:n) if k != t])
                     +sum([inj(k,t,j,u=u) for j in (2:n) for k in (2:n) if j != k && (k,t)!= (s,t)])
                     -sum([u[k,t]*col_sum(t,u) for k in (3:n) if k != s])
                     +sum([ip(k,t,u=u) for k in (3:n)])
                     +sum([row_sum(k,u) for k in (2:n) if k !=s])
                     -(n-2)*col_sum(t,u))











