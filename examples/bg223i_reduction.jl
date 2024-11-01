t= 7

reduction_string(G0, bg(2,2,3,t,u=u)+rwel(2,4,u=u)*u[t,3]
                     +sum(bg(2,2,i,t,u=u) for i in (4:n) if i!=t)
                     -sum(u[j,2]*wel(t,i,3,u=u) for i in (4:n) for j in (2:n))
                     +sum(rwel(2,i,u=u)*u[t,3] for i in (5:n))
                     +sum(u[j,k]*col_sum(i,u)*u[t,3] for i in (4:n) for j in (3:n) for k in (2:n))
                     +sum(bg(2,j,i,t,u=u) for i in (2:n) for j in (3:n) if i!=t && i!=j)
                     -sum(wel(i,j,k,u=u)*u[t,3] for i in (3:n) for j in (2:n) for k in (4:n) if j!=k && (i,j)!=(t,2))
                     -sum(wel(i,4,3,u=u)*u[t,3] for i in (2:n) if (i,4)!=(t,2))
                     -sum(wel(i,j,2,u=u)*u[t,3] for i in (3:n) for j in (2:n) if j!= 2 && (i,j)!=(t,2))
                     -sum(u[j,k]*row_sum(i,u)*u[t,3] for i in (2:n) for j in (2:n) for k in (3:n) if i!=t && i!=j)
                     +sum(rwel(3,i,u=u)*u[t,3] for i in (2:n) if i!=3)
                     +sum(rwel(k,i,u=u)*u[t,3] for i in (2:n) for k in (5:n) if i!=k)
                     +sum(inj(k,j,i,u=u)*u[t,3] for i in (2:n) for k in (2:n) for j in (3:n) if i!=k)
                     -sum(u[k,j]*wel(t,i,3,u=u) for i in (2:n) for j in (3:n) for k in (2:n) if i!=3 && k !=t)
                     -sum(u[k,3]*ip(t,3,u=u) for k in (2:n) if k != t)
                     -sum(u[k,j]*ip(t,3,u=u) for k in (2:n) for j in (5:n) )
                     +sum(u[i,j]*col_sum(2,u)*u[t,3] for i in (3:n) for j in (2:n))
                     +sum(u[i,j]*col_sum(3,u)*u[t,3] for i in (3:n) for j in (5:n))
                     -sum(u[i,k]*col_sum(k,u)*u[t,3] for i in (3:n) for k in (2:n) if k!=4 && k!=3)
                     -sum(u[i,4]*col_sum(4,u)*u[t,3] for i in (3:n))
                     +rwel(4,2,u=u)*u[t,3]
                     +sum(rwel(4,i,u=u)*u[t,3] for i in (5:n))
                     +sum(u[i,4]*inj(j,3,t,u=u) for i in (2:n) for j in (2:n) if j!=t)
                     +sum(u[t,i]*ip(t,3,u=u) for i in (4:n))
                     -sum(wel(k,j,3,u=u)*u[t,3] for k in (3:n) for j in (5:n))
                     +(n-2)*sum(row_sum(i,u)*u[t,3] for i in (2:n) if i != t)
                     -(n-2)*sum(col_sum(i,u)*u[t,3] for i in (2:n))
                     +2*col_sum(3,u)*u[t,3]
                     -2*sum(inj(i,3,t,u=u) for i in (2:n) if i!= t)
                     +(n-2)*sum(wel(t,i,3,u=u) for i in (2:n) if i !=3)
                     -wel(t,2,3,u=u)
                     +4*ip(t,3,u=u))


