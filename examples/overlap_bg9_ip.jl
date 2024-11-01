r = 5
s = 6
t = 7



reduction_string(G1,bg(8,t,s,r,u=u)*u[3,r]+u[2,t]*u[4,s]*ip(3,r,u=u))

reduction_string(G1,bg(8,t,s,r,u=u)*u[3,r]+u[2,t]*u[4,s]*ip(3,r,u=u)
       +sum(u[2,t]*u[i,s]*ip(3,r,u=u) for i in (5:n))
       +sum(u[j,t]*u[i,s]*ip(3,r,u=u) for i in (2:n) for j in (3:n))
       -sum(u[i,t]*col_sum(s,u)*u[3,r]*u[3,r] for i in (3:n))
       +col_sum(s,u)*u[3,r]*u[3,r]
       +sum(u[i,t]*col_sum(s,u)*u[3,r] for i in (3:n))
       -sum(u[i,s]*ip(3,r,u=u) for i in (2:n))
       -sum(u[i,t]*ip(3,r,u=u) for i in (2:n))
       -bg(8,t,s,r,u=u)
       +ip(3,r,u=u)
       -col_sum(s,u)*u[3,r])

reduction_string(G1,u[2,t]*bg(8,t,s,r,u=u)+ip(2,t,u=u)*u[4,s]*u[3,r])
     
