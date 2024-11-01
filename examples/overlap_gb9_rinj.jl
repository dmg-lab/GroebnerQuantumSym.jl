s = 5
t = 6
r = 7



reduction_string(G1,bg(8,t,s,2,u=u)*u[r,3]+u[2,t]*u[4,s]*rinj(3,r,u=u))






reduction_string(G1,u[s,2]*bg(8,3,t,r,u=u)+rinj(s,2,u=u)*u[4,t]*u[3,r]
+sum(rinj(s,2,u=u)*u[i,t]*u[3,r] for i in (5:n))
-sum(rinj(s,3,u=u)*u[i,t]*u[2,r] for i in (4:n))
-sum(rinj(s,k,u=u)*u[i,t]*u[j,r] for j in (4:n) for i in (2:n) for k in (3:n) if i != 3 && k != s)
-sum(rinj(s,k,u=u)*u[3,t]*u[j,r] for j in (2:n) for k in (4:n) if  j != 3 && k!= s)
-sum(rinj(s,k,u=u)*u[i,t]*u[2,r] for i in (5:n) for k in (4:n) if s!= k)
-sum(u[s,2]*u[i,3]*col_sum(t,u)*u[3,r] for i in (3:n))
-sum(u[s,2]*wel(3,3,t,u=u)*u[i,r] for i in (2:n) if i != 3)
+sum(u[s,2]*u[j,3]*u[i,t]*col_sum(r,u) for i in (3:n) for j in (3:n))
-sum(u[s,2]*u[5,3]*u[i,t]*col_sum(r,u) for i in (3:n))
+sum(u[s,2]*u[i,3]*rwel(t,r,u=u) for i in (3:4))
+sum(wel(s,2,3,u=u)*u[i,t]*u[3,r] for i in (3:n))
)





