s = 5
t = 6
r = 7



reduction_string(G1,bg(8,t,s,2,u=u)*u[r,3]+u[2,t]*u[4,s]*rinj(3,r,u=u))




# u[i,a]*u[j,b]*u[k,c]*u[m,d]

# graded-4 component
gbrep_gr4 = (u[s,2]*bg(8,3,t,r,u=u)+rinj(s,2,u=u)*u[4,t]*u[3,r]
+sum(rinj(s,2,u=u)*u[i,t]*u[3,r] for i in (5:n))
-sum(rinj(s,3,u=u)*u[i,t]*u[2,r] for i in (4:n))
-sum(rinj(s,j,u=u)*u[i,t]*u[m,r] for m in (4:n) for i in (2:n) for j in (3:n) if i != 3 && j != s)
-sum(rinj(s,j,u=u)*u[3,t]*u[m,r] for m in (2:n) for j in (4:n) if  m != 3 && j!= s)
-sum(rinj(s,j,u=u)*u[i,t]*u[2,r] for i in (4:n) for j in (4:n) if s!= j)
-sum(inj(s,a,b,u=u)*u[c,t]*u[3,r] for a in (3:n) for b in (2:n) for c in (2:n) if s!=b)
-sum(u[s,2]*inj(b,t,c,u=u)*u[3,r] for b in (2:n) for c in (2:n) if c != 3 && b != c)
-sum(u[s,a]*inj(b,t,c,u=u)*u[3,r] for a in (3:n) for b in (2:n) for c in (2:n) if c != 3 && b != c && t!=a)
-sum(u[s,a]*u[j,2]*col_sum(t,u)*u[3,r] for a in (3:n) for j in (3:n))
-sum(u[s,a]*u[j,b]*col_sum(t,u)*u[3,r] for a in (2:n) for j in (3:n) for b in (3:n) if b != t && a!=b)
+sum(u[s,a]*u[i,b]*u[j,t]*col_sum(r,u) for a in (2:n) for i in (3:n) for b in (2:n) for j in (3:n) if i!=s && (a,b)!=(2,2))
+sum(wel(s,a,b,u=u)*u[i,t]*u[3,r] for a in (2:n) for i in (2:n) for b in (2:n) if b != a)
-sum(u[s,2]*wel(3,b,t,u=u)*u[m,r] for b in (3:n) for m in (2:n) if b != t && m != 3)
-sum(u[s,a]*wel(3,b,t,u=u)*u[m,r] for a in (3:n) for b in (2:n) for m in (2:n) if b != t && m != 3)
+sum(u[s,a]*wel(2,a,t,u=u)*u[3,r] for a in (3:n))
-sum(u[s,a]*u[b,t]*wel(3,t,r,u=u) for a in (3:n) for b in (4:n))
+sum(u[s,a]*u[2,a]*wel(3,t,r,u=u) for a in (3:n))
+sum(u[s,t]*u[b,t]*wel(3,t,r,u=u) for b in (4:n) if b!=s)
-sum(u[s,2]*u[b,t]*wel(3,t,r,u=u) for b in (4:n))
+sum(u[s,2]*u[j,b]*rwel(t,r,u=u) for j in (3:n) for b in (3:n) if j!= s)
+sum(u[s,a]*u[j,b]*rwel(t,r,u=u) for a in (3:n) for j in (3:n) for b in (2:n) if j!= s)
+sum(u[s,2]*bg(8,j,t,r,u=u) for j in (4:n) if j!= t)
+sum(u[s,a]*bg(8,j,t,r,u=u) for a in (3:n) for j in (2:n) if j!=a && j!=t)
+ip(s,t,u=u)*u[3,t]*u[3,r]
-sum(u[s,a]*ip(3,t,u=u)*u[d,r] for a in (2:n) for d in (2:n))
+sum(u[s,t]*ip(b,t,u=u)*u[3,r] for b in (3:n))
-sum(u[s,a]*ip(b,t,u=u)*u[3,r] for a in (2:n) for b in (4:n))
+sum(u[s,a]*row_sum(2,u)*u[k,t]*u[3,r] for a in (3:n) for k in (4:n))
-sum(u[s,a]*row_sum(j,u)*u[k,t]*u[2,r] for a in (3:n) for j in (3:n) for k in (3:n) if j !=s && (j,k)!=(3,3))
-sum(u[s,a]*row_sum(j,u)*u[k,t]*u[m,r] for a in (3:n) for j in (3:n) for k in (2:n) for m in (4:n) if j != s && (j,k)!=(3,3)) 
)


reduction_string(G0,gbrep_gr4
-sum(row_sum(2,u)*u[b,t]*u[3,r] for b in (4:n))
+sum(row_sum(3,u)*u[b,t]*u[c,r] for b in (2:n) for c in (4:n) if b != 3)
+sum(row_sum(3,u)*u[b,t]*u[2,r] for b in (4:n) if b != 3)
+sum(row_sum(a,u)*u[b,t]*u[c,r] for a in (4:n) for b in (2:n) for c in (4:n) if a!= s)
+sum(row_sum(a,u)*u[b,t]*u[2,r] for a in (4:n) for b in (3:n) if a!= s)
+sum(u[s,j]*row_sum(2,u)*u[i,r] for j in (3:n) for i in (4:n))
-sum(u[s,3]*row_sum(a,u)*u[3,r] for a in (3:n))
+sum(u[s,j]*row_sum(a,u)*u[i,t] for j in (3:n) for a in (3:n) for i in (2:n) if (a,i)!=(3,3) && a!=s) # >6k step
+2*sum(u[s,j]*row_sum(a,u)*u[i,r] for j in (3:n) for a in (3:n) for i in (2:n) if a!=s) # >4k step
-sum(bg(8,a,t,r,u=u) for a in (2:n) if a != t)   ## from here on G0 suffices
+sum(u[a,2]*col_sum(t,u)*u[3,r] for a in (3:n))
+sum(u[a,3]*col_sum(t,u)*u[3,r] for a in (3:n))
+5*sum(u[s,a]*col_sum(t,u)*u[3,r] for a in (2:n))
-sum(u[a,i]*u[b,t]*col_sum(r,u) for a in (3:n) for i in (2:n) for b in (3:n)) # 3k step
+sum(u[s,2]*u[s,j]*col_sum(r,u) for j in (3:n))
-2*sum(u[s,a]*u[b,j]*col_sum(r,u) for a in (2:n) for b in (3:n) for j in (2:n) if (j,a)!=(2,2))
-4*sum(u[s,j]*u[b,t]*col_sum(r,u) for j in (2:n) for b in (3:n))
+sum(u[s,2]*u[b,r]*col_sum(r,u) for b in (3:n))
-sum(u[b,i]*rwel(t,r,u=u) for b in (3:n) for i in (2:n) if b != s) # 4k step
-sum(u[s,i]*rwel(j,r,u=u) for i in (2:n) for j in (2:n) if j != r && (i,j)!=(2,2)) # 4k step 
-5*sum(u[s,a]*rwel(t,r,u=u) for a in (2:n))
+sum(rinj(s,2,u=u)*u[c,r] for c in (4:n))
+sum(rinj(s,3,u=u)*u[c,t] for c in (2:n) if c!=3)
+sum(rinj(s,a,u=u)*u[b,t] for a in (4:n) for b in (2:n) if a != s) # 6k step
-sum(rinj(s,a,u=u)*u[3,r] for a in (4:n) if a != s)
+2*sum(rinj(s,a,u=u)*u[b,r] for a in (3:n) for b in (2:n) if a != s)
+sum(u[s,2]*wel(3,b,t,u=u) for b in (3:n) if b!=t)
-sum(u[s,2]*wel(3,b,r,u=u) for b in (3:n) if b != r)
+u[s,2]*ip(3,t,u=u)
-u[s,2]*ip(3,r,u=u)
-sum(u[s,j]*ip(a,r,u=u) for j in (3:n) for a in (3:n))
-sum(u[s,j]*inj(a,r,b,u=u) for j in (2:n) for a in (2:n) for b in (2:n) if a!=b)
+sum(inj(s,j,a,u=u)*u[3,r] for j in (3:n) for a in (2:n) if a!= s)
+sum(inj(i,t,a,u=u)*u[3,r] for i in (3:n) for a in (2:n) if a!= i)
+sum(wel(s,j,a,u=u)*u[b,r] for j in (2:n) for a in (3:n) for b in (2:n) if b!=3 && j!=a) #5k
+sum(wel(3,i,t,u=u)*u[c,r] for i in (2:n) for c in (2:n) if c != 3 && i!=t)
-u[s,2]*wel(3,t,r,u=u)
)





