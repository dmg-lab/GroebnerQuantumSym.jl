using Oscar
n = 8
G1 = g1_named(n)
u = magic_unitary(n)

G0 = g0_named(n)


### dated:: (see below for new rep)
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
#s=3,4,5,6,....,  t=2 -> "almost" zero, inj_{s,2,2} remains
#s=3,4,5,6,. , t>= 4  -> ok

##########....
##########....
#a#ooooooo....
#a#ooooooo....
#a#ooooooo....
#a#ooooooo....
#a#ooooooo....
#a#ooooooo....


s = 5
t = 6

reduction_string(G0, bg(9,s,t,u=u)
+sum(rinj(s,j,u=u)*u[i,t] for i in (2:n) for j in (2:n) if j!= s && (j!=2 || (j==2 && i!=2 && i!=3)))
-sum(u[s,k]*u[i,j] for i in (3:n) for j in (2:n) for k in (2:n) if j != t && (k!=2||(k==2 && j!=2)))*col_sum(t,u)
+sum(wel(s,k,j,u=u)*u[i,t] for i in (2:n) for j in (2:n) for k in (2:n) if (j!=t || k==3 || k==2) && j!=k)
-sum(u[s,k]*inj(j,t,i,u=u) for i in (2:n) for j in (2:n) for k in (2:n) if i!=j && ((k!=2 && k!=3 && j!=s)|| (k==2 && t!=2) || k==3))
-sum(u[s,j]*ip(i,t,u=u) for i in (3:n) for j in (2:n) if j==3 || (j!=3 && j!=2 && i !=s) || (j==2 && t!=2))
-sum([u[s,k]*rwel(i,t,u=u) for i in (2:n) for k in (2:n) if i != t && (k!=2 || (k==2 && i != 2 && i!=3))])
+sum([u[s,i]*row_sum(j,u)*u[k,t] for i in (3:n) for j in (2:n) for k in (2:n) if j!=s && (k,j)!=(2,2)])
+sum([ip(s,j,u=u)*u[i,t] for i in (2:n) for j in (3:n) if j!=t])
-sum(row_sum(i,u)*u[j,t] for i in (2:n) for j in (2:n) if i !=s &&(j!=2||(j==2 && i!=2)))
+sum([rwel(j,t,u=u) for j in (2:n) if j != t])
-sum([rinj(s,j,u=u) for j in (2:n) if j != s])
-sum([wel(s,k,j,u=u) for j in (2:n) for k in (2:n) if j != k && (j!=t || k ==2 || k==3)])
+sum([inj(k,t,j,u=u) for j in (2:n) for k in (2:n) if j != k && k !=s])
+sum([u[j,k]*col_sum(t,u) for j in (3:n) for k in (2:n) if k!=t])
+sum([u[s,k]*col_sum(t,u) for k in (3:n)])
+(n-3)*sum([u[s,k]*col_sum(t,u) for k in (2:n)])
-sum([u[s,k]*row_sum(j,u) for j in (2:n) for k in (3:n) if j != s])
-sum([ip(s,k,u=u) for k in (3:n)])
+sum([ip(k,t,u=u) for k in (3:n) if (k,t)!=(s,2)])
+sum([row_sum(k,u) for k in (2:n) if k !=s])
-(n-2)*col_sum(t,u))


#### dated: see above
##### before cleaning : 
#initial reduction string has 35.000 bytes, with this it is already zero (tested for n<=9:
reduction_string(G0, bg(9,s,t,u=u)
                     +rinj(s,2,u=u)*sum(u[4:n,t])
                     +sum([rinj(s,i,u=u)*sum(u[2:n,t]) for i in (3,n) if i!= s])
                     -u[s,2]*sum([u[i,j] for i in (3:n) for j in (3:n) if j != t])*col_sum(t,u)
                     +wel(s,2,3,u=u)*sum(u[2:n,t])
                     +sum([rinj(s,j,u=u)*u[i,t] for i in (2:n) for j in (4:n-1) if j!=s])
                     +sum([wel(s,2,i,u=u)*u[j,t] for i in (4:n) for j in (2:n)])
                     -u[s,2]*sum(append!([inj(i,t,j,u=u) for i in (2:n) for j in (2:n) if i != j && t!=2],[0*ip(1,1)]))
                     -u[s,2]*sum(append!([ip(i,t,u=u) for i in (3:n) if t!=2],[0*ip(1,1)]))
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
                     +sum([u[j,k]*col_sum(t,u) for j in (3:n) for k in (3:n)])
                     -sum([u[s,3]*row_sum(j,u) for j in (2:n) if j!=s])
                     -sum([ip(s,k,u=u) for k in (3:n)])
                     +sum([u[s,k]*col_sum(t,u) for k in (3:n) if k !=t])
                     -sum([u[s,k]*row_sum(j,u) for j in (2:n) for k in (4:n) if j != s])
                     +sum([wel(s,k,t,u=u) for k in (4:n) if k != t])
                     +sum([inj(k,t,j,u=u) for j in (2:n) for k in (2:n) if j != k && k !=s])
                     -sum([u[k,t]*col_sum(t,u) for k in (3:n) if k != s])
                     +sum([ip(k,t,u=u) for k in (3:n) if (k,t)!=(s,2)])
                     +sum([row_sum(k,u) for k in (2:n) if k !=s])
                     -(n-2)*col_sum(t,u))

#### dated: 
#### all leading monomials are distinct from each other and smaller theen lm(bg(9,s,t))
vec_of_leading_mons=[lm(rinj(s,2,u=u)*sum(u[4:n,t])),
                     lm(sum([rinj(s,i,u=u)*sum(u[2:n,t]) for i in (3,n) if i!= s])),
                     lm(u[s,2]*sum([u[i,j] for i in (3:n) for j in (3:n) if j != t])*col_sum(t,u)),
                     lm(wel(s,2,3,u=u)*sum(u[2:n,t])),
                     lm(sum([rinj(s,j,u=u)*u[i,t] for i in (2:n) for j in (4:n-1) if j!=s])),
                     lm(sum([wel(s,2,i,u=u)*u[j,t] for i in (4:n) for j in (2:n)])),
                     lm(u[s,2]*sum(append!([inj(i,t,j,u=u) for i in (2:n) for j in (2:n) if i != j && t!=2],[0*ip(1,1)]))),
                     lm(u[s,2]*sum(append!([ip(i,t,u=u) for i in (3:n) if t!=2],[0*ip(1,1)]))),
                     lm(u[s,2]*sum([rwel(i,t,u=u) for i in (4:n) if i!=t])),
                     lm(sum([u[s,i]*row_sum(j,u)*u[k,t] for i in (3:n) for j in (2:n) for k in (3:n) if j!=s])),
                     lm(u[s,3]*sum([row_sum(j,u) for j in (3:n) if j!=s])*u[2,t]),
                     lm(u[s,3]*sum([rwel(i,t,u=u) for i in (2:n) if i!=t])),
                     lm(sum([wel(s,3,j,u=u)*u[i,t] for i in (2:n) for j in (2:n) if j !=3])),
                     lm(sum([wel(s,k,j,u=u)*u[i,t] for i in (2:n) for j in (2:n) for k in (4:n) if j!=t && j!=k])),
                     lm(sum([u[s,k]*u[i,j] for i in (3:n) for j in (2:n) for k in (3:n) if j != t])*col_sum(t,u)),
                     lm(sum([ip(s,3,u=u)*u[i,t] for i in (2:n)])),
                     lm(u[s,3]*sum([inj(j,t,i,u=u) for i in (2:n) for j in (2:n) if i != j])),
                     lm(u[s,3]*sum([ip(i,t,u=u) for i in (3:n)])),
                     lm(sum([u[s,j]*row_sum(i,u) for i in (3:n) for j in (4:n) if  i != s])*u[2,t]),
                     lm(sum([u[s,k]*rwel(i,t,u=u) for i in (2:n) for k in (4:n) if i != t])),
                     lm(sum([ip(s,j,u=u)*u[i,t] for i in (2:n) for j in (4:n) if j!=t])),
                     lm(sum([u[s,k]*inj(j,t,i,u=u) for i in (2:n) for j in (2:n) for k in (4:n) if j!=s && i!=j])),
                     lm(sum([u[s,j]*ip(i,t,u=u) for i in (3:n) for j in (4:n) if i !=s])),
                     lm(sum([row_sum(i,u)*u[j,t] for i in (2:n) for j in (3:n) if i !=s])),
                     lm(sum([row_sum(i,u)*u[2,t] for i in (3:n) if i != s])),
                     lm(sum([rwel(j,t,u=u) for j in (2:n) if j != t])),
                     lm(sum([u[j,2]*col_sum(t,u) for j in (3:n)])),
                     lm(sum([rinj(s,j,u=u) for j in (2:n) if j != s])),
                     lm(sum([wel(s,k,j,u=u) for j in (2:n) for k in (2:n) if j != k])),
                     lm((n-3)*sum([u[s,j]*col_sum(t,u) for j in (2:n)])),
                     lm(sum([u[j,k]*col_sum(t,u) for j in (3:n) for k in (3:n)])),
                     lm(sum([u[s,3]*row_sum(j,u) for j in (2:n) if j!=s])),
                     lm(sum([ip(s,k,u=u) for k in (3:n)])),
                     lm(sum([u[s,k]*col_sum(t,u) for k in (3:n) if k !=t])),
                     lm(sum([u[s,k]*row_sum(j,u) for j in (2:n) for k in (4:n) if j != s])),
                     lm(sum([wel(s,k,t,u=u) for k in (4:n) if k != t])),
                     lm(sum([inj(k,t,j,u=u) for j in (2:n) for k in (2:n) if j != k && k !=s])),
                     lm(sum([u[k,t]*col_sum(t,u) for k in (3:n) if k != s])),
                     lm(sum([ip(k,t,u=u) for k in (3:n) if (k,t)!=(s,2)])),
                     lm(sum([row_sum(k,u) for k in (2:n) if k !=s])),
                     lm((n-2)*col_sum(t,u))]




n_m1 = n - 1
G1_m1 = g1_named(n_m1)
u_m1 = magic_unitary(n_m1)
G0_m1 = g0_named(n_m1)



################### try to prove eq directly.

# projection on degree 3 of bg_9(s,t)
phi3_bg9st = sum(u[s,2]*u[2,b]*u[3,t] for b in (4:n))-
sum(u[s,a]*u[2,1]*u[3,t] for a in (3:n))-
sum(u[s,2]*u[2,3]*u[k,t] for k in (4:n))+
sum(u[s,2]*u[j,3]*u[1,t] for j in (3:n))
## test 
#bg(9,s,t,u=u)-phi3_bg9st
#u₂₁*u₃₆ + u₅₂*u₂₃ - u₅₂*u₁₆ - u₅₂*u₃₆
#i.e., no deg 3 terms


# projection on degree 3 of term1
phi3_term1 = sum(u[s,2]*u[j,b]*u[k,t] for j in (2:n) for b in (3:n) for k in (2:n) if j!= s && (j!=2 || (j==2 && k!=2 && k!=3)))-
sum(u[s,a]*u[j,1]*u[k,t] for a in (3:n) for j in (2:n) for k in (2:n) if j!= s && (j!=2 || (j==2 && k!=2 && k!=3)))
## test 
#sum(rinj(s,j,u=u)*u[i,t] for i in (2:n) for j in (2:n) if j!= s && (j!=2 || (j==2 && i!=2 && i!=3)))-phi3_term1
#---> no deg 3 terms


# projection on degree 3 of term2
phi3_term2 = -sum(u[s,a]*u[j,b]*u[k,t] for a in (2:n) for j in (3:n) for b in (2:n) for k in (1:n) if b != t && (a!=2||(a==2 && b!=2)))
## test
-sum(u[s,k]*u[i,j] for i in (3:n) for j in (2:n) for k in (2:n) if j != t && (k!=2||(k==2 && j!=2)))*col_sum(t,u)-phi3_term2
#---> no deg 3 terms


# projection on degree 3 of term2
phi3_term3 = 

+sum(wel(s,k,j,u=u)*u[i,t] for i in (2:n) for j in (2:n) for k in (2:n) if (j!=t || k==3 || k==2) && j!=k)
-sum(u[s,k]*inj(j,t,i,u=u) for i in (2:n) for j in (2:n) for k in (2:n) if i!=j && ((k!=2 && k!=3 && j!=s)|| (k==2 && t!=2) || k==3))
-sum(u[s,j]*ip(i,t,u=u) for i in (3:n) for j in (2:n) if j==3 || (j!=3 && j!=2 && i !=s) || (j==2 && t!=2))
-sum([u[s,k]*rwel(i,t,u=u) for i in (2:n) for k in (2:n) if i != t && (k!=2 || (k==2 && i != 2 && i!=3))])
+sum([u[s,i]*row_sum(j,u)*u[k,t] for i in (3:n) for j in (2:n) for k in (2:n) if j!=s && (k,j)!=(2,2)])
+sum([ip(s,j,u=u)*u[i,t] for i in (2:n) for j in (3:n) if j!=t])



