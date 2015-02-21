##Dynamic Test calculations
##Author: Zach Riddle
##Code for Bayesian A/B Testing using counts of Data


library(RODBC);

strDBserver = "ESSBODB"
dbhandle <-odbcConnect(strDBserver)
cat(paste("Connected to Database",strDBserver, "\n"))

start_date<-'2014-06-30'
end_date<-'2014-07-15'
test_name<-'%FAQ%'

query<-paste("SELECT CASE WHEN test_grp LIKE '%Default%' THEN 'Default' ELSE test_grp END test_grp,
             COUNT(*) Visitors,
             SUM(CASE WHEN purchaseID IS NOT NULL THEN 1 ELSE 0 END) reserves
             FROM(
             SELECT visid_compound,
             post_evar26 test_grp,
             max(PURCHASEID) purchaseid
             FROM   ESSCLICKSTREAM
             WHERE  date_time >= \'",start_date,"\'
             AND date_time <= \'",end_date,"\'
             AND post_evar26 like \'",test_name,"\'
             AND VISID_COMPOUND IN (SELECT distinct VISID_COMPOUND
             FROM   ESSCLICKSTREAM 
             WHERE  date_time >= \'",start_date,"\' AND date_time<= \'",end_date,"\'
             AND evar26 LIKE \'",test_name,"\'
             GROUP  BY VISID_COMPOUND)
             GROUP BY visid_compound, POST_EVAR26) x
             GROUP BY test_grp
             ORDER BY CASE WHEN test_grp = 'Default' THEN 1 ELSE 2 END",sep="")


inputs <- sqlQuery(dbhandle, query);
inputs;

##create function to evaluate probability
probability_B_beats_A<-function(A_total,B_total,A_conversions,B_conversions){
  total = 0.0
  ##Assign successes and failures for each
  AS<-A_conversions+1
  AF<-A_total-A_conversions+1
  BS<-B_conversions+1
  BF<-B_total-B_conversions+1
  
  for (i in 0:(BS-1)){
    total <- total + (exp(lbeta(AS+i, BF+AF) 
                          - log(BF+i) - lbeta(1+i, BF) - lbeta(AS, AF)))
  }
  
  return(total)
}

t1<-84783
t2<-84794
a<-1757
b<-1818

probability_B_beats_A(t1,t2,a,b)

t1<-inputs[1,2]
a<-inputs[1,3]
  
for(i in 2:length(inputs)){
  t2<-inputs[i,2]
  b<-inputs[i,3]
  cat('Probability of',inputs[i,1],'winning =',probability_B_beats_A(t1,t2,a,b))
}


##For Variance Test
t1<-c(3132, 10840, 18487, 25452, 31373, 39893, 48001,	56509,	64644,	74091,
      84783,  92118,	101409,	109824,	118220,	126032,	134568,	141085)
t2<-c(3163,  10883,	18507,	25495,	31435,	39964,	48031,	56573,	64682,	
      74145,84794,  92120,	101398,	109806,	118188,	126019,	134262,	141066)
a<-c(65,  253,	439,	572,	679,	886,	1058,	1237,	1422,	1605,
     1757,  1842,	2013,	2151,	2305,	2458,	2642,	2781)
b<-c(51,  232,	410,	580,	700,	906,	1090,	1297,	1490,	1687,
     1818,  1913,	2065,	2205,	2338, 2516,	2692,	2835)

for(i in 1:length(a)){
  cat(probability_B_beats_A(t1[i],t2[i],a[i],b[i]),"\n")
}



f<-function(x){
  pbeta(x,232+1,10883-232+1)*dbeta(x,253+1,10840-253+1)
}

integrate(f,0,1)


###################################
### Multi-Armed Bandit Algorithms
###################################
##Computing optimal probabilities for each arm

##Need logbeta probably
compute.probopt<-function(y,n){
  k<-length(y)
  ans<-numeric(k)
  for(i in 1:k){
    indx<-(1:k)[-i]
    f<-function(x){
      r<-dbeta(x,y[i]+1,n[i]-y[i]+1)
      for(j in indx){
        r<-r * pbeta(x,y[j]+1,n[j]-y[j]+1)
      }
      return(r) }
    
    ans[i] = integrate(f,0,1)$value
  }
  return(ans)
}

out<-c()
k<-1
for(i in 1:length(a)){
  out<-c(out,compute.probopt(c(a[i],b[i]),c(t1[i],t2[i])))
  cat(compute.probopt(c(a[i],b[i]),c(t1[i],t2[i])),"\n")
  cat(out[k+1]/(out[k]+out[k+1]),"\n")
  k<-k+2
}

##Same as earlier funcion
##compute each individually with other formula




##################
##MCMC - Markov Chain Monte Carlo Sampling
sim.post<-function(y,n,ndraws){
  k<-length(y)
  ans<-matrix(nrow=ndraws,ncol=k)
  no<-n-y
  for(i in 1:k){
    ans[,i]<-rbeta(ndraws,y[i]+1,no[i]+1)
  }
  return(ans)
}
prob.winner<-function(post){
  k<-ncol(post)
  w<-table(factor(max.col(post),levels=1:k))
  return(w/sum(w))
}
compute.win.prob<-function(y,n,ndraws){
  return(prob.winner(sim.post(y,n,ndraws)))
}

for(i in 1:length(a)){
  cat(compute.win.prob(c(a[i],b[i]),c(t1[i],t2[i]),10000),"\n")
}


a1<-c(2315,2169,2192)
n1<-c(25252,25265,25312)

compute.win.prob(a1,n1,10000)


a1<-c(528,518)
n1<-c(4182,4195)

compute.win.prob(a1,n1,10000)






