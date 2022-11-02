mydata2=ToMR(mydata) #Transfer to MR

#File for MRfunction
df_mr=df_data(mydata2) 
ldrho=mydata2[[1]]@correlation 

ivw=MendelianRandomization::mr_ivw(mydata2[[1]], 
                                   robust = T, 
                                   penalized = T, 
                                   correl= T) #IVW Screening
ivw 

res=run_MR_all(df_mr,ldrho) 
