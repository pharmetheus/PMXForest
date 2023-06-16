# getCovStats works

    Code
      getCovStats(data, covariates, missVal = -99)
    Output
      $WT
           5%     95% 
       58.514 126.000 
      
      $HT
         5%   95% 
      152.4 185.4 
      
      $SEX
      [1] 1 2
      

---

    Code
      getCovStats(data, covariates, missVal = -99, probs = c(0.1, 0.9))
    Output
      $WT
         10%    90% 
       64.29 116.10 
      
      $HT
        10%   90% 
      156.0 182.8 
      
      $SEX
      [1] 1 2
      

