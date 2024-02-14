# getCovStats works

    Code
      getCovStats(data, covariates, missVal = -99)
    Output
      $WT
         5%   95% 
       58.5 126.0 
      
      $HT
       5% 95% 
      152 185 
      
      $SEX
      [1] 1 2
      

---

    Code
      getCovStats(data, covariates, missVal = -99, probs = c(0.1, 0.9))
    Output
      $WT
        10%   90% 
       64.3 116.0 
      
      $HT
      10% 90% 
      156 183 
      
      $SEX
      [1] 1 2
      

---

    Code
      getCovStats(data, covariates, missVal = -99)
    Output
      $WT
         5%   95% 
       58.5 126.0 
      
      $HT
       5% 95% 
      152 185 
      
      $SEX
      [1] 1 2
      
      $FakeWT
          5%    95% 
      0.0585 0.1260 
      

---

    Code
      getCovStats(data, covariates, missVal = -99, nsig = 2)
    Output
      $WT
       5% 95% 
       59 130 
      
      $HT
       5% 95% 
      150 190 
      
      $SEX
      [1] 1 2
      
      $FakeWT
         5%   95% 
      0.059 0.130 
      

