library(nimble)

# Function to standardize covariates
standardize <- nimbleFunction(
  run= function(covariate = double(1)){
    mean_cov = mean(covariate)
    sd_cov = sd(covariate)
    new_cov = (covariate - mean_cov)/sd_cov
    returnType(double(1))
    return(new_cov)
  }
)

Cstandardize = compileNimble(standardize)

# Similarity Function
similarity1 = nimbleFunction(
  run = function(data = double(1), mu = double(0, default = 0), m0 = double(0, default = 0), v= double(0, default = 1), 
                 v0 = double(0, default = 1)){
    
    data1 = (data) - mean(data)
    m_star <- (1/v * sum(data1) + 1/v0*m0)/(length(data1)/v + 1/v0)
    v_star <- 1/(length(data1)/v + 1/v0)
    post_ans = dnorm(mu, m_star, sqrt(v_star))
    
    lk_ans <- dnorm(data1, mu, sqrt(v))
    lk_ans1 = prod(lk_ans)
    
    prior_ans <- dnorm(mu, m0, sqrt(v0))
    
    ans1 = prior_ans*lk_ans1/post_ans
    
    # Get without the last one
    data2 = data[1:(length(data)-1)]
    data2 = data2 - mean(data2)
    m_star2 <- (1/v * sum(data2) + 1/v0*m0)/(length(data2)/v + 1/v0)
    v_star2 <- 1/(length(data2)/v + 1/v0)
    post_ans2 = dnorm(mu, m_star2, sqrt(v_star2))
    
    lk_ans <- dnorm(data2, mu, sqrt(v))
    lk_ans2 = prod(lk_ans)
    
    prior_ans2 <- dnorm(mu, m0, sqrt(v0))
    
    ans2 = prior_ans2*lk_ans2/post_ans2
    
    
    # This may be an issue, come back to this if it aborts again
    if(ans2 == 0){
      return(0)
    }
    
    ans3 = ans1/ans2
    
    
    returnType(double(0))
    return(ans3)
  }
)


ppmxPrediction <- nimbleFunction(
  run = function(data = double(2), new = double(1), 
                 clusters = double(1), numClusters = double(0), nobs_in_clus = double(1),
                 v =  double(0), v0 = double(0), mu = double(0), m0 = double(0), 
                 y = double(1))
  {
    # Number of covariates. The number of rows in the matrix
    n_iters = length(data[1,])
    # The number of observations. "n". Number of columns in the matrix
    n_obs = length(data[,1])
    # Number of responses This should be the same as n_iters. The length of the "y" vector.
    n_covariates = length(new)
    
    #print(data)
    # Make sure that the two are the same length
    if(n_iters != n_covariates){
      #print("number of covariates in new and data differ")
    }
    
    # For each covariate
    for(i in 1:n_iters){
      #Standardize the covariate values with the corresponding new vector observation.
      vec = standardize(c(data[,i],new[i]))
      # Separate the new value back out
      new[i] = vec[n_obs+1]
      # Separate the covariate values back out
      data[,i] = vec[1:(n_obs)]
    }
    
    # Make a matrix with n(col) = number of covariates and n(row) = number of Clusters
    sim_matrix = nimMatrix(value = NA, nrow = numClusters, ncol = n_iters)
    #print(data)
    

    # For each Cluster i
    for(i in 1:numClusters){
      # For each covariate j
      for(j in 1:n_iters){
        # Populate the similarity matrix at i,j with the similarity function evaluated at the covariate values that belong to cluster i
        # and the jth new value
        sim_matrix[i,j] = similarity1(c(data[ ,j][clusters == i], new[j]), mu, m0, v, v0)
        #print(sim_matrix[i,j])
      }
    }
    #print(sim_matrix)
    # It looks like the similarity matrix is generated how I hoped.
    
    # Make a vector that will hold the product of the cluster scores.
    sim_product = nimNumeric(length = numClusters)
    # For each covariate k
    for(k in 1:numClusters){
      # Populate sim_product vector at the kth spot with the product of similarity scores for the cluster.
      # Multiplies across the row. 
      # Answers the question: How much weight should this covariate have?
      sim_product[k] = (nobs_in_clus[k]*prod(sim_matrix[k,1:n_iters]))
    }
    sim_product <- sim_product/sum(sim_product)
    # Make a prediction vector that holds numClusters objects
    pred_vector = nimNumeric(length = numClusters)
    # For each cluster
    for(l in 1:numClusters){
      # Populate the prediction with the weight of the cluster multiplied by the average of the cluster's y values
      pred_vector[l] = sim_product[l]*mean(y[clusters == l])
    }
    # Sum the products in l
    sum_preds = sum(pred_vector)
    # Sum the similarities
    sum_sims = sum(sim_product)
    # Divide the summed predictions by the summed similarities to get back on original scale
#    return(sum_preds/sum_sims)
    return(sum_preds)
    returnType(double(0))
  }
)

ppmxPrediction(matrix(c(1,1.1,2.1,2,3,3.1,4,4.1,5,5.1,6,6.1,7,7.1,8,8.1,9,9.1,10,10.1,11,11.1,12,12.1), nrow = 8, ncol = 3),
             new = c(4,8,12), 
             clusters = c(1,2,2,2,3,3,3,4), numClusters = 4, nobs_in_clus=c(1,3,3,1),
             v = 0.1, v0 = 1, mu = 0, m0 = 0, y = c(1, 10,10,11, 20,21,24,80))

CppmxPrediction <- compileNimble(ppmxPrediction)

CppmxPrediction(matrix(c(1,1.1,2.1,2,3,3.1,4,4.1,5,5.1,6,6.1,7,7.1,8,8.1,9,9.1,10,10.1,11,11.1,12,12.1), nrow = 8, ncol = 3),
               new = c(4,8,12), 
               clusters = c(1,2,2,2,3,3,3,4), numClusters = 4, nobs_in_clus=c(1,3,3,1),
               v = 0.1, v0 = 1, mu = 0, m0 = 0, y = c(1, 10,10,11, 20,21,24,80))
