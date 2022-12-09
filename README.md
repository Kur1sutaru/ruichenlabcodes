## Some tips to run cellqc pipeline


* Create simbolic links to the files
  ln -s
  
* unzip all the .tar.gz files, and create a folder to all of them
  analysis.tar.gz   filtered_feature_bc_matrix.tar.gz  raw_feature_bc_matrix.tar.gz

* Optional configuration for dropkick [to be tested] 
  dropkick:
  skip: false
  method: multiotsu
  numthreads: 2  #INSTEAD 1
*
