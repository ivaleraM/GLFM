# Demo toyImages
source("GLFM_infer.R")
source("GLFM_computeMAP.R")
source("generate_toy_images.R")

# Generative model
N<-1000
s2x<-0.5

#Initialisation
data_gen<- generate_toy_images(N,s2x)

# Inference

# Plot