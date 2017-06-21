#var1<-c("HOLA")
# sapply(var1, tolower)

GLFM_infer<-function(data,varargin){
  x <- "red"
  > switch(x, red="cloth", size=5, name="table")
  varargin_size<-length(varargin)
  switch(varargin_size,"0"=hidden <-c(),"1"="bu","2"="bu1")
  
  switch length(varargin)
  case 0
  hidden = [];
  % initialize default values for params
  params = init_default_params(data, []);
  
  case 1, hidden = varargin{1};
  params = init_default_params(data, []);
  
  case 2, hidden = varargin{1}; params = varargin{2};
  params = init_default_params(data, params); % eventually complete params structure
  
  otherwise
  error('Incorrect number of input parameters: should be 1, 2 or 3');
  end
  
}