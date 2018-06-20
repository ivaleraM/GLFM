disp('Script to check call to GLFMR library')
addpath( genpath('../../src/') );
D=5;
N = 3;
disp(['The data matrix is: \n']) 
data.X= [1.0, 1, -0.3, 1, 1; 6.3, 2, 3.8, 23, 1; 11, 3, 4.1, 4, 2]
data.C=['p','o','G','N','c'];
disp(['The hidden matrix Z is:\n'])
hidden.Z= [1.0 , 0 ; 1 , 1; 1 , 1]
hidden = GLFM_infer(data, hidden);
disp(['The inferred matrix Z is:\n'])
hidden.Z
disp(['The inferred matrix B is:\n'])
hidden.B
disp('Succesful')