% data: dataset structure
% data = 
%  struct with fields:
%
%    continuous: [5×392 double]
%      discrete: [3×392 double]
%         names: {'mpg'  'displacement'  'horsepower'  ...}
%        nClass: [3×1 double]
% Dz: number of latent features
function res = wrapper_mixed_FA(data,Dz)
  %MIXED-DATA FA 
  setSeed(3); % reproducibility
  % 1 of M encoding
  data.categorical = encodeDataOneOfM(data.discrete, data.nClass);
  % initialize
  opt=struct('Dz', Dz, 'nClass', data.nClass);
  [params0, data] = initMixedDataFA(data, [], opt);
  params0.a = 1;
  params0.b = 1;
  % learn theta with EM algorithm
  options = struct('maxNumOfItersLearn', 200,...
                    'lowerBoundTol', 1e-6,...
                    'estimateBeta',1,...% estimate loading factos
                    'estimateMean', 1,...% estimate prior mean (which is equivalent to estimating bias)
                    'estimateCovMat',0);
  funcName = struct('inferFunc', @inferMixedDataFA_miss, 'maxParamsFunc', @maxParamsMixedDataFA);
  [params, logLik] = learnEm(data, funcName, params0, options);
  % Obtain p(z|y,\theta)
  params.psi = randn(size(data.categorical));% initialize variational parameters
  [res.ss, res.logLik, res.postDist] = inferMixedDataFA_miss(data, params, []);

  res.X_pred_mixedFA = complete_matrix(params,res.postDist)
