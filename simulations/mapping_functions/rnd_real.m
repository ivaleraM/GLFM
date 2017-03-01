function x = rnd_real(Zn,Bd,numSamples, s2y, s2u)
    % function to get random samples from the distribution
    x = sqrt(s2u+s2y) .* randn(1,numSamples) + Zn * Bd;
end