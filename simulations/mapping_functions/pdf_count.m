function pdf = pdf_count(X,Zn,Bd,w,s2y)

pdf = normcdf(fpos_1(X+1,w), Zn*Bd, sqrt(s2y)) - ...
    normcdf(fpos_1(X,w), Zn*Bd, sqrt(s2y));