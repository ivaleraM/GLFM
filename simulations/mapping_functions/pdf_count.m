function pdf = pdf_count(X,Zn,Bd,w,s2y, fpos_1_handler)

pdf = normcdf(fpos_1_handler(X+1,w), Zn*Bd, sqrt(s2y)) - ...
    normcdf(fpos_1_handler(X,w), Zn*Bd, sqrt(s2y));