function TF = get_TF(k,rmat,obj)
            TF=exp(-1i*k*rmat.^2/(2*obj.focal));
end