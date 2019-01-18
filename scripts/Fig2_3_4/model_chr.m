function dy=model_chr(t,y,p)

dy=zeros(18,1);

% Xist
% regulation by cXA

    tXA = (0.5*(p(15)))^p(11)/((0.5*(p(15)))^p(11)+p(12)^p(11));

    cXR1 = (1 - y(7)^p(13) / (y(7)^p(13) + p(14)^p(13)));

dy(1) = tXA * cXR1  - y(1);



% cXR
    dy(7) = (1 - y(1)^p(5) / (y(1)^p(5) + p(6)^p(5))) - y(7);

end
