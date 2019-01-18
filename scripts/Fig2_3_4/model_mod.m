function dy=model_mod(t,y,p)

dy=zeros(18,1);

% The parameters 25-32 describe whether a specific regulator type is
% present, they can be either 0 (not present) or 1 (present)
if p(25) == 1
    cXA1 = y(3)^p(9)/(y(3)^p(9) + p(10)^p(9));
    cXA2 = y(4)^p(9)/(y(4)^p(9) + p(10)^p(9));
else
    cXA1=1;
    cXA2=1;
end

if p(26) == 1
    tXA = (0.5*(y(5)+y(6)))^p(11)/((0.5*(y(5)+y(6)))^p(11)+p(12)^p(11));
elseif p(26) ==2 & p(28)~=1
     tXA = (0.5*(p(15)))^p(11)/((0.5*(p(15)))^p(11)+p(12)^p(11));
else
    tXA = 1;
end

if p(27) == 1
    cXR1 = (1 - y(7)^p(13) / (y(7)^p(13) + p(14)^p(13)));
    cXR2 = (1 - y(8)^p(13) / (y(8)^p(13) + p(14)^p(13)));
else
    cXR1 = 1;
    cXR2 = 1;
end

if p(28) == 1
    tXR = (1 - (0.5*(y(9)+y(10)))^p(15) / ((0.5*(y(9)+y(10)))^p(15) + p(16)^p(15)));
else
    tXR = 1;
end

if p(29) == 1
    ecXA1 = y(11)^p(17)/(y(11)^p(17) + p(18)^p(17));
    ecXA2 = y(12)^p(17)/(y(12)^p(17) + p(18)^p(17));
else
    ecXA1=1;
    ecXA2=1;
end


if p(30) == 1
    etXA = (0.5*(y(13)+y(14)))^p(19) / ((0.5*(y(13)+y(14)))^p(19) + p(20)^p(19));
else
    etXA = 1;
end

if p(31) == 1
    ecXR1 = (1 - y(15)^p(21) / (y(15)^p(21) + p(22)^p(21)));
    ecXR2 = (1 - y(16)^p(21) / (y(16)^p(21) + p(22)^p(21)));
else
    ecXR1 = 1;
    ecXR2 = 1;
end

if p(32) == 1
    etXR = (1 - (0.5*(y(17)+y(18)))^p(23) / ((0.5*(y(17)+y(18)))^p(23) + p(24)^p(23)));
else
    etXR = 1;
end

% equations for y1 and y2 (=Xist transcribed from the two chromosomes
% p33 is set to 0 for simulating male cells

dy(1) = cXA1 * tXA * cXR1 * tXR * etXA * etXR *ecXA1 *ecXR1 - y(1);
dy(2) = p(33) * cXA2 * tXA * cXR2 * tXR * etXA * etXR *ecXA2 *ecXR2 - y(2);

% equations for cXA
    dy(3) = (1 - y(1)^p(1) / (y(1)^p(1) + p(2)^p(1))) - y(3);
    
    dy(4) = p(33) * (1 - y(2)^p(1) / (y(2)^p(1) + p(2)^p(1))) - y(4);
% equations for tXA
    dy(5) = (1 - y(1)^p(3) / (y(1)^p(3) + p(4)^p(3))) - y(5);

    dy(6) = p(33) * (1 - y(2)^p(3) / (y(2)^p(3) + p(4)^p(3))) - y(6);

% equations for cXR
    dy(7) = (1 - y(1)^p(5) / (y(1)^p(5) + p(6)^p(5))) - y(7);

    dy(8) = p(33) * (1 - y(2)^p(5) / (y(2)^p(5) + p(6)^p(5))) - y(8);

% equations for tXR
    dy(9) = (1 - y(1)^p(7) / (y(1)^p(7) + p(8)^p(7))) - y(9);

    dy(10) = p(33) * (1 - y(2)^p(7) / (y(2)^p(7) + p(8)^p(7))) - y(10);

% equations for ecXA
    dy(11) = 1 - y(11);

    dy(12) = p(33) - y(12);

% equations for etXA
    dy(13) = 1 - y(13);

    dy(14) = p(33) - y(14);
    
% equations for ecXR
    dy(15) = 1 - y(15);

    dy(16) = p(33) - y(16);

% equations for etXR
    dy(17) = 1 - y(17);

    dy(18) = p(33) - y(18);

end
