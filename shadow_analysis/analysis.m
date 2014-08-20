function analysis
%     [r1, g1, b1, r2, g2, b2] = textread('rgb.txt', '%d%d%d%d%d%d');
%     rr1 = r1./b1; gg1 = g1./b1;
%     rr2 = r2./b2; gg2 = g2./b2;
%     plot(rr2-rr1, gg2-gg1);
    
% 1:shadow 2:unshadow

    [l1, a1, b1, l2, a2, b2] = textread('lab.txt', '%d%d%d%d%d%d');
    %plot(a2-a1, b2-b1);
    
    drawHist((b2-b1)./(a2-a1));
    
end

function drawHist(x)
    [a,b]=hist(x,100); 
    a=a/length(x); 
    bar(b,a);
end