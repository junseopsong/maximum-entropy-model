load P1.txt; % independent MEM
load P2.txt; % pairwise MEM
load PN.txt; % empirical probability

M = 2^9;

D1 = 0; D2 = 0; S1 = 0; S2 = 0; SN = 0;
for i = 1:M
    if PN(i) > 0
        if P1(i) > 0
            D1 = D1 + PN(i)*log2(PN(i)/P1(i));
        end
        if P2(i) > 0
            D2 = D2 + PN(i)*log2(PN(i)/P2(i));
        end
        SN = SN - PN(i)*log(PN(i));
    end
    if P1(i) > 0
        S1 = S1 - P1(i)*log(P1(i));
    end
    if P2(i) > 0
        S2 = S2 - P2(i)*log(P2(i));
    end
end

RD = (D1-D2)/D1
RS = (S1-S2)/(S1-SN);
Reliability = RS/RD

figure;
scatter(PN, P2, 10, 'k', 'filled'); axis square;
hold on; plot(10e-6:10e-6:1,10e-6:10e-6:1,'k');
set(gca,'FontName','Arial','FontSize',12);
set(gca,'XScale','log','YScale','log'); axis([10e-6 1 10e-6 1]);
