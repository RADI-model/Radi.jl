year = 5;

figure(1)
clf

subplot(331)
plot([0:0.01:max(max(dO2f))],[0:0.01:max(max(dO2f))])
hold on
scatter(dO2f(:,year),dO2(:,year))
title('O2')

subplot(332)
plot([0:0.01:max(max(dtCO2f))],[0:0.01:max(max(dtCO2f))])
hold on
scatter(dtCO2f(:,year),dtCO2(:,year))
title('tCO2')

subplot(333)
plot([0:0.01:max(max(dtNO3f))],[0:0.01:max(max(dtNO3f))])
hold on
scatter(dtNO3f(:,year),dtNO3(:,year))
title('tNO3')

subplot(334)
plot([29.8:0.01:max(max(dtSO4f))],[29.8:0.01:max(max(dtSO4f))])
hold on
scatter(dtSO4f(:,year),dtSO4(:,year))
title('tSO4')

subplot(335)
plot([0:0.01:max(max(dtPO4f))],[0:0.01:max(max(dtPO4f))])
hold on
scatter(dtPO4f(:,year),dtPO4(:,year))
title('tPO4')

subplot(336)
plot([0:0.01:max(max(dtNH4f))],[0:0.01:max(max(dtNH4f))])
hold on
scatter(dtNH4f(:,year),dtNH4(:,year))
title('tNH4')

subplot(337)
plot([0:0.01:max(max(dtH2Sf))],[0:0.01:max(max(dtH2Sf))])
hold on
scatter(dtH2Sf(:,year),dtH2S(:,year))
title('tH2S')

subplot(338)
plot([0:1e-5:max(max(dFef))],[0:1e-5:max(max(dFef))])
hold on
scatter(dFef(:,year),dFeII(:,year))
title('Fe')

subplot(339)
plot([0:1e-7:max(max(dMnf(:,year)))],[0:1e-7:max(max(dMnf(:,year)))])
hold on
scatter(dMnf(:,year),dMnII(:,year))
title('Mn')

figure(2)
clf

subplot(331)
plot([10:0.01:max(max(dCaf))],[10:0.01:max(max(dCaf))])
hold on
scatter(dCaf(:,year),dCa(:,year))
title('Ca')

subplot(332)
plot([2:0.01:max(max(dtalkf))],[2:0.01:max(max(dtalkf))])
hold on
scatter(dtalkf(:,year),dalk(:,year))
title('talk')

subplot(333)
plot([200:0.01:max(max(pfocf))],[200:0.01:max(max(pfocf))])
hold on
scatter(pfocf(:,year),pfoc(:,year))
title('pfoc')

subplot(334)
plot([2900:0.01:max(max(psocf))],[2900:0.01:max(max(psocf))])
hold on
scatter(psocf(:,year),psoc(:,year))
title('psoc')

subplot(335)
plot([2:0.01:max(max(pFeOH3f))],[2:0.01:max(max(pFeOH3f))])
hold on
scatter(pFeOH3f(:,year),pFeOH3(:,year))
title('pFeOH3')

subplot(336)
plot([2:0.01:max(max(pMnO2f))],[2:0.01:max(max(pMnO2f))])
hold on
scatter(pMnO2f(:,year),pMnO2(:,year))
title('pMnO2')

subplot(337)
plot([2:0.01:max(max(pcalcitef))],[2:0.01:max(max(pcalcitef))])
hold on
scatter(pcalcitef(:,year),pcalcite(:,year))
title('pcalcite')

subplot(338)
plot([2:0.01:max(max(paragonitef))],[2:0.01:max(max(paragonitef))])
hold on
scatter(paragonitef(:,year),paragonite(:,year))
title('paragonite')

subplot(339)
plot([29000:0.01:max(max(procf))],[29000:0.01:max(max(procf))])
hold on
scatter(procf(:,year),proc(:,year))
title('proc')