function h=ploterrs(xin,yin,err)
%plots the error bar vec and fit for the err returned to the base workspace
%by matrix_expfit(.....,'errpoke')


h=errorbar(err.taxis,squeeze(err.data(yin,xin,:)),squeeze(err.noise(yin,xin,:)),'b');
hold on;
plot(err.taxis,squeeze(err.fit(yin,xin,:)),'k');